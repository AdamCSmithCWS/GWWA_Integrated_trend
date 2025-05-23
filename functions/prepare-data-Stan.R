#' Wrangle data to use for modelling input
#'
#' \code{prepare_data} subsets raw BBS data by selected species and
#'    and wrangles stratified data for use as input for models.
#'
#' @param strat_data Large list of stratified data returned by \code{stratify()}
#' @param species_to_run Character string of the English name of the species to run
#' @param model Character string of model to be used.
#'   Options are "slope", "firstdiff", "gam", "gamye.
#' @param heavy_tailed Logical indicating whether the extra-Poisson error distribution should be modeled as a t-distribution, with heavier tails than the standard normal distribution. Default is currently FALSE, but recent results suggest users should strongly consider setting this to TRUE, even though it requires much longer convergence times
#' @param n_knots Number of knots to be used in GAM function
#' @param basis Which version of the basis-function to use for the GAM smooth, the default is "original" the same basis used in Smith and Edwards 2020 and "mgcv" is an alternate that uses the "tp" basis from the packages mgcv (also used in brms, and rstanarm). If using the "mgcv" option, the user may want to consider adjusting the prior distributions for the parameters and their precision
#' @param min_year Minimum year to keep in analysis
#' @param max_year Maximum year to keep in analysis
#' @param min_n_routes Minimum routes per strata where species has been observed.
#'   Defaults to 3
#' @param min_max_route_years Minimum number of years with non-zero observations
#'   of species on at least 1 route. Defaults to 3
#' @param min_mean_route_years Minimum average of years per route with the
#'   species observed. Defaults to 1.
#' @param strata_rem Strata to remove from analysis. Defaults to NULL
#' @param quiet Should progress bars be suppressed?
#' @param sampler Which MCMC sampling software to use. Currently bbsBayes only
#'   supports "jags".
#' @param ... Additional arguments
#'
#' @return List of data to be used for modelling, including:
#'   \item{model}{The model to be used}
#'   \item{heavy_tailed}{Logical indicating whether the extra-Poisson error distribution should be modeled as a t-distribution}
#'   \item{min_nu}{if heavy_tailed is TRUE, minimum value for truncated gamma on DF of t-distribution noise default is 0 and user must change manually after function is run}
#'   \item{ncounts}{The number of counts containing useful data for the species}
#'   \item{nstrata}{The number of strata used in the analysis}
#'   \item{ymin}{Minimum year used}
#'   \item{ymax}{Maximum year used}
#'   \item{nonzeroweight}{Proportion of routes in each strata with species obervation}
#'   \item{count}{Vector of counts for the species}
#'   \item{strat}{Vector of strata to be used in the analysis}
#'   \item{obser}{Vector of unique observer-route pairings}
#'   \item{year}{Vector of years for each count}
#'   \item{firstyr}{Vector of indicator variables as to whether an observer was a first year}
#'   \item{month}{vector of numeric month of observation}
#'   \item{day}{vector of numeric day of observation}
#'   \item{nobservers}{Total number of observer-route pairings}
#'   \item{fixedyear}{Median of all years (ymin:ymax), included only with slope and firstdiff models}
#'   \item{nknots}{Number of knots to use for smooting functions, included only with GAM}
#'   \item{X.basis}{Basis function for n smoothing functions, included only with GAM}
#'
#' @importFrom stats median
#' @importFrom progress progress_bar
#' @importFrom mgcv s
#' @importFrom mgcv smoothCon
#' @export
#'
#' @examples
#' # Toy example with Pacific Wren sample data
#' # First, stratify the sample data
#'
#' strat_data <- stratify(by = "bbs_cws", sample_data = TRUE)
#'
#' # Prepare the stratified data for use in a model. In this
#' #   toy example, we will set the minimum year as 2009 and
#' #   maximum year as 2018, effectively only setting up to
#' #   model 10 years of data. We will use the "first difference
#' #   model.
#' model_data <- prepare_data(strat_data = strat_data,
#'                            species_to_run = "Pacific Wren",
#'                            model = "firstdiff",
#'                            min_year = 2009,
#'                            max_year = 2018)
#'
#' # You can also specify the GAM model, with an optional number of
#' # knots to use for the GAM basis.
#' # By default, the number of knots will be equal to the floor
#' # of the total unique years for the species / 4
#' model_data <- prepare_data(strat_data = strat_data,
#'                            species_to_run = "Pacific Wren",
#'                            model = "gam",
#'                            n_knots = 9)
#'
#'

prepare_data <- function(strat_data = NULL,
                         species_to_run = NULL,
                         model = NULL,
                         heavy_tailed = TRUE,
                         n_knots = NULL,
                         min_year = NULL,
                         max_year = NULL,
                         min_n_routes = 3,
                         min_max_route_years = 3,
                         min_mean_route_years = 1,
                         strata_rem = NULL,
                         quiet = FALSE,
                         sampler = "Stan",
                         basis = "mgcv",
                         use_pois = FALSE,
                         calculate_nu = FALSE,
                         calculate_log_lik = FALSE,     
                         calculate_CV = FALSE,
                         ...)
{
  x <- NULL; rm(x)
  if (is.null(strat_data))
  {
    stop("No data supplied to prepare_data()."); return(NULL)
  }
  if (is.null(species_to_run))
  {
    stop("No species specified."); return(NULL)
  }
  if (is.null(model))
  {
    stop("No model specified."); return(NULL)
  }
  if (isFALSE(is.element(model, c("slope", "firstdiff", "gam", "gamye"))))
  {
    stop("Invalid model specified"); return(NULL)
  }
  if (isFALSE(heavy_tailed) & isFALSE(use_pois))
  {
    stop("Heavy-tailed models are implied with Negative Binomial (use_pois == FALSE)"); return(NULL)
  }
  

  birds <- strat_data$bird_strat
  route <- strat_data$route_strat
  species <- strat_data$species_strat


  ##################################################################
  # Previous function from Bugs Data Prep
  ##################################################################
  sp_eng <- species_to_run
  sp_aou <- get_species_aou(species, species_to_run)

  spsp <- birds[which(birds$AOU == sp_aou),] # Gets all observations of the species of interest
  names(spsp)[which(names(spsp) == "SpeciesTotal")] <- "TotalInd" # change header to totalInd

  spsp.c <- merge(route,spsp[,-which(names(spsp) %in% c("countrynum",
                                                        "rt.uni.y",
                                                        "rt.uni",
                                                        "RPID"))],
                  by = c("statenum",
                         "Route",
                         "Year",
                         "BCR",
                         "RouteDataID"),all.x = TRUE)

  spsp.c[which(is.na(spsp.c$TotalInd)),"TotalInd"] <- 0

  if (!is.null(strata_rem))
  {
    spsp.c <- spsp.c[-which(spsp.c$strat_name %in% strata_rem),]
  }

  if (!is.null(min_year))
  {
    spsp.c <- spsp.c[which(spsp.c$Year >= min_year), ]
  }

  if (!is.null(max_year))
  {
    spsp.c <- spsp.c[which(spsp.c$Year <= max_year), ]
  }

  if (!isTRUE(quiet))
  {
    message("Preparing data")
    pb <- progress::progress_bar$new(
      format = "\r[:bar] :percent eta: :eta",
      clear = FALSE,
      total = length(unique(spsp.c$strat_name)) + 26,
      width = 80)
    pb$tick(0)
  }

  spsp_routes_ever <- unique(spsp.c[which(spsp.c$TotalInd != 0),c("strat_name","rt.uni")]) #routes which have had at least 1 species counted
  spsp_routes_never <- unique(spsp.c[-which(spsp.c$rt.uni %in% spsp_routes_ever$rt.uni),c("strat_name","rt.uni")]) #routes that have not had this species before
  spsp.c2 <- spsp.c[which(spsp.c$rt.uni %in% spsp_routes_ever$rt.uni),] #All data counts for routes which has had seen
  if (!isTRUE(quiet)){pb$tick()}

  # first year each given route was run
  miny.rt <- tapply(spsp.c2[,"Year"],spsp.c2[,"rt.uni"],min)
  miny.df <- data.frame(rt.uni = names(miny.rt),fyr_rt_run = as.integer(miny.rt))
  if (!isTRUE(quiet)){pb$tick()}

  # number of times the species was seen on the given route since it has run
  n.yr.rt <- tapply(spsp.c[which(spsp.c$TotalInd != 0),"Year"],spsp.c[which(spsp.c$TotalInd != 0),"rt.uni"],length)
  n.yr.df <- data.frame(rt.uni = names(n.yr.rt),nyr_rt_run = as.integer(n.yr.rt))
  if (!isTRUE(quiet)){pb$tick()}

  if(nrow(spsp_routes_ever) > 0)
  {
    spsp_routes_ever <- merge(spsp_routes_ever,miny.df,by = "rt.uni")
    spsp_routes_ever <- merge(spsp_routes_ever,n.yr.df,by = "rt.uni")
    if (!isTRUE(quiet)){pb$tick()}

    # this will give some ever/never stats PER STRATUM
    pR <- data.frame(strat = unique(spsp.c$strat_name), nr.ever = NA, nr.never = NA, p.r.ever = NA,nry.ever = NA,meanry.ever = NA)
    for(p in unique(spsp.c$strat_name))
    {
      if (!isTRUE(quiet)){pb$tick()}
      # routes ever observed in this strata
      pR[pR$strat == p,"nr.ever"] <- nrow(spsp_routes_ever[spsp_routes_ever$strat_name == p,])
      # routes never observed in this strata
      pR[pR$strat == p,"nr.never"] <- nrow(spsp_routes_never[spsp_routes_never$strat_name == p,])

      # proportion of routes that have ever observed by stratum
      pR[pR$strat == p,"p.r.ever"] <-  pR[pR$strat == p,"nr.ever"]/(pR[pR$strat == p,"nr.ever"]+ pR[pR$strat == p,"nr.never"])

      # total counts that have included this species (ex: 5 routes in year 1, 4 in year 2 = 9 overall)
      pR[pR$strat == p,"nry.ever"] <-  length(spsp.c[which(spsp.c$strat_name == p & spsp.c$TotalInd > 0),"rt.uni.y"])

      pR[pR$strat == p,"fy.wspecies"] <-  suppressWarnings(min(spsp_routes_ever[spsp_routes_ever$strat_name == p,"fyr_rt_run"]))
      pR[pR$strat == p,"max.nry"] <-  suppressWarnings(max(spsp_routes_ever[spsp_routes_ever$strat_name == p,"nyr_rt_run"]))

      pR[pR$strat == p,"meanry.ever"] <-  length(spsp.c[which(spsp.c$strat_name == p & spsp.c$TotalInd > 0),"rt.uni.y"])/pR[pR$strat == p,"nr.ever"]
    }

    pR[,"strat"] <- as.character(pR[,"strat"])

    #gets rid of infinite values
    pR[which(pR$fy.wspecies > 2100),"fy.wspecies"] <- NA
    pR[which(pR$max.nry < 0),"max.nry"] <- NA
    if (!isTRUE(quiet)){pb$tick()}

    spsp.c <- spsp.c[which(spsp.c$rt.uni %in% unique(spsp_routes_ever$rt.uni)),]
    #spsp.c.drop <- spsp.c[-which(spsp.c$rt.uni %in% unique(spsp_routes_ever$rt.uni)),]
    if (!isTRUE(quiet)){pb$tick()}

    
    spsp.2 <- spsp.c
    spsp.2 <- spsp.2[which(spsp.2$strat_name %in% pR[which(pR$nr.ever >= min_n_routes & pR$max.nry >= min_max_route_years & pR$meanry.ever >= min_mean_route_years),"strat"]),]
    if (!isTRUE(quiet)){pb$tick()}

    rts.used <- unique(spsp.2[,c("rt.uni","strat_name")])

    #get number of routes used per strata
    #incidentally this ends up being all the strata that are used
    rts.summary <- tapply(rts.used$rt.uni,rts.used$strat_name,length)
    nrts.used <- data.frame(strat_name = names(rts.summary),nroutes.used = as.integer(rts.summary))
    if (!isTRUE(quiet)){pb$tick()}

    spsp.temp.2 <- merge(spsp.2,pR[,c("strat","p.r.ever","meanry.ever","max.nry")], by.x = "strat_name", by.y = "strat",all.x = TRUE)
    if (!isTRUE(quiet)){pb$tick()}

    spsp.2 <- spsp.temp.2

    spsp.2[,"yr"] <- (spsp.2[,"Year"]+1)-min(spsp.2[,"Year"])

    spsp.2[,"count"] <- spsp.2[,c("TotalInd")]

    spsp_f <- spsp.2
    if (!isTRUE(quiet)){pb$tick()}
  }else
  {
    for(p in length(unique(spsp.c$strat_name)) + 6)
    {
      if (!isTRUE(quiet)){pb$tick()}
    }
    spsp_f <- spsp.c[-c(1:nrow(spsp.c)),]
    if (!isTRUE(quiet)){pb$tick()}
  }

  spsp_f[,"stratcode"] <- spsp_f[,"strat_name"]
  spsp_f[,"stratum"] <- as.numeric(factor(spsp_f$strat_name)); if (!isTRUE(quiet)){pb$tick()}

  spsp_f <- spsp_f[order(spsp_f$stratum,spsp_f$Route,spsp_f$Year),]; if (!isTRUE(quiet)){pb$tick()}

  strat.list <- unique(spsp_f[,c("statenum","countrynum","BCR","stratum","strat_name","stratcode","State")])
  strat.use <- strat.list[,"stratcode"]

  pR2 <- merge(pR[which(pR$strat %in% strat.list[,"strat_name"]),],strat.list,by.x = "strat", by.y = "strat_name")
  if (!isTRUE(quiet)){pb$tick()}

  pR2 <- pR2[order(pR2$stratum),]
  if (!isTRUE(quiet)){pb$tick()}

  spsp_f <- spsp_f[order(spsp_f$stratum,spsp_f$Route,spsp_f$Year),]
  if (!isTRUE(quiet)){pb$tick()}

  spsp_f[,"rt.uni.ob"] <- paste(spsp_f$statenum,spsp_f$Route,spsp_f$ObsN,sep = "-")

  strata.adj <-  pR2[,"stratcode"]

  pR.wts <- unique(spsp_f[,c("strat_name","p.r.ever")])

  tmp <- unique(spsp_f[,c("stratum","rt.uni.ob")])

  # number of observers per stratum
  nobservers <- tapply(tmp$rt.uni.ob,tmp$stratum,length)
  if (!isTRUE(quiet)){pb$tick()}

  for (s in 1:length(nobservers)) {
    sel1 <- which(spsp_f$stratum == s)

    tmp <- unique(spsp_f[sel1,c("stratum","rt.uni.ob")])

    tmp[,"obser"] <- as.numeric(factor(tmp[,"rt.uni.ob"]))
    if(s == 1)
    {
      obsfact <- tmp
    }else
    {
      obsfact <- rbind(obsfact,tmp)
    }
  }
  if (!isTRUE(quiet)){pb$tick()}

  spsp_ft <- merge(spsp_f,obsfact,by = c("stratum","rt.uni.ob"))
  if (!isTRUE(quiet)){pb$tick()}

  spsp_f <- spsp_ft

  fyears <- tapply(spsp_f$Year,spsp_f$rt.uni.ob,min)
  fyr.df <- data.frame(rt.uni.ob = names(fyears),Year = as.integer(fyears),firstyear = 1)
  spsp_ft <- merge(spsp_f,fyr.df,all.x = TRUE,by = c("rt.uni.ob","Year"))
  if (!isTRUE(quiet)){pb$tick()}

  spsp_f <- spsp_ft

  spsp_f[which(is.na(spsp_f$firstyear)),"firstyear"] <- 0
  spsp_f <- spsp_f[,c("count","stratum","obser","yr","firstyear","strat_name","rt.uni","Year","Month","Day","ObsN","Latitude","Longitude")]
  names(spsp_f) <- c("count","strat","obs_route","year","firstyr","strat_name","route","rYear","Month","Day","ObsN","Latitude","Longitude")
  if (!isTRUE(quiet)){pb$tick()}

  ### supporting new Stan models
  spsp_f$observer <- as.integer(factor(spsp_f$ObsN))
  spsp_f$site <- as.integer(factor(spsp_f$route))
  
  

# Observer route combinations for retransformations -----------------------

  
  obs_routes_df <- unique(spsp_f[,c("site","route","strat","strat_name","observer","ObsN")])
  obs_routes_df <- obs_routes_df[order(obs_routes_df$strat,obs_routes_df$observer),]
  obs_routes_df$obs_route <- as.integer(factor(paste(obs_routes_df$route,obs_routes_df$ObsN)))  
  n_obs_routes_strata <- obs_routes_df %>% arrange(strat,
                                          obs_route) %>% 
    group_by(strat) %>% 
    summarise(nobsroutes = n())
  
  nobs_routes_st <- tapply(obs_routes_df$obs_route,obs_routes_df$strat,length)
  nobs_routes_strata <- data.frame(strat = as.integer(names(nobs_routes_st)),nobs_sites = as.integer(nobs_routes_st))
  
  
  nobs_sites_strata <- as.integer(nobs_routes_strata$nobs_sites)
  maxnobs_sites_strata <- max(nobs_sites_strata)
  nstrata <- max(nobs_routes_strata$strat)
  ste_mat <- matrix(data = 0,
                    nrow = nstrata,
                    ncol = maxnobs_sites_strata)
  obs_mat <- matrix(data = 0,
                    nrow = nstrata,
                    ncol = maxnobs_sites_strata)
  for(i in 1:nstrata){
    ste_mat[i,1:nobs_sites_strata[i]] <- obs_routes_df[which(obs_routes_df$strat == i),"site"]
    obs_mat[i,1:nobs_sites_strata[i]] <- obs_routes_df[which(obs_routes_df$strat == i),"observer"]
  }
  
  
  pR_wts <- pR.wts
  n_observers = max(spsp_f$observer)
  nrts_used = nrts.used
  nsites <- max(spsp_f$site)
  ncounts <- nrow(spsp_f)
  
  if (!isTRUE(quiet)){pb$tick()}

  ####################### END BUGS DATA PREP #######################

  ymin = 1
  ymax = max(spsp_f$year)
  nyears = length(ymin:ymax)
  years = c(ymin:ymax)
  if (!isTRUE(quiet)){pb$tick()}
  
  if(heavy_tailed){
    hvt <- 1
    if(calculate_nu){
      calc_nu <- 1
    }else{
      calc_nu <- 0
    }
  }else{
    hvt <- 0
  }
 if(use_pois){
   use_p <- 1
 }else{
   use_p <- 0
 }
  
  if(calculate_log_lik){
    calc_log_lik <- 1
  }else{
    calc_log_lik <- 0
  }
  if(calculate_CV){
    calc_CV <- 1
  }else{
    calc_CV <- 0
  }
    
  
  
  to_return <- list(model = model,
                    
                    nsites = nsites,
                    nstrata = nstrata,
                    ncounts = ncounts,
                    nyears = nyears,
                    
                    #basic data
                    count = spsp_f$count,
                    strat = spsp_f$strat,
                    year = spsp_f$year,
                    site = spsp_f$site,
                    
                    # #spatial structure
                    # N_edges = N_edges,
                    # node1 = node1,
                    # node2 = node2,
                    
                    
                    #Observer information
                    nobservers = n_observers,
                    observer = spsp_f$observer,
                    firstyr = spsp_f$firstyr,
                    
                    #Ragged array information to link sites and observers to strata 
                    nobs_sites_strata = nobs_sites_strata,
                    maxnobs_sites_strata = maxnobs_sites_strata,
                    ste_mat = ste_mat,
                    obs_mat = obs_mat,
                    
                    #weights
                    nonzeroweight = pR_wts$p.r.ever,
                    
                    # extra Poisson variance options
                    calc_nu = calc_nu, # do not calculate df for the t-distributed noise use df == 3 (if == 1 nu ~ gamma(2,0.1))
                    heavy_tailed = hvt, # use a heavy-tailed t-distributed noise (if == 0 use normal)
                    use_pois = use_p, # if 1 then use over-dispersed poisson else use Negatvei Binomial
                    
                    #loo and Cross validation options
                    calc_log_lik = calc_log_lik, #calculate log_lik point wise log-likelihood of data given model
                    calc_CV = calc_CV, # do data and model include cross-validation training and test data? TRUE = 1, FALSE = 0
                    train = as.integer(1:ncounts), # indices of observations in the training dataset
                    test = as.integer(1), # indices of observations in the test dataset
                    ntrain = ncounts, # number of training data - must == ncounts if calc_CV == 0
                    ntest = 1,  # number of testing data - ignored if calc_CV == 0

                    stratify_by = strat_data$stratify_by,
                    r_year = spsp_f$rYear,
                    alt_data = list(
                      full_data_frame = spsp_f,
                    strat_name = spsp_f$strat_name,
                    obs_route = as.integer(spsp_f$obs_route),
                    nobs_routes = max(spsp_f$obs_route),
                    month = spsp_f$Month,
                    day = spsp_f$Day,
                    Latitude = spsp_f$Latitude,
                    Longitude = spsp_f$Longitude)
                    )
  if (!is.null(model))
  {
    if (tolower(model) %in% c("slope","firstdiff"))
    {
      to_return <- c(to_return,
                     list(fixedyear = floor(stats::median(unique(spsp_f$year)))))
    }

    if (tolower(model) %in% c("gam", "gamye"))
    {
      if (is.null(n_knots))
      {
        n_knots <- floor(length(unique((spsp_f$year)))/4)
      }

      if(tolower(basis) == "mgcv"){
        data_pred <- data.frame(x = years)
        smooth_basis = mgcv::smoothCon(mgcv::s(x,k = n_knots+1, bs = "tp"),data = data_pred,
                                       absorb.cons=TRUE,#this drops the constant and absorbs the identifiability constraints into the basis
                                       diagonal.penalty=TRUE) ## If TRUE then the smooth is reparameterized to turn the penalty into an identity matrix, with the final diagonal elements zeroed (corresponding to the penalty nullspace).

        year_basis = smooth_basis[[1]]$X
        #to_return$model <- paste0(to_return$model,"_new")

      }else{
        recenter = floor(diff(c(1,ymax))/2)
        rescale = ymax # this generates a year variable with range = 1, this rescaling helps the convergence for the GAM beta parameters
        spsp_f$yearscale = (spsp_f$year-recenter)/ymax
        if (!isTRUE(quiet)){pb$tick()}

        scaledyear = seq(min(spsp_f$yearscale),max(spsp_f$yearscale),length = nyears)
        names(scaledyear) <- ymin:ymax
        if(ymin != 1)
        {
          newys = 1:(ymin-1)
          newyscale = (newys-recenter)/rescale
          names(newyscale) <- newys
          scaledyear = c(newyscale,scaledyear)
        }
        if (!isTRUE(quiet)){pb$tick()}

        yminsc = scaledyear[as.character(ymin)]
        ymaxsc = scaledyear[as.character(ymax)]
        if(ymin != 1)
        {
          yminpred = 1
          yminscpred = scaledyear[as.character(1)]
        }
        if (!isTRUE(quiet)){pb$tick()}


        knotsX<- seq(yminsc,ymaxsc,length=(n_knots+2))[-c(1,n_knots+2)]
        X_K<-(abs(outer(seq(yminsc,ymaxsc,length = nyears),knotsX,"-")))^3
        X_OMEGA_all<-(abs(outer(knotsX,knotsX,"-")))^3
        X_svd.OMEGA_all<-svd(X_OMEGA_all)
        X_sqrt.OMEGA_all<-t(X_svd.OMEGA_all$v  %*% (t(X_svd.OMEGA_all$u)*sqrt(X_svd.OMEGA_all$d)))
        year_basis<-t(solve(X_sqrt.OMEGA_all,t(X_K)))


      }

      to_return <- c(to_return, 
                     list(nknots_year = n_knots,
                          year_basis = year_basis))
    }

  }
  if (!isTRUE(quiet)){pb$tick()}

  return(to_return)
}
