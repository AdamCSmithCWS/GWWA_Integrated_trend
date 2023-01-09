## an iCAR model for integrated trend analysis with BBS and species specific monitoring for GWWA


library(bbsBayes)
library(tidyverse)
library(cmdstanr)
# library(rstan)
# rstan_options(auto_write = TRUE, javascript = FALSE)
# library(shinystan)
library(sf)
library(spdep)
library(ggforce)
source("functions/neighbours_define_alt.R") ## function to define neighbourhood relationships
source("functions/prepare-jags-data-alt.R") ## small alteration of the bbsBayes function
source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output
## changes captured in a commit on Nov 20, 2020


# load and stratify GWWA data from BBS ---------------------------------------------
species = "Golden-winged Warbler"
strat = "bbs_usgs"
model = "slope"
scope = "integrated"

strat_data = stratify(by = strat)

firstYear = 2009
lastYear = 2021

#output_dir <- "G:/BBS_iCAR_route_trends/output"
output_dir <- "output"


laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system



# load the GWWA monitoring data -------------------------------------------

gwwa_dat <- read.csv("data/2009_2021_GWWA_Data_final.csv")

gwwa_dat <- gwwa_dat %>% 
  mutate(State = str_trim(str_to_upper(State)),
         GWWA = as.integer(GWWA))  %>% 
  filter(!is.na(GWWA))


## make spatial point layer of gwwa monitoring sites -----------------------

gwwa_map = st_as_sf(gwwa_dat,coords = c("Longitude","Latitude"))
st_crs(gwwa_map) <- 4326 #WGS 84??? - this is a guess
#load strata map
gwwa_map <- st_transform(gwwa_map,crs = laea)

ttp <- ggplot()+
  geom_sf(data = gwwa_map,aes(colour = Year))

print(ttp)

quads <- st_read(dsn = "data",
                 layer = "delorme_conus") 
st_crs(quads) <- 4326 #WGS 84
quads <- st_transform(quads,crs = laea)

## bounding box for plotting quads
gwwa_bounds <- st_union(gwwa_map) #union to provide a simple border of the sampled region
bb = st_bbox(gwwa_bounds)
xlms = as.numeric(c(bb$xmin,bb$xmax))
ylms = as.numeric(c(bb$ymin,bb$ymax))

str_del = function(x){
  y = gsub(x,pattern = "US-",replacement = "",fixed = TRUE)
  y = str_trim(str_sub(y,start = 1,end = -6))
  return(y)
}
quads <- quads %>% 
  mutate(State = str_del(DELORME_ED)) %>% 
  filter(State %in% unique(gwwa_map$State))

#unique(gwwa_map$State)


ttq <- ggplot()+
  geom_sf(data = quads,alpha = 0,aes(colour = State))+
  geom_sf(data = gwwa_map, aes(colour = State))+
  coord_sf(xlim = xlms,ylim = ylms)

print(ttq)



## join points with quads --------------------------------------------------
#quads <- st_transform(quads,crs = st_crs(gwwa_map))
out_strat <- NULL
for(st in unique(gwwa_map$State)){
  tmpq = quads %>% filter(State == st)
  tmpd = gwwa_map %>% filter(State == st)
  
  tmpj <- st_join(tmpd,tmpq,join = st_is_within_distance,dist = 5, left = TRUE)
  
  out_strat <- bind_rows(out_strat,tmpj)
}
for(i in 1:nrow(out_strat)){
  out_strat$page[i] <- ifelse(grepl(out_strat$ID.Code[i],pattern = out_strat$PAGE_NUM[i]),
                              paste0(out_strat$State.x[i],as.character(out_strat$PAGE_NUM[i])),
                              paste0(out_strat$State.y[i],as.character(out_strat$PAGE_NUM[i])))
  
}

ttq <- ggplot()+
  geom_sf(data = quads,alpha = 0)+
  geom_sf(data = out_strat, aes(colour = factor(page)))+
  coord_sf(xlim = xlms,ylim = ylms)

print(ttq)


span_func = function(x){
  diff(range(x))
}
span_strat <- out_strat %>%
  as.data.frame() %>% 
  filter(Year < 2022) %>% 
  group_by(page) %>% 
  summarise(span_years = span_func(Year),
            start_year = min(Year),
            end_year = max(Year),
            n_obs = n())



# Setting Minimum Span for inclusion --------------------------------------
insufficient_span <- 3

strat_w_sufficient_data <- span_strat %>% 
  filter(span_years > insufficient_span,
         !is.na(page)) %>% 
  select(page) %>% 
  as.data.frame()

# gwwa_ag = data to merge with the BBS ----------------

gwwa_ag <- out_strat %>% 
  as.data.frame() %>% 
  filter(page %in% strat_w_sufficient_data$page,
         Year < 2022) %>% 
  group_by(page,Year,Observer) %>% 
  summarise(count = sum(GWWA),
            n_survey = n(),
            raw_mean_count = mean(GWWA),
            observer_gwwa = as.integer(factor(Observer)))


simple_raw_plot <- ggplot(data = gwwa_ag,aes(x = Year,y = raw_mean_count+0.1,group = page))+
  geom_point(aes(colour = page))+
  geom_smooth(method = "lm",se = FALSE)+
  scale_y_continuous(trans = "log10")
print(simple_raw_plot)




# spatial points for page centres -----------------------------------------

gwwa_strats <- quads %>% 
  mutate(page = paste0(State,PAGE_NUM)) %>% 
  filter(page %in% strat_w_sufficient_data$page) %>% 
  select(page) %>%
  group_by(page) %>% #this and the next line performs a union function
  summarise()

strat_map_gwwa <- ggplot()+
  geom_sf(data = gwwa_strats,alpha = 0)

print(strat_map_gwwa)

site_centres_gwwa <- st_centroid(gwwa_strats) %>% 
  mutate(site_gwwa = as.integer(factor(page))) %>% 
  rename(site_orig = page)


# link strata to count data -----------------------------------------------

gwwa_data <- site_centres_gwwa %>% 
  as.data.frame() %>% 
  select(site_orig,site_gwwa) %>% 
  right_join(.,gwwa_ag,by = c("site_orig" = "page")) %>% 
  rename(YEAR = Year) %>% 
  mutate(survey = 0)


# bbs_data <- data.frame(count = jags_data$count,
#                        strat_name = jags_data$strat_name,
#                        YEAR = jags_data$r_year,
#                        route = jags_data$route,
#                        Latitude = jags_data$Latitude,
#                        Longitude = jags_data$Longitude,
#                        obser = jags_data$obser,
#                        firstyr = jags_data$firstyr,
#                        survey = 1,
#                        site_bbs = as.integer(factor(jags_data$route)))
# 


# strat_map_gwwa <- ggplot()+
#   geom_sf(data = gwwa_strats,alpha = 0)+
#   geom_sf_text(data = site_centres_gwwa,aes(label = site_orig,colour = site_orig))
# 
# print(strat_map_gwwa)

# GWWA spatial neighbourhood define --------------------------------------------

# 
# ## returns the adjacency data necessary for the stan model
# ## also exports maps and saved data objects to plot_dir
# car_stan_dat_gwwa <- neighbours_define(real_strata_map = site_centres_gwwa,
#                                   strat_link_fill = 50000,
#                                   plot_neighbours = TRUE,
#                                   species = species,
#                                   plot_dir = "route_maps/",
#                                   plot_file = paste0("_GWWA_",scope,"_route_maps.pdf"),
#                                   save_plot_data = TRUE,
#                                   voronoi = TRUE,
#                                   alt_strat = "site_gwwa",
#                                   add_map = gwwa_strats)
# 



# Load BBS data -----------------------------------------------------------


species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)

sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_with_observer_",firstYear,"_",lastYear,"_slope_iCAR.RData")


jags_data = try(prepare_jags_data(strat_data = strat_data,
                                  species_to_run = species,
                                  model = model,
                                  #n_knots = 10,
                                  min_year = firstYear,
                                  max_year = lastYear,
                                  min_n_routes = 1,
                                  min_max_route_years = 2),silent = TRUE) # 

#create a dataframe of the jags_data results
bbs_data <- data.frame(count = jags_data$count,
                       strat_name = jags_data$strat_name,
                       YEAR = jags_data$r_year,
                       route = jags_data$route,
                       Latitude = jags_data$Latitude,
                       Longitude = jags_data$Longitude,
                       obser = jags_data$obser,
                       observer_bbs = as.integer(factor(jags_data$obser)),
                       firstyr = jags_data$firstyr,
                       survey = 1,
                       site_bbs = as.integer(factor(jags_data$route)))


# spatial neighbourhood define --------------------------------------------

# strata map of one of the bbsBayes base maps
# helps group and set boundaries for the route-level neighbours
strata_map  <- load_map(stratify_by = strat)
strata_map <- st_transform(strata_map,laea)

realized_strata_map = filter(strata_map,ST_12 %in% unique(bbs_data$strat_name))

# Spatial boundaries set up --------------------



route_map_bbs = unique(data.frame(route = bbs_data$route,
                                  site_bbs = bbs_data$site_bbs,
                                  #strat = bbs_data$strat_name,
                                  Latitude = bbs_data$Latitude,
                                  Longitude = bbs_data$Longitude))


# reconcile duplicate spatial locations -----------------------------------
# adhoc way of separating different routes with the same starting coordinates
# this shifts the starting coordinates of teh duplicates by ~1.5km to the North East 
# ensures that the duplicates have a unique spatial location, but remain very close to
# their original location and retain the correct neighbourhood relationships
# these duplicates happen when a "new" route is established because some large proportion
# of the end of a route is changed, but the start-point remains the same
dups = which(duplicated(route_map_bbs[,c("Latitude","Longitude")]))
while(length(dups) > 0){
  route_map_bbs[dups,"Latitude"] <- route_map_bbs[dups,"Latitude"]+0.01 #=0.01 decimal degrees ~ 1km
  route_map_bbs[dups,"Longitude"] <- route_map_bbs[dups,"Longitude"]+0.01 #=0.01 decimal degrees ~ 1km
  dups = which(duplicated(route_map_bbs[,c("Latitude","Longitude")]))
  
}
dups = which(duplicated(route_map_bbs[,c("Latitude","Longitude")])) 
if(length(dups) > 0){stop(paste(spec,"ERROR - At least one duplicate route remains"))}


site_centres_bbs = st_as_sf(route_map_bbs,coords = c("Longitude","Latitude"))
st_crs(site_centres_bbs) <- 4269 #NAD83 commonly used by US federal agencies
#load strata map




site_centres_bbs = st_transform(site_centres_bbs,crs = laea) %>% 
  rename(site_orig = route)



# Merge the data ----------------------------------------------------------

bbs_merge <- bbs_data %>% 
  select(count,YEAR,route,observer_bbs,firstyr,survey,site_bbs) %>% 
  rename(site_orig = route,
         observer = observer_bbs) %>% 
  mutate(inds_bbs = 1:nrow(bbs_data),
         offset = log(50),
         inds_gwwa = 0,
         site_gwwa = 1)

## offset is set to a proportion of the designed number of counts
## this should help to estimate the GWWA intercept, which is currently very small
## if the offset scales the values to individual counts
gwwa_merge <- gwwa_data %>% 
  select(count,YEAR,site_orig,survey,site_gwwa,n_survey,Observer) %>% 
  rename(offset = n_survey) %>% 
  mutate(inds_gwwa = 1:nrow(gwwa_data),
         offset = log(offset),
         observer = as.integer(factor(Observer)),
         firstyr = 0,
         inds_bbs = 0,
         site_bbs = 1)


data_all <- bind_rows(bbs_merge,gwwa_merge) %>% 
  mutate(year = (YEAR-min(YEAR))+1)



site_centres <- bind_rows(site_centres_bbs,site_centres_gwwa) %>% 
  mutate(site = as.integer(factor(site_orig)))


# join site centres to BBS strata -----------------------------------------
site_centres <- st_join(site_centres,
                        strata_map,
                        join = st_nearest_feature)

## returns the adjacency data necessary for the stan model
## also exports maps and saved data objects to plot_dir
car_stan_dat <- neighbours_define(real_strata_map = site_centres,
                                  strat_link_fill = 50000,
                                  plot_neighbours = TRUE,
                                  species = species,
                                  plot_dir = "route_maps/",
                                  plot_file = paste0("_merged_new_",scope,"_route_maps.pdf"),
                                  save_plot_data = TRUE,
                                  voronoi = TRUE,
                                  strat_indicator = "site",
                                  add_map = strata_map)

site_centres2 <- site_centres %>% 
  st_transform(.,crs = 4269) %>% 
  st_coordinates()
  
site_list_temp <- site_centres %>% 
  mutate(Longitude = site_centres2[,1],
         Latitude = site_centres2[,2]) %>% 
  as.data.frame() %>% 
  select(site_orig,site,BCR,ST_12,AREA_1,PROVSTATE,COUNTRY,
         Longitude,Latitude)


data_all <- data_all %>% 
  left_join(.,site_list_temp,by = c("site_orig"))
site_list <- data_all %>% 
  as.data.frame() %>% 
  select(site_orig,site,site_bbs,site_gwwa,survey,ST_12,BCR,COUNTRY,Latitude, Longitude) %>% 
  distinct() %>% 
  arrange(site)


nsites = max(data_all$site)
ncounts = nrow(data_all)
ncounts_bbs = max(data_all$inds_bbs)
ncounts_gwwa = max(data_all$inds_gwwa)
nyears = max(data_all$year)
nobservers_bbs = max(data_all[which(data_all$survey == 1),"observer"])
nobservers_gwwa = max(data_all[which(data_all$survey == 0),"observer"])

site_bbs <- site_list$site_bbs
site_gwwa <- site_list$site_gwwa
survey_sites <- site_list$survey

fixedyear = floor(nyears/2)

stan_data <- list(nsites = nsites,
                  ncounts = ncounts,
                  ncounts_bbs = ncounts_bbs,
                  ncounts_gwwa = ncounts_gwwa,
                  nyears = nyears,
                  nobservers_bbs = nobservers_bbs,
                  nobservers_gwwa = nobservers_gwwa,
                  
                  count = data_all$count,
                  inds_bbs = data_all$inds_bbs,
                  inds_gwwa = data_all$inds_gwwa,
                  year = data_all$year,
                  site = data_all$site,
                  survey = data_all$survey,
                  firstyr = data_all$firstyr,
                  observer = data_all$obser,
                  off_set = data_all$offset,
                  
                  survey_sites = survey_sites,
                  site_bbs = site_bbs,
                  site_gwwa = site_gwwa,
                  
                  fixedyear = fixedyear,
                  
                  N_edges = car_stan_dat$N_edges,
                  node1 = car_stan_dat$node1,
                  node2 = car_stan_dat$node2)




if(car_stan_dat$N != stan_data[["nsites"]]){stop("Some routes are missing from adjacency matrix")}

mod.file = "models/slope_iCAR_integrated_simple.stan"



## compile model
slope_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))

init_def <- function(){ list(noise_raw_bbs = rnorm(stan_data$ncounts_bbs,0,0.1),
                             noise_raw_gwwa = rnorm(stan_data$ncounts_gwwa,0,0.1),
                             alpha_raw = rnorm(stan_data$nsites,0,0.1),
                             ALPHA_bbs = 0,
                             ALPHA_gwwa = 0,
                             BETA = 0,
                             eta = 0,
                             obs_raw_gwwa = rnorm(stan_data$nobservers_gwwa,0,0.1),
                             obs_raw_bbs = rnorm(stan_data$nobservers_bbs,0,0.1),
                             sdnoise_bbs = 0.2,
                             sdnoise_gwwa = 0.2,
                             sdobs_gwwa = 0.1,
                             sdobs_bbs = 0.1,
                             sdbeta_space = runif(1,0.01,0.1),
                             #sdbeta_rand = runif(1,0.01,0.1),
                             sdalpha = runif(1,0.01,0.1),
                             #beta_raw_rand = rnorm(stan_data$nsites,0,0.01),
                             beta_raw_space = rnorm(stan_data$nsites,0,0.01))
  }

out_base <- paste0(species_f,"_",scope,"_",firstYear)

slope_stanfit <- slope_model$sample(
  data=stan_data,
  refresh=500,
  chains=3, 
  iter_warmup=2000,
  iter_sampling=4000,
  parallel_chains = 3,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 15,
  seed = 123,
  init = init_def,
  output_dir = output_dir,
  output_basename = out_base)

tmp <- slope_stanfit$summary()


slope_stanfit$save_object(file = paste0(output_dir,"/",out_base,"_gamye_iCAR.RDS"))



shiny_explore <- FALSE
if(shiny_explore){
  shinystan::launch_shinystan(shinystan::as.shinystan(slope_stanfit))
  
  #loo_stan = loo(slope_stanfit)
}







save(list = c("out_base",
              "stan_data",
              "site_list",
              "strata_map",
              "firstYear",
              "sp_file",
              "species_f",
              "output_dir",
              "car_stan_dat",
              "site_centres",
              "data_all"),
     file = sp_file)





#stopCluster(cl = cluster)

# post analysis summary ---------------------------------------------------



load(sp_file)
slope_stanfit <- readRDS(paste0(output_dir,"/",out_base,"_gamye_iCAR.RDS"))






# PLOTTING and trend output -----------------------------------------------

# library(tidybayes)


route_trajectories <- TRUE #set to FALSE to speed up mapping
# 
# maps = vector(mode = "list",length = 400)
# maps2 = vector(mode = "list",length = 400)
# 
# maps3 = vector(mode = "list",length = 400)
# 
# maps_rand = vector(mode = "list",length = 400)
# maps_space = vector(mode = "list",length = 400)
# 
# trends_out <- NULL
# trends_out_space <- NULL
# trends_out_rand <- NULL
# sdbeta_dif <- NULL
# sdbeta_space_rand <- NULL


LC = 0.05
UC = 0.95
library(HDInterval)


# site-level betas and alphas ---------------------------------------------


betas <- posterior_samples(slope_stanfit,
                                  "beta",
                                  dims = c("site")) %>% 
  posterior_sums(quantiles = NULL,
                 ci = 0.9,#c(0.025,0.5,0.975),
                 dims = "site") %>% 
  select(site,mean,median,lci,uci) %>% 
  rename(beta_mean = mean,
         beta_median = median,
         beta_lci = lci,
         beta_uci = uci) %>% 
  mutate(trend = (exp(beta_median)-1)*100,
         trend_lci = (exp(beta_lci)-1)*100,
         trend_uci =(exp(beta_uci)-1)*100) 
  

alphas <- posterior_samples(slope_stanfit,
                                  "alpha",
                                  dims = c("site"))%>% 
  posterior_sums(quantiles = NULL,
                 ci = 0.9,#c(0.025,0.5,0.975),
                 dims = "site") %>% 
  select(site,mean,median,lci,uci) %>% 
  rename(alpha_mean = mean,
         alpha_median = median,
         alpha_lci = lci,
         alpha_uci = uci) %>% 
  mutate(abundance = exp(alpha_median),
         abundance_lci = exp(alpha_lci),
         abundance_uci =exp(alpha_uci)) 

param_summary <- inner_join(betas,alphas,
                            by = "site") %>% 
  left_join(.,site_list,
            by = "site") %>% 
  rename(country_stateprov_bcr = ST_12) %>% 
  select(-site)

write.csv(param_summary,
          "trends/site_level_model_parameter_estimates.csv")

# obs_mean_counts <- data_all %>% 
#   ungroup() %>% 
#   mutate(ncounts = ifelse(survey == 1,1,exp(offset)/5)) %>% 
#   group_by(site_orig,site,site_bbs,site_gwwa,
#            survey,ST_12,BCR,COUNTRY,Longitude, Latitude) %>% 
#   summarise(obs_mean = mean(count/ncounts),
#             .groups = "keep")






#if(route_trajectories){
  
  
  n_samples <- posterior_samples(slope_stanfit,
                                 "indices",
                                 dims = c("site","y"))
  
  
  # obs_means
  
  obs_means <- data_all %>% 
    ungroup() %>% 
    mutate(ncounts = ifelse(survey == 1,1,exp(offset)/5)) %>% 
    group_by(site_orig,site,site_bbs,site_gwwa,
             survey,ST_12,BCR,COUNTRY,YEAR) %>% 
    summarise(obs_mean = mean(count/ncounts),
              .groups = "keep")
  
  
  n_samples_t <- n_samples %>% 
    mutate(YEAR = y+(firstYear-1)) %>% 
    full_join(.,obs_means,by = c("site","YEAR"))
  
  ## figure out how to plot these trajectories with observed means                                                   
  
  
site_indices <- n_samples_t %>% 
  group_by(site_orig,survey,BCR,ST_12,YEAR) %>%
  summarise(ind = median(.value),
            lci = hdi(.value,0.90)[1],
            uci = hdi(.value,0.90)[2],
            obs_means = mean(obs_mean),
            .groups = "keep") %>% 
  rename(year = YEAR)
  

nsites_strat <- site_indices %>% 
  ungroup() %>% 
  select(site_orig,ST_12,survey) %>% 
  distinct() %>% 
  group_by(ST_12) %>% 
  summarise(nsites = n())

pdf("Figures/site_trajectories.pdf",
    height = 8.5,
    width = 11)
for(s in nsites_strat$ST_12){
  if(is.na(s)){next}
  tmp <- site_indices %>% 
    filter(ST_12 == s)
  
  tmppp <- ggplot(data = tmp, aes(x = year,y = ind))+
    geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.2)+
    geom_line()+
    geom_point(aes(x = year,y = obs_means),colour = "lightblue")+
 #   geom_point(aes(x = year,y = obs_means2),colour = "pink")+
    scale_y_continuous(limits = c(0,NA))+
    scale_x_continuous(limits = c(2009,2021))+
    facet_wrap(facets = vars(site_orig),scales = "free_y")
    
  print(tmppp)
}
dev.off()



# plotting indices by BCR -------------------------------------------------

n_samples <- posterior_samples(slope_stanfit,
                               "indices_plot",
                               dims = c("site","y"))


# obs_means

sites_list <- data_all %>% 
  ungroup() %>% 
  select(site_orig,site,site_bbs,site_gwwa,
           survey,ST_12,BCR,COUNTRY) %>% 
  distinct()

nsites_bcr <- sites_list %>% 
  group_by(survey,BCR) %>% 
  summarise(n_sites = n()) %>% 
  mutate(survey = ifelse(survey == 1,"BBS","GWWA")) %>% 
  pivot_wider(id_cols = BCR,
              names_from = survey,
              values_from = n_sites,
              names_prefix = "n_sites_",
              values_fill = 0) 

n_samples_t <- n_samples %>% 
  mutate(year = y+(firstYear-1)) %>% 
  left_join(.,sites_list,by = c("site"))


  
  
  BCR_indices <- n_samples_t %>% 
    group_by(.draw,BCR,year) %>% 
    summarise(comp_ind = mean(.value), # mean across all sites for each draw
              nsites = n(),
              .groups = "drop") %>% 
    group_by(BCR,year) %>% 
    summarise(ind = median(comp_ind),
              lci = hdi(comp_ind,0.90)[1],
              uci = hdi(comp_ind,0.90)[2],
              nsites = min(nsites),
              nsites2 = max(nsites),
              .groups = "keep") %>% 
    left_join(.,nsites_bcr,by = "BCR")
  
  bcr_plot <- ggplot(data = BCR_indices,aes(x = year, y = ind))+
    geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.2)+
    geom_line(aes(alpha = nsites))+
    scale_alpha_continuous(trans = "sqrt")+
    facet_wrap(~BCR,nrow = 4,scales = "free_y")
  
  
  print(bcr_plot)
  
  posterior_trends <- function(n_samples = n_samples_t,
                               startyear = firstYear,
                               endyear = lastYear,
                               region = "Overall",
                               lq = 0.025,
                               uq = 0.975){
    
    texp <- function(x,ny = 10){
      (x^(1/ny)-1)*100
    }
    
    
    
    
    chng <- function(x){
      (x-1)*100
    }
    
    prob_dec <- function(ch,thresh){
      
      length(which(ch < thresh))/length(ch)
    }
    
    
    nyrs <- endyear-startyear
    syr = paste(startyear)
    eyr = paste(endyear)
    
    if(region == "Overall"){
    n_samples <- n_samples %>% 
      mutate(Overall = "Overall")
    }
    
    out_trends <- n_samples %>% 
      rename_with(~gsub(pattern = region,replacement = "region",.x,fixed = TRUE))  %>% 
      filter(year %in% c(startyear,endyear)) %>% 
      group_by(.draw,region,year) %>% 
      summarise(comp_ind = mean(.value),
                .groups = "drop") %>%  
      pivot_wider(names_from = year,
                  values_from = comp_ind) %>%
      rename_with(~gsub(syr,"startyear",.x,fixed = TRUE)) %>% 
      rename_with(~gsub(eyr,"endyear",.x,fixed = TRUE)) %>% 
      group_by(region,.draw) %>% 
      summarise(t = texp(endyear/startyear,ny = nyrs),
                ch = chng(endyear/startyear),
                .groups = "drop") %>% 
      group_by(region) %>% 
      summarise(trend = mean(t),
                lci = quantile(t,lq,names = FALSE),
                uci = quantile(t,uq,names = FALSE),
                percent_change = median(ch),
                p_ch_lci = quantile(ch,lq,names = FALSE),
                p_ch_uci = quantile(ch,uq,names = FALSE),
                prob_decline = prob_dec(ch,0),
                prob_decline_GT30 = prob_dec(ch,-30),
                prob_decline_GT50 = prob_dec(ch,-50),
                prob_decline_GT70 = prob_dec(ch,-70),
                .groups = "keep") %>% 
      rename_with(~gsub(pattern = "region",replacement = region,.x,fixed = TRUE))
    
    
    return(out_trends)
  }
  
  # BCR_trends <- posterior_trends(region = "BCR") %>%
  #   left_join(.,nsites_bcr,by = "BCR")
  # 
  #strata_trends <- posterior_trends(region = "ST_12") 
  
  survey_trends <- posterior_trends(region = "survey") 
  
  overall_trends <- posterior_trends()
  
  site_trends <- posterior_trends(region = "site_orig")
  
  trends_out <- bind_rows(overall_trends,survey_trends,site_trends)
  write.csv(trends_out,paste0("trends/",species_f,"_","site_trends_site_model_",firstYear,".csv"),
            row.names = FALSE)
  
  
  
  
  
  





  
  #sp_file <- paste0("output/",species,"Canadian_",firstYear,"_",lastYear,"_slope_route_iCAR2.RData")

   # if(species == "Northern Cardinal"){next
    #   
    # #csv_files <- dir(output_dir,pattern = out_base,full.names = TRUE)
    # }
shinycheck <- FALSE
    if(shinycheck){
    ### may be removed after re-running     launch_shinystan(slope_stanfit)
    slope_stanfit1 <- rstan::read_stan_csv(csv_files)
    launch_shinystan(as.shinystan(slope_stanfit1))
    ###
    }


    slope_stanfit <- slope_stanfit
     ####
    # add trend and abundance ----------------------------------------
    
    beta_samples = posterior_samples(slope_stanfit,"beta",
                                     dims = "s")
    # beta_samples2 = posterior_samples(slope_stanfit,"beta",
    #                                  dims = "s")
    
    slopes = beta_samples %>% group_by(s) %>% 
      summarise(b = mean(.value),
                lci = quantile(.value,LC),
                uci = quantile(.value,UC),
                sd = sd(.value),
                prec = 1/var(.value),
                trend = mean((exp(.value)-1)*100),
                lci_trend = quantile((exp(.value)-1)*100,LC),
                uci_trend = quantile((exp(.value)-1)*100,UC),
                .groups = "keep")
    
    alpha_samples = posterior_samples(slope_stanfit,"alpha",
                                      dims = "s")
    interc = alpha_samples %>% group_by(s) %>% 
      summarise(abund = mean(exp(.value)),
                lci_i = quantile(exp(.value),LC),
                uci_i = quantile(exp(.value),UC),
                sd_i = sd(exp(.value)),
                prec_i = 1/var(.value),
                .groups = "keep")
    
    #plot(log(interc$i),slopes$b)
    slops_int = inner_join(slopes,interc,by = "s")
    slops_int$site = slops_int$s
    
    
    
    # random effect plus mean component of slope ----------------------------------------
    
    # BETA_samples = posterior_samples(slope_stanfit,BETA) %>% 
    #   rename(BETA = .value) %>% 
    #   ungroup() %>% 
    #   select(BETA,.draw)
    # 
    # beta_rand_samples = posterior_samples(slope_stanfit,beta_rand[s]) %>% 
    #   rename(beta_rand = .value) %>% 
    #   ungroup() %>% 
    #   select(beta_rand,.draw,s)
    # 
    # beta_rand_samples <- inner_join(beta_rand_samples,BETA_samples,by = c(".draw"))
    # 
    # slopes_rand_full = beta_rand_samples %>% group_by(s) %>% 
    #   summarise(b = mean(beta_rand + BETA),
    #             lci = quantile(beta_rand + BETA,LC),
    #             uci = quantile(beta_rand + BETA,UC),
    #             sd = sd(beta_rand + BETA),
    #             prec = 1/var(beta_rand + BETA),
    #             .groups = "keep")
    # 
    # slopes_rand_full_int = inner_join(slopes_rand_full,interc,by = "s")
    # slopes_rand_full_int$site = slopes_rand_full_int$s
    
    # 
    # beta_rand_samples = posterior_samples(slope_stanfit,
    #                                       "beta_rand",
    #                                       dims = "s")
    # 
    # slopes_rand = beta_rand_samples %>% group_by(s) %>% 
    #   summarise(b = mean(.value),
    #             lci = quantile(.value,LC),
    #             uci = quantile(.value,UC),
    #             sd = sd(.value),
    #             prec = 1/var(.value),
    #             trend = mean((exp(.value)-1)*100),
    #             lci_trend = quantile((exp(.value)-1)*100,LC),
    #             uci_trend = quantile((exp(.value)-1)*100,UC),
    #             .groups = "keep")
    # 
    # slops_rand_int = inner_join(slopes_rand,interc,by = "s")
    # slops_rand_int$site = slops_rand_int$s
    
    
    # spatial component of slope ----------------------------------------
    
    
    # beta_space_samples = posterior_samples(slope_stanfit,"beta_space",
    #                                        dims = "s")
    # 
    # slopes_space = beta_space_samples %>% group_by(s) %>% 
    #   summarise(b = mean(.value),
    #             lci = quantile(.value,LC),
    #             uci = quantile(.value,UC),
    #             sd = sd(.value),
    #             prec = 1/var(.value),
    #             trend = mean((exp(.value)-1)*100),
    #             lci_trend = quantile((exp(.value)-1)*100,LC),
    #             uci_trend = quantile((exp(.value)-1)*100,UC),
    #             .groups = "keep")
    # 
    # slops_space_int = inner_join(slopes_space,interc,by = "s")
    # slops_space_int$site = slops_space_int$s
    # 
    
    # Compare spatial and random variation ------------------------------------
    # sdbeta_rand_tmp_samples <- posterior_samples(slope_stanfit,
    #                                              "sdbeta_rand")
    # sdbeta_space_tmp_samples <- posterior_samples(slope_stanfit,
    #                                               "sdbeta_space")
    # 
    # sdbeta_space_rand_tmp_samples <- bind_rows(sdbeta_rand_tmp_samples,
    #                                            sdbeta_space_tmp_samples)
    # 
    # 
    # sdbeta_space_rand_tmp <- sdbeta_space_rand_tmp_samples %>% 
    #   group_by(.variable) %>%
    #   summarise(mean = mean((.value)),
    #             lci = quantile((.value),LC),
    #             uci = quantile((.value),UC),
    #             sd = sd((.value)),
    #             .groups = "keep") %>% 
    #   mutate(species = species)
    # #combines all species estimates
    # sdbeta_space_rand <- bind_rows(sdbeta_space_rand,sdbeta_space_rand_tmp)
    # 
    # 
    # 
    # # difference rand-spatial -------------------------------------------------
    # 
    # sdbeta_space_tmp_samples <- sdbeta_space_tmp_samples %>% 
    #   rename(sd_space = .value) %>% 
    #   ungroup() %>% 
    #   select(-.variable)
    # 
    # sdbeta_rand_tmp_samples <- sdbeta_rand_tmp_samples %>% 
    #   rename(sd_rand = .value)%>% 
    #   ungroup() %>% 
    #   select(-.variable)
    # 
    # sdbeta_tmp_samples <- inner_join(sdbeta_rand_tmp_samples,sdbeta_space_tmp_samples)
    # 
    # sdbeta_tmp_dif <- sdbeta_tmp_samples %>% 
    #   group_by(.draw) %>%
    #   summarise(dif = sd_rand-sd_space) %>% 
    #   ungroup() %>% 
    #   summarise(mean = mean((dif)),
    #             lci = quantile((dif),LC),
    #             uci = quantile((dif),UC),
    #             sd = sd((dif))) %>% 
    #   mutate(species = species)
    # 
    # 
    # 
    
    
    route_trajectories <- FALSE
    
    
    # Route-level trajectories ------------------------------------------------
  
    # ind_samples <- posterior_samples(slope_stanfit,
    #                                     "indices",
    #                                  dims = c("s","y")) %>% 
    #   rename(site = s)
    # 
    # I_samples <- posterior_samples(slope_stanfit,
    #                                "I",
    #                                dims = c("y"))
    # 
    # BETA <- posterior_samples(slope_stanfit,
    #                                   "BETA") %>% 
    #   summarise(trend = mean((exp(.value)-1)*100),
    #             lci = quantile((exp(.value)-1)*100,LC),
    #             uci = quantile((exp(.value)-1)*100,UC),
    #             sd = sd((exp(.value)-1)*100))
    # 
    # I_all <- I_samples %>% group_by(y) %>% 
    #   summarise(index = mean(.value),
    #             lci = quantile(.value,LC),
    #             uci = quantile(.value,UC),
    #             sd = sd(.value),
    #             .groups = "keep") %>% 
    #   mutate(scale = "SurveyWide",
    #          YEAR = y+(min(data_all$YEAR)-1))
    # 
    # overall_traj <- ggplot(data = I_all,aes(x = YEAR,y = index))+
    #   geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.1)+
    #   geom_line()+
    #   scale_y_continuous(limits = c(0,NA))
    # print(overall_traj)
    # 
    # inds_all <- ind_samples %>% left_join(.,site_list,by = "site") %>% 
    #   group_by(site, site_orig,survey,y) %>% 
    #   summarise(index = mean(.value),
    #             lci = quantile(.value,LC),
    #             uci = quantile(.value,UC),
    #             sd = sd(.value),
    #             .groups = "keep")%>% 
    #   mutate(scale = "Site",
    #          YEAR = y+(min(data_all$YEAR)-1)) 
    # 
    # 
    # inds_survey <- ind_samples %>% left_join(.,site_list,by = "site") %>% 
    #   group_by(survey,y,.draw) %>% 
    #   summarise(v = mean(.value),
    #             .groups = "drop") %>% 
    #   group_by(survey,y) %>% 
    #   summarise(index = mean(v),
    #             lci = quantile(v,LC),
    #             uci = quantile(v,UC),
    #             sd = sd(v),
    #             .groups = "keep")%>% 
    #   mutate(scale = "Survey",
    #          YEAR = y+(min(data_all$YEAR)-1),
    #          Survey = ifelse(survey == 1,"BBS","GWWA")) %>% 
    #   filter(!(Survey == "BBS" & YEAR >2019))
    # 
    # 
    # survey_trajs <- ggplot(data = inds_survey,aes(x = YEAR,y = index))+
    #   geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.1)+
    #   geom_line()+
    #   scale_y_continuous(limits = c(0,NA))+
    #   facet_wrap(~Survey)
    # print(survey_trajs)
    # # connect trends to original route names ----------------------------------
    # 
    route_map_out = left_join(site_centres,slops_int,by = "site")
    route_map_out$species <- species
    

    
    # route_map_out_rand = left_join(site_centres,slops_rand_int,by = "site")
    # route_map_out_rand$species <- species
    # 
    # 
    # 
    # route_map_out_space = left_join(site_centres,slops_space_int,by = "site")
    # route_map_out_space$species <- species
    # 
    # 
    
    
    ### setting up boundaries for plots
    # load(paste0("route_maps/",species_f,"_route_data.RData"))
    
    site_bounds <- st_union(site_centres) #union to provide a simple border of the realised strata
    bb = st_bbox(site_bounds)
    xlms = as.numeric(c(bb$xmin,bb$xmax))
    ylms = as.numeric(c(bb$ymin,bb$ymax))
    
    bcr_map <- load_map(stratify_by = "bcr")
    bcr_map <- st_transform(bcr_map,crs = laea)
    
    
    # add mapping of trends ---------------------------------------------------
    
    plot_trend <- TRUE #set to false to plot the slopes
    if(plot_trend){
      breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
      lgnd_head <- "Trend\n"
      trend_title <- "trends"
      labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
      labls = paste0(labls, " %/year")
      route_map_out$Tplot <- cut(route_map_out$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
      route_map_out <- route_map_out %>% 
        mutate(h_ci = (uci_trend-lci_trend)/2)
      
      #route_map_out_space$Tplot <- cut(route_map_out_space$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
      #route_map_out_rand$Tplot <- cut(route_map_out_rand$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
      
      
    }else{
      breaks <- c(-0.07, -0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04, 0.07)
      lgnd_head <- "slope\n"
      trend_title <- "trend-slopes"
      labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
      labls = paste0(labls, " slope")
      route_map_out$Tplot <- cut(route_map_out$b,breaks = c(-Inf, breaks, Inf),labels = labls)
      route_map_out <- route_map_out %>% 
        mutate(h_ci = (uci-lci)/2)
      
      route_map_out_space$Tplot <- cut(route_map_out_space$b,breaks = c(-Inf, breaks, Inf),labels = labls)
      route_map_out_rand$Tplot <- cut(route_map_out_rand$b,breaks = c(-Inf, breaks, Inf),labels = labls)
      
    }
    map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                     "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
    names(map_palette) <- labls
    
    route_map_out <- route_map_out %>% 
      mutate(survey = ifelse(is.na(site_gwwa),"BBS","GWWA"),
             log_abund = log(abund),
             rel_abund = abund-mean(abund))
    
    tmap = ggplot(route_map_out)+
      #geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
      geom_sf(data = bcr_map,colour = gray(0.8),fill = NA)+
      geom_sf(aes(colour = Tplot,size = abund,shape = survey))+
      scale_size_continuous(range = c(0.5,2.5),
                            name = "Relative abundance",
                            trans = "log10")+
      scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                          guide = guide_legend(reverse=TRUE),
                          name = paste0(lgnd_head,firstYear,"-",lastYear))+
      coord_sf(xlim = xlms,ylim = ylms)+
      labs(title = paste("DRAFT ",species,trend_title,"by BBS route"),
           subtitle = "Route-level trends from a spatial iCAR model, using Stan")
    
    
    png(filename = paste0("Figures/images/",species_f,"_Trends_w_observer_",firstYear,".png"),
        res = 600,
        width = 20,
        height = 15,
        units = "cm")
    print(tmap)
    dev.off()
    
    
    abunddist <- ggplot(data = route_map_out,aes(x = abund))+
      geom_histogram()+
      facet_wrap(~survey,nrow = 1)
    
    print(abunddist)
    
  
    print(species)
 
  

