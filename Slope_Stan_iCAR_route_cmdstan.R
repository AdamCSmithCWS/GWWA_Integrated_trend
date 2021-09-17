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
source("functions/neighbours_define.R") ## function to define neighbourhood relationships
source("functions/prepare-jags-data-alt.R") ## small alteration of the bbsBayes function
source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output
## changes captured in a commit on Nov 20, 2020


# load and stratify GWWA data from BBS ---------------------------------------------
species = "Golden-winged Warbler"
strat = "bbs_usgs"
model = "slope"
scope = "integrated"

strat_data = stratify(by = strat)

firstYear = 2000
lastYear = 2021

#output_dir <- "G:/BBS_iCAR_route_trends/output"
output_dir <- "output"


laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system



# load the GWWA monitoring data -------------------------------------------

gwwa_dat <- read.csv("data/All_stacked_GWWA_data_2009-2021.csv")
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

strat_w_sufficient_data <- span_strat %>% 
  filter(span_years > 0,
         !is.na(page)) %>% 
  select(page) %>% 
  as.data.frame()

# gwwa_ag = data to merge with the BBS ----------------

gwwa_ag <- out_strat %>% 
  as.data.frame() %>% 
  filter(page %in% strat_w_sufficient_data$page,
         Year < 2022) %>% 
  group_by(page,Year) %>% 
  summarise(count = sum(GWWA),
            n_survey = n(),
            raw_mean_count = mean(GWWA))


simple_raw_plot <- ggplot(data = gwwa_ag,aes(x = Year,y = raw_mean_count,group = page))+
  geom_point(aes(colour = page))+
  geom_smooth(method = "lm",se = FALSE)
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

sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")


jags_data = try(prepare_jags_data(strat_data = strat_data,
                                  species_to_run = species,
                                  model = model,
                                  #n_knots = 10,
                                  min_year = firstYear,
                                  max_year = lastYear,
                                  min_n_routes = 1),silent = TRUE) # 

#create a dataframe of the jags_data results
bbs_data <- data.frame(count = jags_data$count,
                       strat_name = jags_data$strat_name,
                       YEAR = jags_data$r_year,
                       route = jags_data$route,
                       Latitude = jags_data$Latitude,
                       Longitude = jags_data$Longitude,
                       obser = jags_data$obser,
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
  select(count,YEAR,route,obser,firstyr,survey,site_bbs) %>% 
  rename(site_orig = route) %>% 
  mutate(inds_bbs = 1:nrow(bbs_data),
         offset = 0,
         inds_gwwa = 0,
         site_gwwa = 1)


gwwa_merge <- gwwa_data %>% 
  select(count,YEAR,site_orig,survey,site_gwwa,n_survey) %>% 
  rename(offset = n_survey) %>% 
  mutate(inds_gwwa = 1:nrow(gwwa_data),
         offset = log(offset),
         obser = 1,
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
                                  plot_file = paste0("_merged_",scope,"_route_maps.pdf"),
                                  save_plot_data = TRUE,
                                  voronoi = TRUE,
                                  alt_strat = "site",
                                  add_map = strata_map)

site_list_temp <- site_centres %>% 
  as.data.frame() %>% 
  select(site_orig,site,BCR,ST_12,AREA_1,PROVSTATE,COUNTRY)


data_all <- data_all %>% 
  left_join(.,site_list_temp,by = c("site_orig"))
site_list <- data_all %>% 
  as.data.frame() %>% 
  select(site_orig,site,site_bbs,site_gwwa,survey,ST_12,BCR,COUNTRY) %>% 
  distinct() %>% 
  arrange(site)


nsites = max(data_all$site)
ncounts = nrow(data_all)
ncounts_bbs = max(data_all$inds_bbs)
ncounts_gwwa = max(data_all$inds_gwwa)
nyears = max(data_all$year)
nobservers = max(data_all$obser)

site_bbs <- site_list$site_bbs
site_gwwa <- site_list$site_gwwa
survey_sites <- site_list$survey

fixedyear = floor(nyears/2)

stan_data <- list(nsites = nsites,
                  ncounts = ncounts,
                  ncounts_bbs = ncounts_bbs,
                  ncounts_gwwa = ncounts_gwwa,
                  nyears = nyears,
                  nobservers = nobservers,
                  
                  count = data_all$count,
                  inds_bbs = data_all$inds_bbs,
                  inds_gwwa = data_all$inds_gwwa,
                  year = data_all$year,
                  site = data_all$site,
                  survey = data_all$survey,
                  firstyr = data_all$firstyr,
                  observer = data_all$obser,
                  offset = data_all$offset,
                  
                  survey_sites = survey_sites,
                  site_bbs = site_bbs,
                  site_gwwa = site_gwwa,
                  
                  fixedyear = fixedyear,
                  
                  N_edges = car_stan_dat$N_edges,
                  node1 = car_stan_dat$node1,
                  node2 = car_stan_dat$node2)




if(car_stan_dat$N != stan_data[["nsites"]]){stop("Some routes are missing from adjacency matrix")}

mod.file = "models/slope_iCAR_integrated.stan"



## compile model
slope_model <- cmdstan_model(mod.file)

init_def <- function(){ list(noise_raw_bbs = rnorm(stan_data$ncounts_bbs,0,0.1),
                             noise_raw_gwwa = rnorm(stan_data$ncounts_gwwa,0,0.1),
                             alpha_raw = rnorm(stan_data$nsites,0,0.1),
                             ALPHA_bbs = 0,
                             ALPHA_gwwa = 0,
                             BETA = 0,
                             eta = 0,
                             obs_raw = rnorm(stan_data$nobservers,0,0.1),
                             sdnoise_bbs = 0.2,
                             sdnoise_gwwa = 0.2,
                             sdobs = 0.1,
                             sdbeta_space = runif(1,0.01,0.1),
                             sdbeta_rand = runif(1,0.01,0.1),
                             beta_raw_space = rnorm(stan_data$nsites,0,0.01),
                             beta_raw_rand = rnorm(stan_data$nsites,0,0.01))}

out_base <- paste0(species_f,"_",scope,"_",firstYear)

slope_stanfit <- slope_model$sample(
  data=stan_data,
  refresh=200,
  chains=3, iter_sampling=2000,
  iter_warmup=2000,
  parallel_chains = 3,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 14,
  seed = 123,
  init = init_def,
  output_dir = output_dir,
  output_basename = out_base)




csv_files <- dir(output_dir,pattern = out_base,full.names = TRUE)
csv_files <- csv_files[1:3]
#slope_stanfit$save_object(file = paste0(output_dir,"/",out_base,"_gamye_iCAR.RDS"))



shiny_explore <- FALSE
if(shiny_explore){
  sl_rstan <- rstan::read_stan_csv(csv_files)
  shinystan::launch_shinystan(shinystan::as.shinystan(sl_rstan))
  
  loo_stan = loo(sl_rstan)
}





# slope_model = stan_model(file=mod.file)
# 
# print(paste(firstYear,species))
# ## run sampler on model, data
# slope_stanfit <- sampling(slope_model,
#                                data=stan_data,
#                                verbose=TRUE, refresh=100,
#                                chains=4, iter=900,
#                                warmup=600,
#                                cores = 4,
#                                pars = parms,
#                                control = list(adapt_delta = 0.8,
#                                               max_treedepth = 15))
# 

save(list = c("slope_stanfit",
              "out_base",
              "stan_data",
              "site_list",
              "strata_map",
              "firstYear",
              "sp_file",
              "species_f",
              "csv_files",
              "output_dir",
              "car_stan_dat",
              "site_centres"),
     file = sp_file)







#stopCluster(cl = cluster)








# post loop analysis ------------------------------------------------------


# 
# launch_shinystan(slope_stanfit) 
# 
# 
# library(loo)
# library(tidyverse)
# 
# log_lik_1 <- extract_log_lik(slope_stanfit, merge_chains = FALSE)
# r_eff <- relative_eff(exp(log_lik_1), cores = 10)
# loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 10)
# print(loo_1)
# 
# doy = ((jags_data$month-4)*30+jags_data$day)
# plot(loo_1$pointwise[,"influence_pareto_k"],log(stan_data$count+1))
# plot(loo_1$pointwise[,"influence_pareto_k"],doy)
# plot(doy,log(stan_data$count+1))
# 
# 
# 
# loo2 = data.frame(loo_1$pointwise)
# 
# loo2$flag = cut(loo2$influence_pareto_k,breaks = c(0,0.5,0.7,1,Inf))
# dts = data.frame(count = stan_data$count,
#                  obser = stan_data$obser,
#                  route = stan_data$route,
#                  year = stan_data$year)
# loo2 = cbind(loo2,dts)
# 
# plot(log(loo2$count+1),loo2$influence_pareto_k)
# 
# obserk = loo2 %>% group_by(obser) %>% 
#   summarise(n = log(n()),
#             mean_k = mean(influence_pareto_k),
#             max_k = max(influence_pareto_k),
#             sd_k = sd(influence_pareto_k),
#             mean_looic = mean(looic),
#             mean_ploo = mean(p_loo))
# plot(obserk$n,obserk$max_k)
# plot(obserk$n,obserk$mean_k)
# plot(obserk$n,obserk$sd_k)
# plot(obserk$n,obserk$mean_looic)
# plot(obserk$n,obserk$mean_ploo)
# 
# 
# yeark = loo2 %>% group_by(year) %>% 
#   summarise(n = n(),
#             mean_k = mean(influence_pareto_k),
#             q90 = quantile(influence_pareto_k,0.9),
#             max_k = max(influence_pareto_k),
#             sd_k = sd(influence_pareto_k),
#             route = mean(route),
#             sd = sd(route))
# plot(yeark$year,yeark$max_k)
# plot(yeark$year,yeark$mean_k)
# plot(yeark$year,yeark$sd_k)
# plot(yeark$year,yeark$q90)
# 
# routek = loo2 %>% group_by(route) %>% 
#   summarise(n = n(),
#             mean_k = mean(influence_pareto_k),
#             q90_k = quantile(influence_pareto_k,0.9),
#             max_k = max(influence_pareto_k),
#             sd_k = sd(influence_pareto_k),
#             route = mean(route),
#             sd = sd(route))
# plot(routek$route,routek$max_k)
# plot(routek$n,routek$mean_k)
# 
# plot(routek$route,routek$mean_k)
# plot(routek$route,routek$sd_k)
# plot(routek$route,routek$q90_k)
# 
# 


# PLOTTING and trend output -----------------------------------------------

# library(tidybayes)


route_trajectories <- FALSE #set to FALSE to speed up mapping
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




  
  #sp_file <- paste0("output/",species,"Canadian_",firstYear,"_",lastYear,"_slope_route_iCAR2.RData")

   # if(species == "Northern Cardinal"){next
    #   
    # #csv_files <- dir(output_dir,pattern = out_base,full.names = TRUE)
    # }

    ### may be removed after re-running     launch_shinystan(slope_stanfit)
    #sl_rstan <- rstan::read_stan_csv(csv_files)
    #launch_shinystan(as.shinystan(sl_rstan))
    sl_rstan <- slope_stanfit
     ####
    # add trend and abundance ----------------------------------------
    
    beta_samples = posterior_samples(sl_rstan,"beta",
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
    
    alpha_samples = posterior_samples(sl_rstan,"alpha",
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
    
    # BETA_samples = posterior_samples(sl_rstan,BETA) %>% 
    #   rename(BETA = .value) %>% 
    #   ungroup() %>% 
    #   select(BETA,.draw)
    # 
    # beta_rand_samples = posterior_samples(sl_rstan,beta_rand[s]) %>% 
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
    
    
    beta_rand_samples = posterior_samples(sl_rstan,
                                          "beta_rand",
                                          dims = "s")
    
    slopes_rand = beta_rand_samples %>% group_by(s) %>% 
      summarise(b = mean(.value),
                lci = quantile(.value,LC),
                uci = quantile(.value,UC),
                sd = sd(.value),
                prec = 1/var(.value),
                trend = mean((exp(.value)-1)*100),
                lci_trend = quantile((exp(.value)-1)*100,LC),
                uci_trend = quantile((exp(.value)-1)*100,UC),
                .groups = "keep")
    
    slops_rand_int = inner_join(slopes_rand,interc,by = "s")
    slops_rand_int$site = slops_rand_int$s
    
    
    # spatial component of slope ----------------------------------------
    
    
    beta_space_samples = posterior_samples(sl_rstan,"beta_space",
                                           dims = "s")
    
    slopes_space = beta_space_samples %>% group_by(s) %>% 
      summarise(b = mean(.value),
                lci = quantile(.value,LC),
                uci = quantile(.value,UC),
                sd = sd(.value),
                prec = 1/var(.value),
                trend = mean((exp(.value)-1)*100),
                lci_trend = quantile((exp(.value)-1)*100,LC),
                uci_trend = quantile((exp(.value)-1)*100,UC),
                .groups = "keep")
    
    slops_space_int = inner_join(slopes_space,interc,by = "s")
    slops_space_int$site = slops_space_int$s
    
    
    # Compare spatial and random variation ------------------------------------
    sdbeta_rand_tmp_samples <- posterior_samples(sl_rstan,
                                                 "sdbeta_rand")
    sdbeta_space_tmp_samples <- posterior_samples(sl_rstan,
                                                  "sdbeta_space")
    
    sdbeta_space_rand_tmp_samples <- bind_rows(sdbeta_rand_tmp_samples,
                                               sdbeta_space_tmp_samples)
    
    
    sdbeta_space_rand_tmp <- sdbeta_space_rand_tmp_samples %>% 
      group_by(.variable) %>%
      summarise(mean = mean((.value)),
                lci = quantile((.value),LC),
                uci = quantile((.value),UC),
                sd = sd((.value)),
                .groups = "keep") %>% 
      mutate(species = species)
    #combines all species estimates
    sdbeta_space_rand <- bind_rows(sdbeta_space_rand,sdbeta_space_rand_tmp)
    
    
    
    # difference rand-spatial -------------------------------------------------
    
    sdbeta_space_tmp_samples <- sdbeta_space_tmp_samples %>% 
      rename(sd_space = .value) %>% 
      ungroup() %>% 
      select(-.variable)
    
    sdbeta_rand_tmp_samples <- sdbeta_rand_tmp_samples %>% 
      rename(sd_rand = .value)%>% 
      ungroup() %>% 
      select(-.variable)
    
    sdbeta_tmp_samples <- inner_join(sdbeta_rand_tmp_samples,sdbeta_space_tmp_samples)
    
    sdbeta_tmp_dif <- sdbeta_tmp_samples %>% 
      group_by(.draw) %>%
      summarise(dif = sd_rand-sd_space) %>% 
      ungroup() %>% 
      summarise(mean = mean((dif)),
                lci = quantile((dif),LC),
                uci = quantile((dif),UC),
                sd = sd((dif))) %>% 
      mutate(species = species)
    
   
    
    
    
    route_trajectories <- TRUE
    
    
    # Route-level trajectories ------------------------------------------------
    if(route_trajectories){
      
      
      
      n_samples <- posterior_samples(sl_rstan,
                                                   "indices",
                                     dims = c("site","y"))
      
      
      n_samples_t <- n_samples %>% 
        left_join(site_list,by = "site") %>% 
        mutate(year = y+(firstYear-1))
      
      
      BCR_indices <- n_samples_t %>% 
        group_by(.draw,BCR,year) %>% 
        summarise(comp_ind = mean(.value),
                  .groups = "drop") %>% 
        group_by(BCR,year) %>% 
        summarise(ind = median(comp_ind),
                  lci = quantile(comp_ind,0.025),
                  uci = quantile(comp_ind,0.975),
                  .groups = "keep")
      
      bcr_plot <- ggplot(data = BCR_indices,aes(x = year, y = ind))+
        geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.2)+
        geom_line()+
        facet_wrap(~BCR,nrow = 4)
      
      
      print(bcr_plot)
      
      posterior_trends <- function(n_samples = n_samples_t,
                                   startyear = firstYear,
                                   endyear = lastYear,
                                   region = "BCR",
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
        
      }
      
      BCR_trends <- posterior_trends()
      
      strata_trends <- posterior_trends(region = "ST_12") 
      
      
      
      
      
      
      
      
      # 
      # 
      # 
      # 
      # 
      # 
      # 
      # 
      # 
      # 
      # nyears = stan_data$nyears
      # fixedyear = stan_data$fixedyear
      # YEARS = c(min(jags_data$r_year):max(jags_data$r_year))
      # 
      # if(length(YEARS) != nyears){stop("years don't match YEARS =",length(YEARS),"nyears =",nyears)}
      # 
      # ind_fxn = function(a,b,sdn,sdob,y,fy){
      #   i = exp(a + b*(y-fy) + (0.5*(sdn^2))+ (0.5*(sdob^2)))
      #   return(i)
      # }
      # 
      # ### this could be simplified to just estimate the start and end-years
      # i_samples = NULL
      # for(yr in 1:nyears){
      #   i_t = ab_samples %>% 
      #     mutate(i = ind_fxn(alpha,beta,sdnoise,sdobs,yr,fixedyear),
      #            y = yr,
      #            year = YEARS[yr])
      #   i_samples <- bind_rows(i_samples,i_t)
      # }
      # 
      # ### this could be tweaked to sum across all routes in the original strata
      # ### just join to the strata-route dataframe - route_map
      # ### then add the non-zero-weights for the strata
      # ### then add the area-weights for the strata
      # ### and change the group-by value
      # indices = i_samples %>% group_by(s,y,year) %>% 
      #   summarise(index = mean(i),
      #             lci_i = quantile(i,LC),
      #             uci_i = quantile(i,UC),
      #             sd_i = sd(i),
      #             .groups = "keep")
      # 
      # raw = data.frame(s = stan_data$route,
      #                  y = stan_data$year,
      #                  count = stan_data$count,
      #                  obs = stan_data$observer)
      # indices = left_join(indices,raw,by = c("y","s"))
      # 
      # rts = route_map %>% tibble() %>% 
      #   select(route,site,strat) %>% 
      #   mutate(s = site) 
      # 
      # 
      # indices = left_join(indices,rts,by = "s")
      # indices$obs <- factor(indices$obs)
      # nroutes = stan_data$nroutes
      # 
      # # setting up the plot dimensions
      # npg = ceiling(nroutes/9)
      # ncl = 3
      # nrw = 3
      # if(npg*9-nroutes < 3){
      #   nrw = 2
      #   npg = ceiling(nroutes/6) 
      #   if(npg*6-nroutes < 3){
      #     ncl = 2
      #     npg = ceiling(nroutes/4)  
      #   }
      # }
      # #### 
      # pdf(paste0("trajectories/",species,"_route_trajectories2.pdf"),
      #     width = 11,
      #     height = 8.5)
      # 
      # for(j in 1:npg){
      #   traj = ggplot(data = indices,aes(x = year,y = index,colour = strat))+
      #     geom_ribbon(aes(ymin = lci_i,ymax = uci_i),alpha = 0.4,fill = grey(0.5))+
      #     geom_line()+
      #     geom_point(aes(x = year,y = count, colour = obs), fill = grey(0.5),alpha = 0.5,inherit.aes = FALSE)+
      #     facet_wrap_paginate(~ strat+route,scales = "free",ncol = ncl,nrow = nrw,page = j)+
      #     theme(legend.position = "none")
      #   try(print(traj),silent = TRUE)
      # }
      # dev.off()
      # 
      
      
    }
    
    
    ind_samples <- posterior_samples(sl_rstan,
                                        "indices",
                                     dims = c("s","y")) %>% 
      rename(site = s)
    
    I_samples <- posterior_samples(sl_rstan,
                                   "I",
                                   dims = c("y"))
    
    BETA <- posterior_samples(sl_rstan,
                                      "BETA") %>% 
      summarise(trend = mean((exp(.value)-1)*100),
                lci = quantile((exp(.value)-1)*100,LC),
                uci = quantile((exp(.value)-1)*100,UC),
                sd = sd((exp(.value)-1)*100))
    
    I_all <- I_samples %>% group_by(y) %>% 
      summarise(index = mean(.value),
                lci = quantile(.value,LC),
                uci = quantile(.value,UC),
                sd = sd(.value),
                .groups = "keep") %>% 
      mutate(scale = "SurveyWide",
             YEAR = y+(min(data_all$YEAR)-1))
    
    overall_traj <- ggplot(data = I_all,aes(x = YEAR,y = index))+
      geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.1)+
      geom_line()+
      scale_y_continuous(limits = c(0,NA))
    print(overall_traj)
    
    inds_all <- ind_samples %>% left_join(.,site_list,by = "site") %>% 
      group_by(site, site_orig,survey,y) %>% 
      summarise(index = mean(.value),
                lci = quantile(.value,LC),
                uci = quantile(.value,UC),
                sd = sd(.value),
                .groups = "keep")%>% 
      mutate(scale = "Site",
             YEAR = y+(min(data_all$YEAR)-1)) 
    
    
    inds_survey <- ind_samples %>% left_join(.,site_list,by = "site") %>% 
      group_by(survey,y,.draw) %>% 
      summarise(v = mean(.value),
                .groups = "drop") %>% 
      group_by(survey,y) %>% 
      summarise(index = mean(v),
                lci = quantile(v,LC),
                uci = quantile(v,UC),
                sd = sd(v),
                .groups = "keep")%>% 
      mutate(scale = "Survey",
             YEAR = y+(min(data_all$YEAR)-1),
             Survey = ifelse(survey == 1,"BBS","GWWA")) %>% 
      filter(!(Survey == "BBS" & YEAR >2019))
    
    
    survey_trajs <- ggplot(data = inds_survey,aes(x = YEAR,y = index))+
      geom_ribbon(aes(ymin = lci,ymax = uci),alpha = 0.1)+
      geom_line()+
      scale_y_continuous(limits = c(0,NA))+
      facet_wrap(~Survey)
    print(survey_trajs)
    # connect trends to original route names ----------------------------------
    
    route_map_out = left_join(site_centres,slops_int,by = "site")
    route_map_out$species <- species
    
    trends_out <- bind_rows(trends_out,route_map_out)
    
    
    
    route_map_out_rand = left_join(site_centres,slops_rand_int,by = "site")
    route_map_out_rand$species <- species
    
    trends_out_rand <- bind_rows(trends_out_rand,route_map_out_rand)
    
    # slopes_rand_full_int
    # route_map_out_rand = left_join(site_centres,slopes_rand_full_int,by = "site")
    # route_map_out_rand$species <- species
    # 
    # trends_out_rand <- bind_rows(trends_out_rand,route_map_out_rand)
    
    
    
    route_map_out_space = left_join(site_centres,slops_space_int,by = "site")
    route_map_out_space$species <- species
    
    trends_out_space <- bind_rows(trends_out_space,route_map_out_space)
    
    
    
    ### setting up boundaries for plots
    # load(paste0("route_maps/",species_f,"_route_data.RData"))
    
    site_bounds <- st_union(site_centres) #union to provide a simple border of the realised strata
    bb = st_bbox(site_bounds)
    xlms = as.numeric(c(bb$xmin,bb$xmax))
    ylms = as.numeric(c(bb$ymin,bb$ymax))
    
    
    
    
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
      
      route_map_out_space$Tplot <- cut(route_map_out_space$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
      route_map_out_rand$Tplot <- cut(route_map_out_rand$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
      
      
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
      geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
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
    
    
    png(filename = paste0("Figures/images/",species_f,"_Trends_",firstYear,".png"),
        res = 600,
        width = 20,
        height = 15,
        units = "cm")
    print(tmap)
    dev.off()
    
    # 
    # tmap2 = ggplot(route_map_out)+
    #   geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
    #   geom_sf(aes(colour = Tplot,size = abund))+
    #   scale_colour_manual(values = map_palette, aesthetics = c("colour"),
    #                       guide = guide_legend(reverse=TRUE),
    #                       name = paste0(lgnd_head,firstYear,"-",lastYear))+
    #   coord_sf(xlim = xlms,ylim = ylms)+
    #   theme(legend.position = "none")+
    #   labs(title = paste(species))
    # 
    # maps2[[jj]] <- tmap2
    # 
    # 
    # 
    # tmap3 = ggplot(route_map_out)+
    #   geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
    #   geom_sf(aes(colour = Tplot,size = 1/h_ci))+
    #   scale_size_continuous(range = c(0.5,3))+
    #   scale_colour_manual(values = map_palette, aesthetics = c("colour"),
    #                       guide = guide_legend(reverse=TRUE),
    #                       name = paste0(lgnd_head,firstYear,"-",lastYear))+
    #   coord_sf(xlim = xlms,ylim = ylms)+
    #   labs(title = paste(species))
    # 
    # maps3[[jj]] <- tmap3
    # 
    # 
    # 
    # 
    # tmap_space = ggplot(route_map_out_space)+
    #   geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
    #   geom_sf(aes(colour = Tplot,size = abund))+
    #   scale_colour_manual(values = map_palette, aesthetics = c("colour"),
    #                       guide = guide_legend(reverse=TRUE),
    #                       name = paste0(lgnd_head,firstYear,"-",lastYear))+
    #   coord_sf(xlim = xlms,ylim = ylms)+
    #   theme(legend.position = "none")+
    #   labs(title = paste("spatial component"))
    # maps_space[[jj]] <- tmap_space
    # 
    # 
    # 
    # 
    # 
    # tmap_rand = ggplot(route_map_out_rand)+
    #   geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
    #   geom_sf(aes(colour = Tplot,size = abund))+
    #   scale_colour_manual(values = map_palette, aesthetics = c("colour"),
    #                       guide = guide_legend(reverse=TRUE),
    #                       name = paste0(lgnd_head,firstYear,"-",lastYear))+
    #   coord_sf(xlim = xlms,ylim = ylms)+
    #   theme(legend.position = "none")+
    #   labs(title = paste("random component"))
    # 
    # maps_rand[[jj]] <- tmap_rand
    # 
    
    print(species)
    # write.csv(route_map_out,
    #           file = paste0("output/",species," ",firstYear," ",lastYear,"_Canadian_trends_and_intercepts.csv"))
    
    
  
}



# overall trend maps and trends -------------------------------------------


pdf(file = paste0("figures/Combined_",firstYear,"_",lastYear,"_",scope,"_trend_map_route2.pdf"),
    height = 8.5,
    width = 11)
for(j in 1:length(maps)){
  if(!is.null(maps[[j]])){print(maps[[j]])}
}
dev.off()



# comparison trend maps and trends -------------------------------------------

library(patchwork)

pdf(file = paste0("figures/Combined_",firstYear,"_",lastYear,"_",scope,"_trend_map_route2_by_half_CI.pdf"),
    height = 8.5,
    width = 11)
for(j in 1:length(maps3)){
  if(!is.null(maps3[[j]])){print(maps3[[j]])}
}
dev.off()





pdf(file = paste0("figures/Combined_",firstYear,"_",lastYear,"_",scope,"_space_all_trend_map_route2.pdf"),
    width = 8.5,
    height = 11)
for(j in 1:length(maps)){
  if(!is.null(maps[[j]])){
    print(maps2[[j]] /(maps_space[[j]]))
  }
}
dev.off()






# add in the original BBS route database start coordinates

starts = unique(strat_data$route_strat[,c("rt.uni","Latitude","Longitude")])

# rename the trend output columns -----------------------------------------

trends_out2 = trends_out %>% 
  data.frame() %>% 
  mutate(h_ci = (uci_trend-lci_trend)/2) %>% 
  select(.,
         species,route,strat,
         trend,lci_trend,uci_trend,h_ci,
         abund,lci_i,uci_i,
         b,lci,uci,sd) %>% 
  relocate(.,
           species,route,strat,
           trend,lci_trend,uci_trend,h_ci,
           abund,lci_i,uci_i,
           b,lci,uci,sd) %>% 
  rename(.,
         english_name = species,BBS_route = route, BBS_stratum = strat,
         Trend = trend,lci95_Trend = lci_trend,uci95_Trend = uci_trend,half_CI_width = h_ci,
         Mean_abundance = abund,lci95_Mean_abundance = lci_i,uci95_Mean_abundance = uci_i,
         slope = b,lci95_slope = lci,uci95_slope = uci,sd_slope = sd) %>% 
  left_join(starts,by = c("BBS_route" = "rt.uni"))


trends_out_space2 = trends_out_space %>% 
  data.frame() %>% 
  mutate(h_ci = (uci_trend-lci_trend)/2) %>% 
  select(.,
         species,route,strat,
         trend,lci_trend,uci_trend,h_ci,
         abund,lci_i,uci_i,
         b,lci,uci,sd) %>% 
  relocate(.,
           species,route,strat,
           trend,lci_trend,uci_trend,h_ci,
           abund,lci_i,uci_i,
           b,lci,uci,sd) %>% 
  rename(.,
         english_name = species,BBS_route = route, BBS_stratum = strat,
         Trend = trend,lci95_Trend = lci_trend,uci95_Trend = uci_trend,half_CI_width = h_ci,
         Mean_abundance = abund,lci95_Mean_abundance = lci_i,uci95_Mean_abundance = uci_i,
         slope = b,lci95_slope = lci,uci95_slope = uci,sd_slope = sd) %>% 
  left_join(starts,by = c("BBS_route" = "rt.uni"))

trends_out_rand2 = trends_out_rand %>% 
  data.frame() %>% 
  mutate(h_ci = (uci_trend-lci_trend)/2) %>% 
  select(.,
         species,route,strat,
         trend,lci_trend,uci_trend,h_ci,
         abund,lci_i,uci_i,
         b,lci,uci,sd) %>% 
  relocate(.,
           species,route,strat,
           trend,lci_trend,uci_trend,h_ci,
           abund,lci_i,uci_i,
           b,lci,uci,sd) %>% 
  rename(.,
         english_name = species,BBS_route = route, BBS_stratum = strat,
         Trend = trend,lci95_Trend = lci_trend,uci95_Trend = uci_trend,half_CI_width = h_ci,
         Mean_abundance = abund,lci95_Mean_abundance = lci_i,uci95_Mean_abundance = uci_i,
         slope = b,lci95_slope = lci,uci95_slope = uci,sd_slope = sd) %>% 
  left_join(starts,by = c("BBS_route" = "rt.uni"))


# Export the trend estimates ----------------------------------------------


write.csv(trends_out2,
          file = paste0("output/combined_",firstYear,"_",lastYear,"_",scope,"_trends_and_intercepts.csv"))


write.csv(trends_out_space2,
          file = paste0("output/combined_",firstYear,"_",lastYear,"_",scope,"_spatial_trends_and_intercepts2.csv"))



write.csv(trends_out_rand2,
          file = paste0("output/combined_",firstYear,"_",lastYear,"_",scope,"_random_trends_and_intercepts2.csv"))



# graph the spatial and random variance comparison ------------------------

var_plot = ggplot()+
  geom_boxplot(data = sdbeta_space_rand,aes(x = .variable,y = mean))+
  theme_classic()

sdbeta_difs <- sdbeta_dif %>% 
  mutate(species = fct_reorder(species,mean))

var_dif_plot = ggplot(data = sdbeta_difs,aes(x = species,y = mean))+
  geom_point(aes(size = 1/(sd^2)))+
  geom_errorbar(aes(ymin = lci, ymax = uci),alpha = 0.3,width = 0)+
  geom_hline(yintercept = 0)+
  scale_size_continuous(range = c(0.5,2))+
  labs(title = "Difference in sd_beta (random - spatial)")+
  theme(axis.text = element_text(size = 5),
        legend.position = "none")+
  coord_flip()
pdf(file = paste0("figures/Combined_",firstYear,"_",lastYear,"_",scope,"_difference_beta_sd.pdf"),
    width = 8.5,
    height = 17)
print(var_dif_plot)
dev.off()



