library(bbsBayes)
library(tidyverse)
library(cmdstanr)

bbs_data <- stratify(by = "latlong")


#setwd("C:/GitHub/bbsStanBayes")

#setwd("C:/Users/SmithAC/Documents/GitHub/bbsStanBayes")


model_sel <- "firstdiff"

species <- "Golden-winged Warbler"
species_f <- gsub(species,pattern = " ",replacement = "_") # species name without spaces



# Stan models -------------------------------------------------------------

## spatial versions of both teh slope and gamye exist for the Stan models and can be fit with this script 
fit_spatial <- TRUE # TRUE = spatial sharing of information and FALSE = non-spatial sharing

## the bbsBayes prepare_data function doesn't create all of the objects required for the Stan versions of the models
## this source() call over-writes the bbsBayes function prepare_data()
source("Functions/prepare-data-Stan.R")
if(fit_spatial){
source("Functions/neighbours_define.R") # function to generate spatial neighbourhoods to add to the spatial applications of the models
}



sp_data <- prepare_data(bbs_data,
                        species_to_run = species,
                        model = model_sel,
                        min_max_route_years = 2,
                        min_n_routes = 1)


stan_data <- sp_data

# Spatial neighbourhoods --------------------------------------------------
if(fit_spatial){
  
base_strata_map <- bbsBayes::load_map(stratify_by = stan_data[["stratify_by"]])

alt_df <- stan_data[["alt_data"]][[1]]
strata_df <- alt_df %>% 
  select(strat,strat_name) %>% 
  distinct() %>% 
  arrange(strat)

realized_strata_map <- base_strata_map %>% 
  inner_join(.,strata_df,by = c("ST_12" = "strat_name"))


neighbours <- neighbours_define(real_strata_map = realized_strata_map, #sf map of strata
                                strat_link_fill = 10000, #distance to fill if strata are not connected
                                buffer = TRUE,
                                convex_hull = FALSE,
                                plot_neighbours = TRUE,
                                species = species,
                                plot_dir = "maps/",
                                plot_file = "_strata_map",
                                save_plot_data = TRUE,
                                voronoi = FALSE,
                                nn_fill = FALSE,
                                add_map = NULL,
                                strat_indicator = "strat",
                                island_link_dist_factor = 1.2 #consider nearest strata neighbours if distances are within this factor of each other, when linking otherwise isolated islands of strata
                                )



stan_data[["N_edges"]] = neighbours$N_edges
stan_data[["node1"]] <- neighbours$node1
stan_data[["node2"]] <- neighbours$node2
}#end of if fit_spatial



# extra list elements not required by Stan
# temporarily save them as objects then return them to the stan_data list after model fitting (below)
tmp_stratify_by <- stan_data[["stratify_by"]]  
tmp_model <- stan_data[["model"]]
tmp_alt_data <- stan_data[["alt_data"]]


stan_data[["stratify_by"]] <- NULL 
stan_data[["model"]] <- NULL
stan_data[["alt_data"]] <- NULL

stan_data[["zero_betas"]] <- rep(0,stan_data$nstrata)
stan_data[["Iy1"]] <- c((stan_data$fixedyear-1):1)
stan_data[["nIy1"]] <- length(stan_data[["Iy1"]])

stan_data[["Iy2"]] <- c((stan_data$fixedyear+1):stan_data$nyears)
stan_data[["nIy2"]] <- length(stan_data[["Iy2"]])



# Add in gwwa data --------------------------------------------------------



gwwa_dat <- read.csv("data/2009_2021_GWWA_Data_final.csv")

gwwa_dat <- gwwa_dat %>% 
  mutate(State = str_trim(str_to_upper(State)),
         GWWA = as.integer(GWWA))  %>% 
  filter(!is.na(GWWA)) %>% 
  mutate(Observer = gsub(pattern = "[[:punct:]]|[[:space:]]",
                         replacement = "",
                         x = Observer))


## make spatial point layer of gwwa monitoring sites -----------------------

gwwa_map = st_as_sf(gwwa_dat,coords = c("Longitude","Latitude"))
st_crs(gwwa_map) <- 4326 #WGS 84??? - this is a guess
#load strata map
gwwa_map <- st_transform(gwwa_map,crs = st_crs(realized_strata_map))

ttp <- ggplot()+
  geom_sf(data = gwwa_map,aes(colour = Year))

print(ttp)

quads <- st_read(dsn = "data",
                 layer = "delorme_conus") 
st_crs(quads) <- 4326 #WGS 84

quads <- st_transform(quads,crs = st_crs(realized_strata_map))

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
            .groups = "drop") %>% 
  mutate(observer_gwwa = as.integer(factor(Observer)))



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



# spatial join gwwa quad centres with latlong BBS strata ------------------
# joins them to the nearest strata because there are 3 gwwa quads that are
# just outside of the strata that have BBS data
# we don't want to estimate trajectories for strata with only GWWA
# so link them to their nearest degreee block

site_centres_gwwa <- sf::st_join(site_centres_gwwa,realized_strata_map,
                                 join = st_nearest_feature)



gwwa_data <- site_centres_gwwa %>% 
  as.data.frame() %>% 
  select(site_orig,site_gwwa,ST_12,strat) %>% 
  right_join(.,gwwa_ag,by = c("site_orig" = "page")) %>% 
  rename(YEAR = Year,
         strat_name = ST_12) %>% 
  mutate(offset = log(n_survey)) 


# tmp <- gwwa_data %>% 
#   select(YEAR,site_gwwa,offset) %>% 
#   distinct() %>% 
#   mutate(first_years_gwwa = YEAR-(min(tmp_alt_data$full_data_frame$rYear)-1))

# first_years_gwwa_df <- gwwa_data %>% 
#   select(site_gwwa,strat,YEAR) %>% 
#   distinct() %>% 
#   group_by(site_gwwa,strat) %>% 
#   summarise(first_years_gwwa = min(YEAR)-(min(tmp_alt_data$full_data_frame$rYear)-1)) %>% 
#   arrange(site_gwwa) 


#strat_vec_gwwa <- unique(gwwa_data$strat)

stan_data[["nsites_gwwa"]] <- max(gwwa_data$site_gwwa)
stan_data[["ncounts_gwwa"]] <- nrow(gwwa_data)
stan_data[["nobservers_gwwa"]] <- max(gwwa_data$observer_gwwa)
#stan_data[["nstrat_gwwa"]] <- length(strat_vec_gwwa)
stan_data[["base_year_gwwa"]] <- min(gwwa_data$YEAR-(min(tmp_alt_data$full_data_frame$rYear)-1))-1
stan_data[["nyears_gwwa"]] <- (stan_data[["nyears"]]-stan_data[["base_year_gwwa"]])

stan_data[["count_gwwa"]] <- gwwa_data$count
stan_data[["strat_gwwa"]] <- gwwa_data$strat
stan_data[["year_gwwa"]] <- gwwa_data$YEAR-(min(tmp_alt_data$full_data_frame$rYear)-1)
stan_data[["site_gwwa"]] <- gwwa_data$site_gwwa
stan_data[["observer_gwwa"]] <- gwwa_data$observer_gwwa
stan_data[["off_set"]] <- gwwa_data$offset












mod.file = "models/first_difference_spatial_gwwa_bbs_CV.stan"
out_base <- paste(species_f,"first_difference_Spatial","bbs_gwwa",sep = "_") # text string to identify the saved output from the Stan process unique to species and model, but probably something the user wants to control


## compiles Stan model (this is only necessary if the model has been changed since it was last run on this machine)
model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))

output_dir <- "output/" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output

### this init_def is something that the JAGS versions don't need. It's a function definition, so perhaps something we'll have to build
### into the fit_model function
### the slightly annoying thing is that it's model-specific, so we'll need to have a few versions of it
init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts*stan_data$use_pois,0,0.1),
                             strata_raw = rnorm(stan_data$nstrata,0,0.1),
                             STRATA = 0,
                             nu = 10,
                             sdstrata = runif(1,0.01,0.1),
                             eta = 0,
                             obs_raw = rnorm(stan_data$nobservers,0,0.1),
                             ste_raw = rnorm(stan_data$nsites,0,0.1),
                             sdnoise = runif(1,0.3,1.3),
                             sdobs = runif(1,0.01,0.1),
                             sdste = runif(1,0.01,0.2),
                             sdbeta = runif(1,0.01,0.1),
                             sdBETA = runif(1,0.01,0.1),
                             BETA_raw = rnorm((stan_data$nyears-1),0,0.1),
                             beta_raw = matrix(rnorm(stan_data$nstrata*(stan_data$nyears-1),0,0.1),nrow = stan_data$nstrata,ncol = stan_data$nyears-1),
                             STE_gwwa = rnorm(stan_data$nsites_gwwa,0,0.1),
                             noise_gwwa_raw = rnorm(stan_data$ncounts*stan_data$use_pois,0,0.1),
                             nu_gwwa = 10,
                             obs_gwwa_raw = rnorm(stan_data$nobservers_gwwa,0,0.1),
                             ste_gwwa_raw = rnorm(stan_data$nsites_gwwa,0,0.1),
                             sdnoise_gwwa = runif(1,0.01,0.1),
                             sdobs_gwwa = runif(1,0.01,0.1),
                             sdste_gwwa = runif(1,0.01,0.1))}



stanfit <- model$sample(
  data=stan_data,
  refresh=200,
  chains=4, 
  iter_sampling=1000,
  iter_warmup=2000,
  parallel_chains = 4,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 12,
  #seed = 123,
  #init = init_def,
  output_dir = output_dir,
  output_basename = out_base)

# shinystan::launch_shinystan(shinystan::as.shinystan(stanfit))


# loo_out <- stanfit$loo()


fit_summary <- stanfit$summary()

stan_data[["stratify_by"]] <- tmp_stratify_by 
stan_data[["model"]] <- tmp_model
stan_data[["alt_data"]] <- tmp_alt_data
stan_data[["strat_name"]] <- tmp_alt_data$strat_name

save(list = c("stanfit","stan_data",
              "out_base",
              "fit_summary"),
     file = paste0(output_dir,"/",out_base,"_Stan_fit.RData"))





