

library(bbsBayes2)
library(tidyverse)
library(patchwork)
library(sf)
library(cmdstanr)

species <- "Golden-winged Warbler"
species_f <- gsub(species,pattern = " ",replacement = "_") # species name without spaces

model = "first_diff"

stratification <- "latlong"

model_variant <- "spatial"


# Fit BBS only model ------------------------------------------------------
s <- stratify(stratification,
              species,
              release = 2022) # ensures only data up to 2021 are included matches GWWA data

p <- prepare_data(s,
                  min_n_routes = 1,
                  min_max_route_years = 2)

strata_map <- load_map(stratification)

ps <- prepare_spatial(p, strata_map)



pm <- prepare_model(ps, model, model_variant)


bbs_strata_map <- strata_map %>% 
  filter(strata_name %in% pm$meta_strata$strata_name)



bbs_summary_pre_1990 <- pm$raw_data  %>% 
  filter(year < 1991) %>% 
  group_by(strata_name) %>% 
  summarise(log_mean_count = log(mean(count)),
            .groups = "keep") %>% 
  mutate(time_period = "Pre-1990")


bbs_summary_post_2000 <- pm$raw_data  %>% 
  filter(year > 2000) %>% 
  group_by(strata_name) %>% 
  summarise(log_mean_count = log(mean(count)),
            .groups = "keep") %>% 
  mutate(time_period = "Post-2000")


bbs_summary <- bind_rows(bbs_summary_pre_1990,
                         bbs_summary_post_2000) %>% 
  mutate(survey = "BBS")


fit_bbs <- run_model(pm,
                     refresh = 500,
                     iter_warmup = 2000,
                     iter_sampling = 2000,
                     max_treedepth = 11,
                     output_basename = "bbs_only",
                     output_dir = "output")

fit_bbs <- readRDS("output/bbs_only.rds")




inds <- bbsBayes2::generate_indices(fit_bbs)

trends <- bbsBayes2::generate_trends(inds)
trends_09 <- bbsBayes2::generate_trends(inds,min_year = 2009)
trends_11 <- bbsBayes2::generate_trends(inds,min_year = 2011)

map_t <- bbsBayes2::plot_map(trends)
map_t_11 <- bbsBayes2::plot_map(trends_11)
map_t_09 <- bbsBayes2::plot_map(trends_09)

maps <- map_t / map_t_05 

pdf(paste0("figures/bbs_long_term_trends.pdf"))
print(map_t)
dev.off()
pdf(paste0("figures/bbs_2011_term_trends.pdf"))
print(map_t_11)
dev.off()
pdf(paste0("figures/bbs_2009_term_trends.pdf"))
print(map_t_09)
dev.off()



bcrs <- bbsBayes2::load_map("bcr") %>% 
  rename(bcr = strata_name)

bcr_latlong <- sf::st_join(x = strata_map,
                       y = bcrs,
                       largest = TRUE) %>% 
  sf::st_drop_geometry()

inds_bcr <- bbsBayes2::generate_indices(fit_bbs,
                             regions_index = bcr_latlong,
                             regions = "bcr")

trends_bcr <- bbsBayes2::generate_trends(inds_bcr)
trends_bcr_2011 <- bbsBayes2::generate_trends(inds_bcr,
                                   min_year = 2011)
trends_bcr_2009 <- bbsBayes2::generate_trends(inds_bcr,
                              min_year = 2009)

trajs_bcr <- bbsBayes2::plot_indices(inds_bcr,
                          min_year = 2005)
trajs_bcr[["BCR28"]]

trends_bcr_2009$trends[6,c("trend","percent_change")]

# trend percent_change
# <dbl>          <dbl>
#   1 -6.07          -52.8






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
gwwa_map <- st_transform(gwwa_map,crs = bbsBayes2::equal_area_crs)

ttp <- ggplot()+
  geom_sf(data = gwwa_map,aes(colour = Year))

print(ttp)

quads <- st_read(dsn = "data",
                 layer = "delorme_conus") 
st_crs(quads) <- 4326 #WGS 84

quads <- st_transform(quads,crs = bbsBayes2::equal_area_crs)

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

# aggregate gwwa data to merge with the BBS ----------------

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

site_centres_gwwa <- sf::st_join(site_centres_gwwa,bbs_strata_map,
                                 join = st_nearest_feature) %>%
  select(-area_sq_km) %>% 
  left_join(pm$meta_strata, by = 'strata_name')



gwwa_data <- site_centres_gwwa %>% 
  as.data.frame() %>% 
  select(site_orig,site_gwwa,strata_name,strata) %>% 
  right_join(.,gwwa_ag,by = c("site_orig" = "page")) %>% 
  rename(YEAR = Year) %>% 
  mutate(offset = log(n_survey)) 



gwwa_summary <- gwwa_data  %>% 
  group_by(strata_name) %>% 
  summarise(log_mean_count = log(mean(count)),
            .groups = "keep") %>% 
  mutate(time_period = "Post-2000",
         survey = "GWWA")

datasets_summary <- bind_rows(bbs_summary,
                              gwwa_summary) %>% 
  mutate(time_period = factor(time_period,
                              levels = c("Pre-1990","Post-2000"),
                              ordered = TRUE))

map_summary <- bbs_strata_map %>% 
  left_join(datasets_summary,
            by = c("strata_name"))

base_state_map <- bbsBayes::load_map(stratify_by = "state")

data_bounding_box <- map_summary %>% 
  sf::st_bbox()


map_gwwa_data <- ggplot()+
  geom_sf(data = base_state_map,
          fill = NA,
          colour = grey(0.8))+
  geom_sf(data = map_summary,
          aes(fill = log_mean_count))+
  coord_sf(xlim = data_bounding_box[c("xmin","xmax")],
           ylim = data_bounding_box[c("ymin","ymax")])+
  facet_grid(cols = vars(time_period),rows = vars(survey))+
  theme_bw() +
  scale_fill_viridis_c(na.value = "grey70")

print(map_gwwa_data)

pdf(file = "Figures/Survey_mean_counts_comparison.pdf",
    width = 8, height = 7)
print(map_gwwa_data)
dev.off()

map_gwwa_data <- ggplot()+
  geom_sf(data = base_state_map,
          fill = NA,
          colour = grey(0.5))+
  geom_sf(data = map_summary,
          aes(fill = log_mean_count))+
  coord_sf(xlim = data_bounding_box[c("xmin","xmax")],
           ylim = data_bounding_box[c("ymin","ymax")])+
  facet_grid(cols = vars(time_period),rows = vars(survey))+
  theme_bw() +
  scale_fill_viridis_c(na.value = "grey70")

png(filename = "Mean_counts_by_survey_and_period.png",
    res = 300,
    width = 8, height = 7, units = "in")
print(map_gwwa_data)
dev.off()



# Compile stan data for gwwa model ----------------------------------------


stan_data <- pm$model_data

stan_data[["calc_cv"]] <- 0
stan_data[["n_train"]] <- stan_data[["n_counts"]]
stan_data[["n_test"]] <- as.integer(1)
stan_data[["train"]] <- 1:stan_data[["n_counts"]]
stan_data[["test"]] <- as.integer(1)


  
  
stan_data[["n_sites_gwwa"]] <- max(gwwa_data$site_gwwa)
stan_data[["n_counts_gwwa"]] <- nrow(gwwa_data)
stan_data[["n_observers_gwwa"]] <- max(gwwa_data$observer_gwwa)
stan_data[["base_year_gwwa"]] <- min(gwwa_data$YEAR-(min(pm$raw_data$year)-1))-1
stan_data[["n_years_gwwa"]] <- (stan_data[["n_years"]]-stan_data[["base_year_gwwa"]])

stan_data[["count_gwwa"]] <- gwwa_data$count
stan_data[["strat_gwwa"]] <- gwwa_data$strata
stan_data[["year_gwwa"]] <- gwwa_data$YEAR-(min(pm$raw_data$year)-1)
stan_data[["site_gwwa"]] <- gwwa_data$site_gwwa
stan_data[["observer_gwwa"]] <- gwwa_data$observer_gwwa
stan_data[["off_set"]] <- gwwa_data$offset



mod.file = "models/first_difference_spatial_gwwa_bbs_CV.stan"
out_base <- paste(species_f,"first_difference_Spatial","bbs_gwwa",sep = "_") # text string to identify the saved output from the Stan process unique to species and model, but probably something the user wants to control


## compiles Stan model (this is only necessary if the model has been changed since it was last run on this machine)
model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))

output_dir <- "output/" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output


# Fit integrated model ----------------------------------------------------


stanfit <- model$sample(
  data=stan_data,
  refresh=200,
  chains=4, 
  iter_sampling=4000,
  iter_warmup=2000,
  parallel_chains = 4,
  thin = 2,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 12,
  #seed = 123,
  #init = init_def,
  output_dir = output_dir,
  output_basename = out_base,
  show_exceptions = FALSE)

# shinystan::launch_shinystan(shinystan::as.shinystan(stanfit))


# loo_out <- stanfit$loo()


fit_summary <- stanfit$summary()

fit_integrated <- pm
fit_integrated[["folds"]] <- NULL
fit_integrated[["init_values"]] <- NULL

fit_integrated[["model_fit"]] <- stanfit

saveRDS(fit_integrated,
        paste0(output_dir,"/integrated.rds"))

save(list = c("stanfit","stan_data",
              "out_base",
              "fit_summary"),
     file = paste0(output_dir,"/",out_base,"_Stan_fit.RData"))






inds <- bbsBayes2::generate_indices(fit_integrated)

trends <- bbsBayes2::generate_trends(inds)
trends_09 <- bbsBayes2::generate_trends(inds,min_year = 2009)
trends_11 <- bbsBayes2::generate_trends(inds,min_year = 2011)

map_t <- bbsBayes2::plot_map(trends)
map_t_05 <- bbsBayes2::plot_map(trends_05)
map_t_11 <- bbsBayes2::plot_map(trends_11)
map_t_09 <- bbsBayes2::plot_map(trends_09)

maps <- map_t / map_t_05 

pdf(paste0("figures/bbs_long_term_trends.pdf"))
print(map_t)
dev.off()
pdf(paste0("figures/bbs_2005_term_trends.pdf"))
print(map_t_05)
dev.off()
pdf(paste0("figures/bbs_2009_term_trends.pdf"))
print(map_t_09)
dev.off()

pdf(paste0("figures/bbs_2011_term_trends.pdf"))
print(map_t_11)
dev.off()


bcrs <- bbsBayes2::load_map("bcr") %>% 
  rename(bcr = strata_name)

bcr_latlong <- sf::st_join(x = strata_map,
                           y = bcrs,
                           largest = TRUE) %>% 
  sf::st_drop_geometry()

inds_bcr <- bbsBayes2::generate_indices(fit_bbs,
                                        regions_index = bcr_latlong,
                                        regions = "bcr")

trends_bcr <- bbsBayes2::generate_trends(inds_bcr)
trends_bcr_2005 <- bbsBayes2::generate_trends(inds_bcr,
                                              min_year = 2005)
trends_bcr_2009 <- bbsBayes2::generate_trends(inds_bcr,
                                              min_year = 2009)

trajs_bcr <- bbsBayes2::plot_indices(inds_bcr,
                                     min_year = 2005)
trajs_bcr[["BCR28"]]

trends_bcr_2009$trends[6,c("trend","percent_change")]

# trend percent_change
# <dbl>          <dbl>
#   1 -6.07          -52.8





