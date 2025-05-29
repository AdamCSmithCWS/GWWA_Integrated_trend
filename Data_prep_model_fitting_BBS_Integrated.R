

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
  summarise(log_mean_count = log(mean(count,na.rm = TRUE)),
            mean_count = mean(count,na.rm = TRUE),
            .groups = "keep") %>% 
  mutate(time_period = "Pre-1990")


bbs_summary_post_2000 <- pm$raw_data  %>% 
  filter(year > 2000) %>% 
  group_by(strata_name) %>% 
  summarise(log_mean_count = log(mean(count,na.rm = TRUE)),
            mean_count = mean(count,na.rm = TRUE),
            .groups = "keep") %>% 
  mutate(time_period = "Post-2000")


bbs_summary <- bind_rows(bbs_summary_pre_1990,
                         bbs_summary_post_2000) %>% 
  mutate(survey = "BBS")

re_fit_bbs <- FALSE
if(re_fit_bbs){
fit_bbs <- run_model(pm,
                     refresh = 500,
                     iter_warmup = 2000,
                     iter_sampling = 2000,
                     max_treedepth = 11,
                     output_basename = "bbs_only",
                     output_dir = "output")
}else{
fit_bbs <- readRDS("output/bbs_only.rds")
}



inds <- bbsBayes2::generate_indices(fit_bbs,
                                    hpdi = TRUE)

trends <- bbsBayes2::generate_trends(inds,
                                     quantiles = c(0.025,0.1,0.9,0.975))
trends_09 <- bbsBayes2::generate_trends(inds,min_year = 2009,
                                        quantiles = c(0.025,0.1,0.9,0.975))
trends_11 <- bbsBayes2::generate_trends(inds,min_year = 2011,
                                        quantiles = c(0.025,0.1,0.9,0.975))

map_t <- bbsBayes2::plot_map(trends)
map_t_11 <- bbsBayes2::plot_map(trends_11)
map_t_09 <- bbsBayes2::plot_map(trends_09)

maps <- map_t / map_t_09 

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
                             regions = "bcr",
                             hpdi = TRUE)

trends_bcr <- bbsBayes2::generate_trends(inds_bcr,
                                         quantiles = c(0.025,0.1,0.9,0.975))
trends_bcr_2011 <- bbsBayes2::generate_trends(inds_bcr,
                                   min_year = 2011,
                                   quantiles = c(0.025,0.1,0.9,0.975))
trends_bcr_2009 <- bbsBayes2::generate_trends(inds_bcr,
                              min_year = 2009,
                              quantiles = c(0.025,0.1,0.9,0.975))

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
  summarise(log_mean_count = log(mean(count,na.rm = TRUE)),
            mean_count = mean(count,na.rm = TRUE),
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
          aes(fill = mean_count))+
  coord_sf(xlim = data_bounding_box[c("xmin","xmax")],
           ylim = data_bounding_box[c("ymin","ymax")])+
  facet_grid(cols = vars(time_period),rows = vars(survey))+
  theme_bw() +
  scale_fill_viridis_c(na.value = "grey70",
                       transform = "log10",
                       name = "Observed\nmean count\n")+
  labs(caption = "Observed mean counts during BBS (top row) and species specific 
       surveys for Golden-winged Warblers (GWWA, bottom row). During the 
       first few decades of the BBS (pre-1990), the species was regularly 
       observed in most of its range. Since 2000, the species is very 
       rarely observed during BBS and completely absent from all surveys 
       in the grey grid-cells. By contrast, the species is observed in 
       much higher numbers during the species-specific surveys",
       alt = "Observed mean counts during BBS (top row) and species specific 
       surveys for Golden-winged Warblers (GWWA, bottom row).During the 
       first few decades of the BBS (pre-1990), the species was regularly 
       observed in most of its range. Since 2000, the species is very 
       rarely observed during BBS and completely absent from all surveys 
       in the grey grid-cells. By contrast, the species is observed in 
       much higher numbers during the species-specific surveys")+
  theme(plot.caption = element_text(hjust = 0))

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
          aes(fill = mean_count))+
  coord_sf(xlim = data_bounding_box[c("xmin","xmax")],
           ylim = data_bounding_box[c("ymin","ymax")])+
  facet_grid(cols = vars(time_period),rows = vars(survey))+
  theme_bw() +
  scale_fill_viridis_c(na.value = "grey70",
                       transform = "log10",
                       name = "Observed\nmean count\n")+
  labs(caption = "Observed mean counts during BBS (top row) and species specific surveys for Golden-winged
  Warblers (GWWA, bottom row). During the first few decades of the BBS (pre-1990), the species was regularly
  observed in most of its range. Since 2000, the species is very rarely observed during BBS and
  completely absent from all surveys in the grey grid-cells. By contrast, the species is observed in
  much higher numbers during the species-specific surveys",
       alt = "Observed mean counts during BBS (top row) and species specific 
       surveys for Golden-winged Warblers (GWWA, bottom row).During the 
       first few decades of the BBS (pre-1990), the species was regularly 
       observed in most of its range. Since 2000, the species is very 
       rarely observed during BBS and completely absent from all surveys 
       in the grey grid-cells. By contrast, the species is observed in 
       much higher numbers during the species-specific surveys")+
  theme(plot.caption = element_text(hjust = 0))

png(filename = paste0("figures/images/Mean_counts_by_survey_and_period.png"),
    res = 300,
    width = 8, height = 7, units = "in")
print(map_gwwa_data)
dev.off()


map_gwwa_data <- ggplot()+
  geom_sf(data = base_state_map,
          fill = NA,
          colour = grey(0.5))+
  geom_sf(data = map_summary,
          aes(fill = mean_count))+
  coord_sf(xlim = data_bounding_box[c("xmin","xmax")],
           ylim = data_bounding_box[c("ymin","ymax")])+
  facet_grid(cols = vars(time_period),rows = vars(survey))+
  theme_bw() +
  scale_fill_viridis_c(na.value = "grey70",
                       transform = "log10",
                       name = "Observed\nmean count\n")
png(filename = paste0("Mean_counts_by_survey_and_period.png"),
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

re_fit_integrated <- FALSE

if(re_fit_integrated){

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

}else{
  load(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))
  fit_integrated <- readRDS(paste0(output_dir,"/integrated.rds"))
}



# Integrated indices and trends -------------------------------------------



inds_int <- bbsBayes2::generate_indices(fit_integrated,
                                        hpdi = TRUE)

trends_int <- bbsBayes2::generate_trends(inds_int,
                                         quantiles = c(0.025,0.1,0.9,0.975))
trends_int_09 <- bbsBayes2::generate_trends(inds_int,min_year = 2009,
                                            quantiles = c(0.025,0.1,0.9,0.975))
trends_int_11 <- bbsBayes2::generate_trends(inds_int,min_year = 2011,
                                            quantiles = c(0.025,0.1,0.9,0.975))

# map_t_int <- bbsBayes2::plot_map(trends_int)
# map_t_int_11 <- bbsBayes2::plot_map(trends_int_11)
# map_t_int_09 <- bbsBayes2::plot_map(trends_int_09)
# 
# maps <- map_t / map_t_09 
# 


bcrs <- bbsBayes2::load_map("bcr") %>% 
  rename(bcr = strata_name)

bcr_latlong <- sf::st_join(x = strata_map,
                           y = bcrs,
                           largest = TRUE) %>% 
  sf::st_drop_geometry()

inds_int_bcr <- bbsBayes2::generate_indices(fit_integrated,
                                        regions_index = bcr_latlong,
                                        regions = "bcr",
                                        hpdi = TRUE)

trends_int_bcr <- bbsBayes2::generate_trends(inds_int_bcr,
                                             quantiles = c(0.025,0.1,0.9,0.975))
trends_int_bcr_2009 <- bbsBayes2::generate_trends(inds_int_bcr,
                                              min_year = 2009,
                                              quantiles = c(0.025,0.1,0.9,0.975))

trajs_bcr_int <- bbsBayes2::plot_indices(inds_int_bcr,
                                     min_year = 2005)
trajs_bcr_int[["BCR28"]]

trajs_bcr <- bbsBayes2::plot_indices(inds_bcr,
                                         min_year = 2005)
trajs_bcr[["BCR28"]]

trends_int_bcr_2009$trends[6,c("trend","percent_change")]
trends_bcr_2009$trends[6,c("trend","percent_change")]

# trend percent_change
# <dbl>          <dbl>
#   1 -6.07          -52.8

# Plotting trajectories ---------------------------------------------------


inds_int_plot <- inds_int$indices %>% 
  mutate(survey = "Integrated")


inds_int_plot_bcr <- inds_int_bcr$indices %>% 
  mutate(survey = "Integrated")

inds_plot_bcr <- inds_bcr$indices %>% 
  mutate(survey = "BBS-only")

inds_plot <- inds$indices %>% 
  mutate(survey = "BBS-only") %>% 
  bind_rows(inds_int_plot,
            inds_int_plot_bcr,
            inds_plot_bcr) %>% 
  filter(region %in% c("continent","BCR28")) %>% 
  mutate(Region = factor(region,
                         levels = c("continent","BCR28"),
                         labels = c("Survey-wide","Appalachian Mountains"),
                         ordered = TRUE))



trajs_l <- ggplot(data = filter(inds_plot,region %in% c("continent","BCR28")))+
  geom_ribbon(aes(x = year,
                  ymin = `index_q_0.05`,
                  ymax = `index_q_0.95`,
                  fill = survey),
              alpha = 0.2)+
  geom_line(aes(x = year,
                y = index,
                colour = survey))+
  facet_grid(rows = vars(Region),
             scales = "free")+
  theme_bw()+
  labs(title = "1966-2021")+
  ylab("Expected count on BBS")+
  xlab("")+
  scale_colour_viridis_d(begin = 0.1,end = 0.7,
                         aesthetics = c("fill","colour"))+
  theme(plot.title = element_text(hjust = 1),
        axis.text = element_text(size = 8))

trajs_s <- ggplot(data = filter(inds_plot,
                                #region %in% c("continent","BCR28"),
                                region %in% c("BCR28"),
                                year > 2005))+
  geom_ribbon(aes(x = year,
                  ymin = `index_q_0.05`,
                  ymax = `index_q_0.95`,
                  fill = survey),
              alpha = 0.2)+
  geom_line(aes(x = year,
                y = index,
                colour = survey))+
  facet_grid(rows = vars(Region),
             scales = "free")+
  theme_bw()+
  labs(title = "2005-2021")+
  xlab("")+
  ylab("Expected count on BBS")+
  scale_colour_viridis_d(begin = 0.1,end = 0.7,
                         aesthetics = c("fill","colour"))+  
  theme(plot.title = element_text(hjust = 1),
        axis.text = element_text(size = 8))

trajs <- trajs_l + trajs_s + plot_layout(guides = "collect")+
  labs(caption = "Estimated population trajectories from the BBS-only and integrated models.The long-term traejctories (left column) 
  both Survey-wide and for the Applachian region only, show relatively little variation between the two models. The more 
  recent population trajectory for the Appalachian mountains where the species-specifc survey contributes data, 
  shows more variation through time than the trajectory from only BBS data (right side).",
       alt = "Estimated population trajectories from the BBS-only and integrated models.
       The long-term traejctories (left column) both Survey-wide and for the Applachian
       region only, show relatively little variation between the two models. The more 
       recent population trajectory for the Appalachian mountains where the species-specifc
       survey exists, shows more variation through time than the trajectory from only
       BBS data (right side).")


png(filename = paste0("figures/images/trajectory_comparison.png"),
    res = 300,
    width = 8, height = 4.9, units = "in")
print(trajs)
dev.off()

trajs <- trajs_l + trajs_s + plot_layout(guides = "collect")


png(filename = paste0("trajectory_comparison.png"),
    res = 300,
    width = 8, height = 4.9, units = "in")
print(trajs)
dev.off()



# Plotting trends ---------------------------------------------------------

breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
labls = paste0(labls, " %")
col_viridis <- FALSE
if (col_viridis)
{
  map_palette <- c("#fde725", "#dce319", "#b8de29", "#95d840", "#73d055", "#55c667",
                   "#238a8d", "#2d708e", "#39568c", "#453781", "#481567")
}else
{
  map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                   "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
}

names(map_palette) <- labls


trends_int_plot <- trends_int$trends %>% 
  mutate(survey = "Integrated")

trends_plots_09 <- trends_09$trends %>% 
  mutate(survey = "BBS-only") 

trends_int_plot_09 <- trends_int_09$trends %>% 
  mutate(survey = "Integrated")

trends_plot <- trends$trends %>% 
  mutate(survey = "BBS-only") %>% 
  bind_rows(trends_int_plot,trends_plots_09,
            trends_int_plot_09) %>% 
  filter(region_type == "stratum") %>% 
  mutate(time_scale = paste(start_year,end_year,sep = "-"),
         Tplot = cut(trend,breaks = c(-Inf, breaks, Inf),labels = labls),
         Tplot_uci = cut(trend_q_0.9,breaks = c(-Inf, breaks, Inf),labels = labls),
         Tplot_lci = cut(trend_q_0.1,breaks = c(-Inf, breaks, Inf),labels = labls),
         halfwidth_of_80_percent_credible_interval = (trend_q_0.9-trend_q_0.1)/2)

trend_maps <- strata_map %>% 
  inner_join(trends_plot,
             by = c("strata_name" = "region"))



bb <- st_bbox(trend_maps)
base_state_map <- sf::st_transform(base_state_map,crs = st_crs(trend_maps))

trend_comp_map_l <- ggplot()+
  geom_sf(data = filter(trend_maps,time_scale == "1966-2021"),
          aes(fill = Tplot))+
  geom_sf(data = base_state_map,
          fill = NA)+
  coord_sf(xlim = data_bounding_box[c("xmin","xmax")],
           ylim = data_bounding_box[c("ymin","ymax")])+
  scale_fill_manual(values = map_palette, aesthetics = c("fill"),
                    guide = ggplot2::guide_legend(reverse=TRUE),
                    name = paste0("Trend\n1966-2021"))+
facet_grid(rows = vars(survey))+
  labs(title = "1966-2021")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 1),
        axis.text = element_text(size = 8))




trend_comp_map_s <- ggplot()+
  geom_sf(data = filter(trend_maps,time_scale == "2009-2021"),
          aes(fill = Tplot))+
  geom_sf(data = base_state_map,
          fill = NA)+
  coord_sf(xlim = data_bounding_box[c("xmin","xmax")],
           ylim = data_bounding_box[c("ymin","ymax")])+
  scale_fill_manual(values = map_palette, aesthetics = c("fill"),
                    guide = ggplot2::guide_legend(reverse=TRUE),
                    name = paste0("Trend\n2009-2021"))+
  facet_grid(rows = vars(survey))+
  labs(title = "2009-2021")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 1),
        axis.text = element_text(size = 8))


trend_comp_map <- trend_comp_map_l + trend_comp_map_s +
  labs(caption = "Spatial pattern in long-term (1966-2021) and more recent (2009-2021) population trends, from the two
       models. The long-term trends from the BBS-only and the integrated models are essentially the same. The
       more recent trends show more variation in space with the integration of data from both the BBS and
       the species-specific survey.",
       alt = "Spatial pattern in long-term (1966-2021) and more recent (2009-2021) population trends, from the two
       models. The long-term trends from the BBS-only and the integrated models are essentially the same. The
       more recent trends show more variation in space with the integration of data from both the BBS and
       the species-specific survey.")



  

png(filename = paste0("figures/images/trend_comparison.png"),
    res = 300,
    width = 8, height = 4.9, units = "in")
print(trend_comp_map)
dev.off()


trend_comp_map <- trend_comp_map_l + trend_comp_map_s 





png(filename = paste0("trend_comparison.png"),
    res = 300,
    width = 8, height = 4.9, units = "in")
print(trend_comp_map)
dev.off()




trend_comp_map_uci <- ggplot()+
  geom_sf(data = trend_maps,
          aes(fill = Tplot_uci))+
  geom_sf(data = base_state_map,
          fill = NA)+
  coord_sf(xlim = data_bounding_box[c("xmin","xmax")],
           ylim = data_bounding_box[c("ymin","ymax")])+
  scale_fill_manual(values = map_palette, aesthetics = c("fill"),
                    guide = ggplot2::guide_legend(reverse=TRUE),
                    name = paste0("Trend"))+
  facet_grid(cols = vars(time_scale),
             rows = vars(survey))+
  theme_bw()
trend_comp_map_lci <- ggplot()+
  geom_sf(data = trend_maps,
          aes(fill = Tplot_lci))+
  geom_sf(data = base_state_map,
          fill = NA)+
  coord_sf(xlim = data_bounding_box[c("xmin","xmax")],
           ylim = data_bounding_box[c("ymin","ymax")])+
  scale_fill_manual(values = map_palette, aesthetics = c("fill"),
                    guide = ggplot2::guide_legend(reverse=TRUE),
                    name = paste0("Trend"))+
  facet_grid(cols = vars(time_scale),
             rows = vars(survey))+
  theme_bw()

trend_comp_map_uci + trend_comp_map_lci 




mean_abund_comp_map <- ggplot()+
  geom_sf(data = trend_maps,
          aes(fill = rel_abundance))+
  geom_sf(data = base_state_map,
          fill = NA)+
  coord_sf(xlim = data_bounding_box[c("xmin","xmax")],
           ylim = data_bounding_box[c("ymin","ymax")])+
  scale_fill_viridis_c(aesthetics = c("fill"),
                    guide = ggplot2::guide_legend(reverse=TRUE),
                    name = paste0("Relative Abundance"),
                    transform = "log10")+
  facet_grid(cols = vars(time_scale),
             rows = vars(survey))+
  theme_bw()

mean_abund_comp_map  


trend_ci_width_comp_map_l <- ggplot()+
  geom_sf(data = filter(trend_maps,time_scale == "1966-2021"),
          aes(fill = halfwidth_of_80_percent_credible_interval))+
  geom_sf(data = base_state_map,
          fill = NA)+
  coord_sf(xlim = data_bounding_box[c("xmin","xmax")],
           ylim = data_bounding_box[c("ymin","ymax")])+
  scale_fill_viridis_c(aesthetics = c("fill"),
                       guide = ggplot2::guide_colourbar(reverse=TRUE,
                                                        nbin = 6),
                       name = paste0("Half-width of \n80% CI trend\n1966-2021"))+
  facet_grid(rows = vars(survey))+
  labs(title = "1966-2021")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 1),
        axis.text = element_text(size = 8))
trend_ci_width_comp_map_s <- ggplot()+
  geom_sf(data = filter(trend_maps,time_scale == "2009-2021"),
          aes(fill = halfwidth_of_80_percent_credible_interval))+
  geom_sf(data = base_state_map,
          fill = NA)+
  coord_sf(xlim = data_bounding_box[c("xmin","xmax")],
           ylim = data_bounding_box[c("ymin","ymax")])+
  scale_fill_viridis_c(aesthetics = c("fill"),
                       guide = ggplot2::guide_colourbar(reverse=TRUE,
                                                        nbin = 6),
                       name = paste0("Half-width of \n80% CI trend\n2009-2021"))+
  facet_grid(rows = vars(survey))+
  labs(title = "2009-2021")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 1),
        axis.text = element_text(size = 8))

trend_ci_width_comp_map <- trend_ci_width_comp_map_l + trend_ci_width_comp_map_s +
  labs(caption = "Spatial pattern in the uncertainty of trend estimates for the long-term (1966-2021) and more 
       recent (2009-2021) population trends, from the two models. The uncertainty long-term trends from the 
       BBS-only and the integrated models are are very similar. The more recent trends are more
       precise when information from the two surveys are integrated.",
       alt = "Spatial pattern in the uncertainty of trend estimates for the long-term (1966-2021) and more 
       recent (2009-2021) population trends, from the two models. The uncertainty long-term trends from the 
       BBS-only and the integrated models are are very similar. The more recent trends are more
       precise when information from the two surveys are integrated.")



png(filename = paste0("figures/images/CI_halfwidth.png"),
    res = 300,
    width = 8, height = 4.5, units = "in")
print(trend_ci_width_comp_map)
dev.off()


trend_ci_width_comp_map <- trend_ci_width_comp_map_l + trend_ci_width_comp_map_s 



png(filename = paste0("CI_halfwidth.png"),
    res = 300,
    width = 8, height = 4.5, units = "in")
print(trend_ci_width_comp_map)
dev.off()


