

library(bbsBayes2)
library(tidyverse)
library(patchwork)


species <- "Golden-winged Warbler"

model = "first_diff"

stratification <- "latlong"

model_variant <- "spatial"



s <- stratify(stratification,
              species)

p <- prepare_data(s,
                  min_n_routes = 1,
                  min_max_route_years = 2)

strata_map <- load_map(stratification)

ps <- prepare_spatial(p, strata_map)



pm <- prepare_model(ps, model, model_variant)

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





