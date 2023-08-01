

library(bbsBayes2)
library(tidyverse)

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








