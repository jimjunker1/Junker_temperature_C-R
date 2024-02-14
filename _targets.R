# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline
here::i_am("_targets.R")
# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tidyverse", "junkR", "qs", "rlist", "furrr") # packages that your targets need to run
  # format = "qs", # Optionally set the default storage format. qs is fast.
)

# tar_make_clustermq() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
options(clustermq.scheduler = "multiprocess")

# tar_make_future() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
future::plan(future.callr::callr)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_target(name = prod_file, here::here("data/prod_full_list.rds"), format = "file"),
  tar_target(name = prod_data, readRDS(prod_file), format = "qs"),
  tar_target(name = prod_seasons, split_prod_seasons(prod_data), format = "qs"),
  tar_target(name = diet_file, here::here("data/modeled-diets.rds"), format = "file"),
  tar_target(name = diet_df, read_diet(diet_file), format = "qs"),
  tar_target(name = diet_matrices, create_diet_matrices(diet_df, prod_seasons), format = "qs"),
  tar_target(name = fluxes, estimate_fluxes(prod_seasons, diet_matrices), format = "qs"),
  tar_target(name = compressed_fluxes, compress_fluxes(fluxes), format = 'qs'),
  tar_target(name = consumption, estimate_consumption(compressed_fluxes), format = 'qs'),
  tar_target(name = bio_data, split_bio_dates(prod_data), format = 'qs'),
  tar_target(name = chla_file, here::here("data/chla.rds"), format = "file"),
  tar_target(name = chla_data, readRDS(chla_file), format = 'qs')
  )
