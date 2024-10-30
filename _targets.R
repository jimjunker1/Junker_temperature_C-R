# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline
here::i_am('_targets.R')
# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c('tidyverse', 'junkR', 'qs',
               'rlist', 'furrr', 'rstan','tidybayes',
               'deSolve','rootSolve','stats','dfoptim', 'FME',
               'rriskDistributions'), # packages that your targets need to run
  format = 'qs', # Optionally set the default storage format. qs is fast.
  garbage_collection = TRUE,
  memory = 'transient')

# tar_make_clustermq() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
options(clustermq.scheduler = 'multiprocess')

# tar_make_future() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
future::plan(future.callr::callr)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source('other_functions.R') # Source other scripts as needed.

# Replace the target list below with your own:
list(
  # summarise and clean all the consumer flux and biomass data
  tar_target(name = prod_file, here::here('data/prod_full_list.rds'), format = 'file'),
  tar_target(name = prod_data, readRDS(prod_file), format = 'qs'),
  tar_target(name = prod_seasons, split_prod_seasons(prod_data), format = 'qs'),
  tar_target(name = diet_file, here::here('data/modeled-diets.rds'), format = 'file'),
  tar_target(name = diet_df, read_diet(diet_file), format = 'qs'),
  tar_target(name = diet_matrices, create_diet_matrices(diet_df, prod_seasons), format = 'qs'),
  tar_target(name = fluxes, estimate_fluxes(prod_seasons, diet_matrices), format = 'qs'),
  tar_target(name = compressed_fluxes, compress_fluxes(fluxes), format = 'qs'),
  tar_target(name = consumption, estimate_consumption(compressed_fluxes), format = 'qs'),
  tar_target(name = bio_data, split_bio_dates(prod_data), format = 'qs'),
  tar_target(name = biomasses, compress_biomass(bio_data), format = 'qs'),
  tar_target(name = bodysizes, compress_bodysizes(prod_data), format ='qs'),
  # summarise the resource data
  tar_target(name = chla_file, here::here('data/chla.rds'), format = 'file'),
  tar_target(name = chla_data, readRDS(chla_file), format = 'qs'),
  # summarise and clean environmental data
  tar_target(name = env_file_1, here::here('data/doy_tempkt_light.rds'), format = 'file'),
  tar_target(name = env_file_2, here::here('data/stream_temps_full.rds'), format = 'file'),
  tar_target(name = env_raw_1, readRDS(env_file_1), format = 'qs'),
  tar_target(name = env_raw_2, readRDS(env_file_2), format = 'qs'),
  tar_target(name = env_data, clean_environmental_data(env_raw_1, env_raw_2), format = 'qs'),
  # summarise parameters to feed as priors 
  tar_target(name = e_estimates, estimate_e(diet_matrices, fluxes), format = 'qs'),
  # clean the consumer, resource, and environmental data
  tar_target(name = FW_inputs, clean_FW_inputs(env_data, chla_data, biomasses, bodysizes, consumption), format = 'qs'),
  tar_target(name = FW_summaries, summarize_FW_inputs(FW_inputs), format = 'qs'),
  tar_target(name = obs_data, create_obs_ts(FW_summaries, trips = 50), format = 'qs'),
  # estimate functional response
  # tar_target(name = Func_resp, estimate_functional_response(FW_inputs), format = 'qs'),
  # run the food web models
  # tar_target(name = RM_sims, simulate_RMmodels(), format = 'qs')#,
  # tar_target(name = RM_plots, plot_RMmodels(R_sims), format = 'qs'),
  # tar_target(name = RM_model_simp, fit_RMsimp_ode(obs_data), format = 'qs'),
  # tar_target(name = RM_model_full, fit_RMfull_ode(FW_inputs), format = 'qs')
  # tar_target(name = RM_summary, summarise_RMmodels(obs_data, RM_models), format = 'qs'),
  tar_target(name = CRratios, summarise_CR(FW_summaries), format = 'qs'),
  tar_target(name = CRplots, plot_CR(CRratios), format = 'qs'),
  NULL
  )
