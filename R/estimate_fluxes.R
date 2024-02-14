#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param prod_seasons
#' @param diet_matrices
#' @return
#' @author jimjunker1
#' @export
estimate_fluxes <- function(prod_seasons, diet_matrices) {

  # ++++ Helper functions ++++ ##
  prod_to_met <- function(production, NPE,...){
    production %>% 
      dplyr::mutate(across(contains("prod"), list(loss = ~.x/NPE), .names = "{.fn}")) %>% 
      dplyr::select(taxon_id, loss) %>% tibble::deframe() -> losses
  }
  
  boot_flux_function = function(mat, losses, resource_effs_vct,...){
    
    colnames(mat) %>% tibble %>%
      setNames(., nm = 'matrix_name') %>%
      dplyr::mutate(efficiencies = dplyr::recode(matrix_name, !!!resource_effs_vct, .default = 0.7)) -> x 
    spp_order = colnames(mat)
    setNames(x$efficiencies, nm = x$matrix_name)-> matrix_efficiencies
    
    resource_losses = setNames(rep(0, length(resource_effs_vct)), nm = names(resource_effs_vct))
    losses = c(losses,resource_losses)
    matrix_efficiencies = matrix_efficiencies[spp_order];losses = losses[spp_order]
    
    flux_obj = fluxing(mat = mat, losses = losses, efficiencies =  matrix_efficiencies,
                       bioms.losses = FALSE, bioms.prefs = FALSE, ef.level = 'prey', method = 'tbp')
    return(flux_obj)
    
  }
  ## ++++ End Helper functions ++++ ##
  # resource date frame with efficiences
  diet_item = list("amorphous_detritus",
                   "cyanobacteria",
                   "diatom",
                   "filamentous",
                   "green_algae",
                   "plant_material",
                   "animal")
  
  efficiencies = list(
    c(0.08, 0.1, 0.12),
    c(0.08, 0.1, 0.12),
    c(0.24, 0.3, 0.36),
    c(0.24, 0.3, 0.36),
    c(0.24, 0.3, 0.36),
    c(0.08, 0.1, 0.12),
    c(0.56, 0.7, 0.84)
  )
  resource_keyval = setNames(efficiencies, nm = unlist(diet_item))
  
  # estimate the beta distributions to draw diet AEs from
  set.seed(1312)
  beta_dist.pars = map(efficiencies, ~rriskDistributions::get.beta.par(p = c(0.025,0.5,0.975), q = unlist(.x), plot = FALSE, show.output = FALSE))
  res_effs = map(beta_dist.pars, ~rbeta(nboot, shape1 = .x[1], shape2 = .x[2])) %>% 
    setNames(., nm = unlist(diet_item)) %>% bind_cols %>% split(., seq(nrow(.))) %>% map(~unlist(.))
  
  #NPE quantile 0.025,0.5,0.975
  set.seed(1312)
  NPE_dist.pars = rriskDistributions::get.beta.par(p = c(0.025,0.5,0.975),q = c(0.4,0.45,0.5), show.output = FALSE, plot = FALSE)
  NPE = rbeta(nboot, shape1 = NPE_dist.pars[1], shape2 = NPE_dist.pars[2])

  #split into dates for each boot
  prod_seasons = prod_seasons |> future_map(~.x |> future_map(~.x |> named_group_split(boot_id)|> map(~.x |> named_group_split(date_id))))
  # convert production to metabolic fluxes
  fluxes = prod_seasons |> future_map(~.x %>% map(., ~.x %>% map2(., NPE, ~.x %>% map2(., .y, ~prod_to_met(.x, NPE = .y))))) %>% rlist::list.subset(names(stream_order_list))
# map all diets to to dates and calculate fluxes with full variability
flux_full = map2(fluxes, diet_matrices, function(x,y){
    map2(x,y, function(a,b){
    pmap(list(a,b,res_effs), ~..1 %>% map(\(x) boot_flux_function(mat = ..2, losses = x, resource_effs_vct = ..3)))})})

  return(flux_full = flux_full)

}
