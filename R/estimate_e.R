#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param diet_matrices
#' @param fluxes
#' @return
#' @author jimjunker1
#' @export
estimate_e <- function(diet_matrices, fluxes) {
  
  diet_item = list("cyanobacteria",
                   "diatom",
                   "filamentous",
                   "green_algae",
                   "plant_material")

  NPE = fluxes[['NPE']] |> data.frame(NPE = _) |> dplyr::mutate(boot_id = paste0("X",1:n())) 
  
  res_effs = fluxes[['res_effs']] |> 
    bind_rows(.id = 'boot_id') |> 
    dplyr::mutate(boot_id = paste0("X",boot_id)) |> 
    pivot_longer(-boot_id, names_to = "prey", values_to = "AE")  |> 
    left_join(NPE) |> 
    dplyr::mutate(GGE = AE*NPE) |> 
    dplyr::filter(prey %in% diet_item);rm(NPE)

  diet_vec = diet_matrices |> 
    future_map(~.x |> map(~.x |> map(~.x |> data.frame() |> dplyr::filter(rowSums(.) > 0) |> rownames_to_column('prey') |> pivot_longer(cols = -prey, names_to = 'consumer', values_to = 'perc')) |> bind_rows(.id = "boot_id")) |> bind_rows(.id = 'yr_third')) |> bind_rows(.id = 'site_id') |> dplyr::filter(prey %in% diet_item) |> dplyr::select(-yr_third)
  
  e_df = diet_vec |>
    left_join(res_effs |> dplyr::select(boot_id, prey, GGE), by = c("boot_id","prey"))|>
    group_by(site_id, boot_id, consumer) |>
    dplyr::summarise(e = weighted.mean(GGE,perc)) |> 
    dplyr::select(-consumer) |> 
    dplyr::filter(!is.na(e))
  
  return(e_df)
}
