#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param compressed_fluxes
#' @return
#' @author jimjunker1
#' @export
estimate_consumption <- function(compressed_fluxes) {

  green_diet_item = list("cyanobacteria",
                   "diatom",
                   "filamentous",
                   "green_algae",
                   "plant_material")
  
  green_consumption = compressed_fluxes |> map(~.x |> map(~.x |> dplyr::filter(prey %in% green_diet_item) |> group_by(boot_id, date_id) |> dplyr::summarise(flux_mg_m_int = sum(flux_mg_m_int, na.rm = TRUE)))) 
  
  return(green_consumption)

}
