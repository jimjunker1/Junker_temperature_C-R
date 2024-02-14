#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param fluxes
#' @return
#' @author jimjunker1
#' @export
compress_fluxes <- function(fluxes) {

flux_list = future_map(fluxes, ~map(.x, ~map(.x, ~map(.x, ~.x |> data.frame() |>  rownames_to_column("prey") |> pivot_longer(cols = -prey, names_to = "consumer", values_to = "flux_mg_m_int") |> dplyr::filter(flux_mg_m_int > 0)) |> bind_rows(.id = 'date_id')) |> bind_rows(.id = "boot_id")) |> bind_rows(.id = "yr_third") |> dplyr::select(-yr_third) |> named_group_split(boot_id))

return(flux_list)

}
