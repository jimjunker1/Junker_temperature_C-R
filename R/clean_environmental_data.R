#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param env_raw
#' @return
#' @author jimjunker1
#' @export
clean_environmental_data <- function(env_raw_1, env_raw_2) {

  full_env = env_raw_1 |> group_by(doy) |> dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
  
  sample_env = env_raw_2 |> bind_rows(.id = 'site_id') |>
    dplyr::select(-c(JULIAN,SITE), date_id = 'DATE' ) |>
    dplyr::mutate(doy = lubridate::yday(date_id),
                  tempC = overkt_to_C(tempkt)) |> 
    left_join(full_env |> dplyr::select(doy, light, Intensity), by = 'doy') |> 
    dplyr::select(site_id, date_id, doy, light, tempC)
  
  return(list(sample_env = sample_env,
              full_env = full_env))

}
