#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param env_data
#' @param chla_data
#' @param biomasses
#' @param consumption
#' @param e_estimates
#' @return
#' @author jimjunker1
#' @export
clean_FW_inputs <- function(env_data, chla_data, biomasses, bodysizes, consumption) {

  sample_data = chla_data |>
    dplyr::rename(date_id = 'DATE', site_id = 'SITE', chla_mg_m = 'chla') |> 
    left_join(env_data[['sample_env']]) 
  
  sample_nas = sample_data|>
    dplyr::filter(is.na(tempC)) |> 
    dplyr::select(-light, -tempC) |> 
    left_join(env_data[['full_env']] |> 
                dplyr::select( doy, contains('tempC')) |> 
                rename_with(.data = _,~tolower(gsub('_tempC','',.x)), contains('tempC')) |> 
                pivot_longer(-doy, names_to = 'site_id', values_to = 'tempC')) |> 
    left_join(env_data[['full_env']] |> dplyr::select(doy, light))
  
  sample_full = bind_rows(sample_data |> dplyr::filter(!is.na(tempC)), sample_nas)
  
  FW_list = future_pmap(list(biomasses, consumption, bodysizes), ~pmap(list(..1,..2,..3), ~plyr::join_all(list(..1, ..2,..3, sample_full)) |>  dplyr::mutate(date_id = as.Date(date_id))))
  
  return(FW_list)
                 
  }
