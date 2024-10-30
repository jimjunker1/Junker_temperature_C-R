#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param obs_data
#' @return
#' @author jimjunker1
#' @export
summarise_CR <- function(FW_summaries) {

  CR_list = FW_summaries |> pluck('FW_dataList') |> 
    map(~.x |>
          group_by(doy) |> 
          dplyr::mutate(CRratio = log10(C_obs)/log10(R_obs)),
                        is = log10((flux_obs/C_obs)/R_obs))
  
  return(CR_list)

}
