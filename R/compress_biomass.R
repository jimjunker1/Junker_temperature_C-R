#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param bio_data
#' @return
#' @author jimjunker1
#' @export
compress_biomass <- function(bio_data) {
  
  predator_lists = list(hver = c("Limnophora riparia"),
                        st6 = c("Limnophora riparia", "Antocha sp.", "Ephydridae.sp.", "Clinocera stagnalis"),
                        st9 = c("Limnophora riparia", "Antocha sp.", "Ephydridae.sp.", "Clinocera stagnalis"),
                        st7 = c("Limnophora riparia", "Macropelopia", "Dicranota", "Clinocera stagnalis"),
                        oh2 = c("Limnophora riparia", "Macropelopia", "Dicranota", "Clinocera stagnalis"),
                        st14 = c("Limnophora riparia", "Antocha sp.")) %>%
    rlist::list.subset(names(stream_order_list))

  biomass = map2(.x = bio_data, .y = predator_lists, ~.x |> map(\(x) dplyr::filter(x, taxon_id %ni% ..2) |> group_by(site_id, boot_id, date_id) |> dplyr::summarise(bio_mg_m = sum(bio_mg_m, na.rm = TRUE))))
    
return(biomass)
}
