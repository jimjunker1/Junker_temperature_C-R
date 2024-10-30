#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param prod_data
#' @return
#' @author jimjunker1
#' @export
compress_bodysizes <- function(prod_data) {

  int_spp_M_boots = prod_data %>% purrr::map(~rlist::list.subset(.,grepl("Bint|Nint", names(.))) %>%
                                                # flatten %>% 
                                                plyr::join_all() %>%
                                                dplyr::mutate(spp_Mmean = bio_mg_m/n_ind_m) %>%
                                                group_by(boot_id, date_id) |> 
                                                dplyr::mutate(relB = bio_mg_m/sum(bio_mg_m)) |> 
                                                dplyr::summarise(M_nmean = sum(bio_mg_m)/sum(n_ind_m),
                                                                 M_bmean = weighted.mean(spp_Mmean, relB)) |> 
                                               named_group_split(boot_id))
  
  return(int_spp_M_boots)
  
}
