#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param prod_data
#' @return
#' @author jimjunker1
#' @export
split_bio_dates <- function(prod_data) {

  int_spp_bio_boots = prod_data %>% purrr::map(~rlist::list.subset(.,grepl("Bint", names(.))) %>%
                                             flatten_df %>%
                                             dplyr::mutate(month_id = lubridate::month(as.Date(date_id))) |> named_group_split(boot_id)) |>
    rlist::list.subset(names(stream_order_list))
  return(int_spp_bio_boots)
}
