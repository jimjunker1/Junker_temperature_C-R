#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param name
#' @param 
#' @return
#' @author jimjunker1
#' @export
split_prod_seasons <- function(prod_data) {

  int_spp_boots = prod_data %>% purrr::map(~rlist::list.subset(.,grepl("Pint", names(.))) %>%
                               flatten_df %>%
                               dplyr::mutate(month_id = lubridate::month(as.Date(date_id)),
                                             yr_third = case_when(month_id %in% 1:4 ~ "first",
                                                                  month_id %in% 5:8 ~ "second",
                                                                  month_id %in% 9:12 ~ "third"),
                                             jul_day = julian(as.Date(date_id), origin = as.Date("2010-01-01")),
                                             y_day = lubridate::yday(date_id)) |> 
                                 named_group_split(yr_third)) |>
    # map(~.x |> named_group_split(date_id)) |> furrr::future_map(~.x |> named_group_split(boot_id))) %>%
    rlist::list.subset(names(stream_order_list))
  return(int_spp_boots)
}
