#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param CRratios
#' @return
#' @author jimjunker1
#' @export
plot_CR <- function(CRratios) {

  CRplotDf = CRratios |> bind_rows(.id = 'site_id') |>
    dplyr::mutate(site_id = factor(site_id, levels = names(stream_temp_labels))) |> 
    ggplot()+
    geom_line(aes(x = doy, y = CRratio, group = site_id, color = site_id))+
    scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos])+
    scale_color_manual(values = ocecolors[['temperature']][oce_temp_pos])+
    facet_wrap(~site_id)
  
  ISplot = CRratios |> bind_rows(.id = 'site_id') %>%
    dplyr::mutate(site_id = factor(site_id, levels = names(stream_temp_labels))) |> 
    ggplot()+
    geom_line(aes(x = doy, y = is, group = site_id, color = site_id), linewidth = 1.2) +
    geom_point(aes(x = doy, y = is, group = site_id, fill = site_id), shape = 21, color = 'black', size = 2)+
    # geom_smooth(aes)
    scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos])+
    scale_color_manual(values = ocecolors[['temperature']][oce_temp_pos])+
    facet_wrap(~site_id)

}
