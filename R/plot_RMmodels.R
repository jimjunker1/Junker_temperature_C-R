#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author jimjunker1
#' @export
plot_RMmodels <- function(R_sims) {

  
  png(filename = here::here("ignore/CRplot.png"), height = 5, width = 5, units = "in",
      res = 350, bg = "transparent")
  x |> dplyr::filter(time > 0) |> 
    ggplot()+
    geom_path(aes(x = R, y = C, color = log10(time)), linewidth = 1.2)+
    # geom_smooth(aes(x = R, y = C, color = after_stat(..y..)), linewidth = 1, se = FALSE)+
    viridis::scale_color_viridis(option = 'inferno', begin = 0, end = 0.8)+
    theme_minimal()+
    theme(legend.position = 'none',
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 16),
          panel.background = element_rect(fill ='transparent', color = NA))
  
dev.off()

}
