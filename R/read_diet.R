#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param diet_file
#' @return
#' @author jimjunker1
#' @export
read_diet <- function(diet_file) {

  readRDS(diet_file) |> pluck("diet_seasonal_boot_split")

}
