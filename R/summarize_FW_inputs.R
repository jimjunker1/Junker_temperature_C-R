#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param FW_inputs
#' @return
#' @author jimjunker1
#' @export
summarize_FW_inputs <- function(FW_inputs) {

  light = readRDS(here::here("data/doy_tempkt_light.rds")) |> dplyr::summarise(light = mean(light), .by = doy) |> dplyr::select(time = 'doy', light_obs = 'light')
  temp_lists = readRDS(here::here("data/doy_tempkt_light.rds")) |> dplyr::select(time = 'doy', contains('tempC')) |> group_by(time) |> dplyr::summarise(across(contains('tempC'), \(x) mean(x, na.rm = TRUE))) |> dplyr::rename_with(\(x) tolower(gsub("_tempC","",x))) |> pivot_longer(-time, names_to = 'site_id', values_to = 'tempC') |> named_group_split(site_id) |> rlist::list.subset(names(stream_temp_labels))
  
  ### summarise cBioLists

  cBio_lists = FW_inputs |> 
    future_map(~.x |>
          map(~.x |> 
                group_by(doy, date_id) |> 
                dplyr::summarise(R_mean = mean(chla_mg_m, na.rm = TRUE),
                                 R_sigma = sd(chla_mg_m, na.rm = TRUE),
                                 C = unique(bio_mg_m),
                                 flux = unique(flux_mg_m_int),
                                 M_nmean = unique(M_nmean),
                                 M_bmean = unique(M_bmean)) |> 
                dplyr::mutate(cBio = flux/C,
                              doy = yday(date_id))) |> 
          bind_rows() |> 
          dplyr::summarise(C_obs = mean(C, na.rm = TRUE),
                           sigma_C = sd(C, na.rm = TRUE),
                           flux_obs = mean(flux, na.rm = TRUE),
                           sigma_flux = sd(flux, na.rm = TRUE),
                           cB_obs = mean(cBio, na.rm = TRUE),
                           sigma_cB = sd(cBio, na.rm = TRUE),
                           M_nmean = mean(M_nmean),
                           M_bmean = mean(M_bmean)))
  
  x_R = FW_inputs |>
    map(~.x |> 
          pluck(1) |> 
          ungroup() |>
          dplyr::group_by(doy, date_id) |>
          dplyr::summarise(R_obs = mean(chla_mg_m, na.rm = TRUE),
                           sigma_R = sd(chla_mg_m))) |> 
    bind_rows(.id = 'site_id') |>
    arrange(date_id) |> ungroup() |> 
    dplyr::mutate(doy = yday(as.Date(date_id)),
                  R_obs = if_else(is.na(R_obs),
                                  approx(doy,R_obs,doy)$y, R_obs),
                  sigma_R = if_else(is.na(sigma_R), 
                                    approx(doy, sigma_R,doy)$y, sigma_R)) |> 
    named_group_split(site_id) |> rlist::list.subset(names(stream_temp_labels))
  
  st6BioDf = FW_inputs[['st6']] |>
    map(~.x |> group_by(doy, date_id) |> 
                 dplyr::summarise(C = unique(bio_mg_m),
                                  flux = unique(flux_mg_m_int),
                                  M_nmean = unique(M_nmean),
                                  M_bmean = unique(M_bmean)) |> 
                 dplyr::mutate(cBio = flux/C,
                               doy = yday(date_id))) |> 
    bind_rows() |> 
    dplyr::summarise(C_obs = mean(C, na.rm = TRUE),
                     sigma_C = sd(C, na.rm = TRUE),
                     flux_obs = mean(flux, na.rm = TRUE),
                     sigma_flux = sd(flux, na.rm = TRUE),
                     cB_obs = mean(cBio, na.rm = TRUE),
                     sigma_cB = sd(cBio, na.rm = TRUE),
                     M_nmean = mean(M_nmean),
                     M_bmean = mean(M_bmean))
  cBio_lists[['st6']] <- st6BioDf
  
  st6_r = FW_inputs['st6'] |>
    map(~.x |> 
          pluck(1) |> 
          ungroup() |>
          dplyr::group_by(doy, date_id) |>
          dplyr::summarise(R_obs = mean(chla_mg_m, na.rm = TRUE),
                           sigma_R = sd(chla_mg_m))) |> 
    bind_rows(.id = 'site_id') |>
    arrange(date_id) |> ungroup() |> 
    dplyr::mutate(doy = yday(as.Date(date_id)),
                  R_obs = if_else(is.na(R_obs),
                                  approx(doy,R_obs,doy)$y, R_obs),
                  sigma_R = if_else(is.na(sigma_R), 
                                    approx(doy, sigma_R,doy)$y, sigma_R)) |> na.omit()
  
  x_R[['st6']] <- st6_r
  
  
  FW_dataList = map2(x_R, cBio_lists, \(x,y){ z = left_join(x, y) |> arrange(doy) |> 
    dplyr::mutate(diff_time= c(NA,diff(doy)),
                  cB_obs = cB_obs/diff_time)
  z$cB_obs[1] <- mean(c(z$cB_obs[2], z$cB_obs[dim(z)[1]]))
  z}) |> 
    rlist::list.subset(names(stream_temp_labels))
  
  return(list(light = light, temp_lists = temp_lists, FW_dataList = FW_dataList))
}
