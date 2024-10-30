#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param FW_summaries
#' @param trips
#' @return
#' @author jimjunker1
#' @export
create_obs_ts <- function(FW_summaries, trips = 5) {
  
  light = FW_summaries$light
  temp_lists = FW_summaries$temp_lists
  obs_data = FW_summaries$FW_dataList

  create_obs_data <- function(df, trips = 1L, lightDf = NULL, tempDf = NULL, mCol = NULL, set.seed = 1312,...){
    doyVec = rep(light$time, times = trips) |> bind_cols(rep(1:trips, 365) |> sort()) |>  dplyr::mutate(doy = ...1+((...2-1)*365)) |> select(doy) |> unlist() 
    light_out = rep(lightDf$light_obs,trips) |> c() |> bind_cols(time = doyVec,light_obs = _)
    temp_out = rep(tempDf$tempC, trips) |> c() |> bind_cols(time = doyVec, tempC = _)
    
    if(1 %in% df$doy){
      df = df |>  rename(time = 'doy') |> 
        bind_rows(data.frame(time = c(365),
                             C_obs = mean(c(df$C_obs[1], df$C_obs[dim(df)[1]])),
                             sigma_C = sqrt(sum(c(df$sigma_C[1], df$sigma_C[dim(df)[1]]))/4),
                             R_obs = mean(c(df$R_obs[1], df$R_obs[dim(df)[1]])),
                             sigma_R = sqrt(sum(c(df$sigma_R[1], df$sigma_R[dim(df)[1]]))/4),
                             cB_obs = mean(c(df$cB_obs[1], df$cB_obs[dim(df)[1]])),
                             sigma_cB = sqrt(sum(c(df$sigma_cB[1], df$sigma_cB[dim(df)[1]]))/4)) |> 
                    dplyr::mutate(!!eval(mCol) := mean(c(unlist(df[mCol][1,]), unlist(df[mCol][dim(df)[1],]))))) |>
        arrange(time)
    } else{
    
    df = df |>  rename(time = 'doy') |> 
      bind_rows(data.frame(time = c(1,365),
                           C_obs = mean(c(df$C_obs[1], df$C_obs[dim(df)[1]])),
                           sigma_C = sqrt(sum(c(df$sigma_C[1], df$sigma_C[dim(df)[1]]))/4),
                           R_obs = mean(c(df$R_obs[1], df$R_obs[dim(df)[1]])),
                           sigma_R = sqrt(sum(c(df$sigma_R[1], df$sigma_R[dim(df)[1]]))/4),
                           cB_obs = mean(c(df$cB_obs[1], df$cB_obs[dim(df)[1]])),
                           sigma_cB = sqrt(sum(c(df$sigma_cB[1], df$sigma_cB[dim(df)[1]]))/4)) |> 
                  dplyr::mutate(!!eval(mCol) := mean(c(unlist(df[mCol][1,]), unlist(df[mCol][dim(df)[1],]))))) |>
      arrange(time)}
    
    
    timeVec = rep(df$time, times = trips) |> bind_cols(rep(1:trips, length(unlist(df$time))) |> sort()) |> setNames(nm = c('time','trip'))
    newTimeVec = c()
    newRVec = c()
    newCVec = c()
    newcBVec = c()
    
    for(i in 1:length(unique(df$time))){
      newTimeVec = c(newTimeVec, rep(unlist(df$time)[i], trips))
      newRVec = c(newRVec, rlnorm(trips, mean = log(unlist(df$R_obs)[i]^2/sqrt(unlist(df$sigma_R)[i]^2 + unlist(df$R_obs)[i]^2)), sd = sqrt(log(1+(unlist(df$sigma_R)[i]^2 / unlist(df$R_obs)[i]^2)))))
      newCVec = c(newCVec, rlnorm(trips, mean = log(unlist(df$C_obs)[i]^2/sqrt(unlist(df$sigma_C)[i]^2 + unlist(df$C_obs)[i]^2)), sd = sqrt(log(1+(unlist(df$sigma_C)[i]^2 / unlist(df$C_obs)[i]^2)))))
      newcBVec = c(newcBVec, rlnorm(trips, mean = log(unlist(df$cB_obs)[i]^2/sqrt(unlist(df$sigma_cB)[i]^2 + unlist(df$cB_obs)[i]^2)), sd = sqrt(log(1+(unlist(df$sigma_cB)[i]^2 / unlist(df$cB_obs)[i]^2)))))
    }
    newDf = do.call('cbind', list(rep(1:trips, length(unlist(df$time))), newTimeVec, newRVec, newCVec, newcBVec)) |> 
      data.frame() |>  
      setNames(nm = c('trip','time','R_obs','C_obs','cB_obs')) |> 
      left_join(rep(unlist(df[mCol]), trips) |> bind_cols(timeVec, M_obs = _)) |> 
      dplyr::mutate(time = time+((trip-1)*365)) |> arrange(time) |> dplyr::select(-trip)
    
    return(list(light = light_out, temp = temp_out,
                obs_data = newDf))
  }
  obs_lists = future_map2(obs_data, temp_lists, ~create_obs_data(.x, trips = trips, lightDf = light, tempDf = .y, mCol = 'M_nmean'))
  
  return(obs_lists)

}
