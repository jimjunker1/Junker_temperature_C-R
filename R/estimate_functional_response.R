#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param FW_inputs
#' @return
#' @author jimjunker1
#' @export
# estimate_function_response <- function(FW_inputs) {
targets::tar_load(FW_inputs)

  # get the biomass specific consumption
  cBio_lists = FW_inputs[1] |> 
    future_map(~.x |> 
          map(~.x |> 
                group_by(doy, date_id) |> 
                dplyr::summarise(R_mean = mean(chla_mg_m, na.rm = TRUE),
                                 R_sigma = sd(chla_mg_m, na.rm = TRUE),
                                 C = unique(bio_mg_m),
                                 flux = unique(flux_mg_m_int)) |> 
                dplyr::mutate(cBio = flux/C)) |>
          bind_rows() |> 
          dplyr::summarise(C_mean = mean(C, na.rm = TRUE),
                           C_sigma = sd(C, na.rm = TRUE),
                           cBio_mean = mean(cBio, na.rm = TRUE),
                           cBio_sigma = sd(cBio, na.rm = TRUE)))
  
  
  cBioRate = function(pars, a0 = 0.01, h0 = 1, q0= 1, K0 = 10){
    derivs <- function(t,c,pars){
      with(as.list(c(pars,c)),{
        dc_dt = ((a*R^1+q)/(1+a*h*R^1+q))*C + r*R*(1-R/K)
        
        return(list(dc_dt = dc_dt))
  })}
  
  times = c(seq(0.5,365,0.5))
  i_0 = with(as.list(pars), a0, h0, q0, K0)
  y = (i_0)
  out <- ode(y = y, parms = pars, times = times, func = derivs)
as.data.frame(out)
  }
  # 
  # pars = list(R = (20*0.5*sin((1/120)*seq(0.5,365,0.5))^2)+10[1],
  #          C = (10*0.5*sin((1/120)*seq(0.5,365,0.5)+10)^2)+5[1], a = 0.01, h =1, r = 4, q = 1, K = 30)
  pars = list(R = 10, C = 5, a = 0.01, h =1, r = 4, q = 1, K = 30)
  debugonce(cBioRate)
  out <- cBioRate(pars)
  
  
}

  