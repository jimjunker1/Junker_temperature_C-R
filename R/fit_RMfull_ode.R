#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param FW_inputs
#' @return
#' @author jimjunker1
#' @export
fit_RMfull_ode <- function(FW_inputs) {
   NULL
  # RMtemp_stan <- stan_model(here::here("R/RM_lighttemp_CRsamp.stan"))
  # 
  # x = FW_inputs[[1]][[1]] |> 
  #   group_by(date_id) |> 
  #   dplyr::summarise(across(where(is.numeric), mean)) |> 
  #   na.omit() |> 
  #   dplyr::rename(R = 'chla_mg_m', C = 'bio_mg_m')
  # 
  # 
  # RM_model <- odemodel(
  #   name = "RM_full",
  #   model = list(
  #     R ~ r*R*(1 - R/K) - ((a*R)/(1 + a*h*R))*C,
  #     C ~ e*((a*R)/(1+a*h*R))*C - m*C
  #   ),
  #   observation = list(
  #     chla ~ dnbinom(mu = R, size = size1),
  #     cons ~ dnbinom(mu = C, size = size2)
  #   ),
  #   initial = list(
  #     R ~ R0,
  #     C ~ C0
  #   ),
  #   par = c('r','K','a','h','e','m', 'R0','C0', 'size1','size2')
  # )
  # 
  # CRstart = c(r = 0.02, K = 30, a = 0.01, h = 0.5, e = 0.15, m = 0.01, R0 = 30, C0 = 300, size1 = 1, size2 = 1)
  # 
  # CRfit = fitode(RM_model, data = x, start = CRstart, tcol = "doy")
  
  
  # stan_data <- list(
  #   T = nrow(data),
  #   ts = seq(1, nrow(data)), # Assuming time steps are 1, 2, 3, ...
  #   prey_initial = data$prey_population[1], # Initial prey population
  #   predator_initial = data$predator_population[1], # Initial predator population
  #   prey_population = data$prey_population, # Prey population data
  #   predator_population = data$predator_population, # Predator population data
  #   temperature = data$temperature, # Temperature data
  #   light = data$light, # Light availability data
  #   prey_consumption = data$prey_consumption # Prey consumption data
  # )

}
