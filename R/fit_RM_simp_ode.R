#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param FW_inputs
#' @return
#' @author jimjunker1
#' @export
fit_RM_ode <- function(obs_data) {

  #Hver ODE model
  # run it in the background after tweaking below
  rstudioapi::jobRunScript(
    path = here::here("ignore/hver_ODEfit.R"),
    name = "hver ODE",
    workingDir = here::here(),
    importEnv = FALSE,
  )
  hver_list = obs_data[['hver']] 
  hver_list$obs_data = hver_list$obs_data |> dplyr::mutate(across(c(R_obs, C_obs, cB_obs), ~.x*10)) 
  
  T_func <- approxfun(hver_list[['temp']]$time, hver_list[['temp']]$tempC, method = "linear", yleft = min(hver_list[['temp']]$tempC), yright = max(hver_list[['temp']]$tempC))
  
  L_func <- approxfun(hver_list[['light']]$time, hver_list[['light']]$light_obs, method = "linear", yleft = min(hver_list[['light']]$light_obs), yright = max(hver_list[['light']]$light_obs))
  
  M_func = approxfun(hver_list[['obs_data']]$time, hver_list[['obs_data']]$M_obs, method = 'linear', yleft = min(hver_list[['obs_data']]$M_obs), yright = max(hver_list[['obs_data']]$M_obs))
  
  ode_system <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      T <- T_func(time)
      L <- L_func(time)
      M <- M_func(time)
      # 
      # r <- r_base * exp(-Er/C_to_overkt(T))
      K <- K_base * L^K_scale #* L#min(c(1,(1-L/800)))#L
      # 
      # a <- a_base + a_scale * T#exp(a_T_coef * (T))
      m <- m_base * exp(-Em/C_to_overkt(T))* (M^-alpha) #
      
      # cB <- a_T * R * C / (1 + a_T * h * R)
      # dRdt <- r_TL * R * (1 - R / K_TL) - cB
      # dCdt <- e * cB - m_T * C
      # dRdt <- r_TL * R * (1 - R / K_TL) - a_T * R * C / (1 + a_T * h * R)
      # dCdt <- e * a_T * R * C / (1 + a_T * h * R) - m_T * C
      dRdt <- r* R * (1 - R / K) - a * R * C #/ (1 + a * h * R)
      dCdt <- e * a * R * C - m*C#/ (1 + a * h * R) - m * C
      
      # return(list(c(dRdt, dCdt, cB)))
      return(list(c(dRdt, dCdt)))
    })
  }
  
  loss_function <- function(parameters) {
    # Initial state
    state <- c(R = hver_list[['obs_data']]$R_obs[1], C = hver_list[['obs_data']]$C_obs[1])#, cB = obs_data$cB_obs[1]) 
    # Solve ODE
    out <- ode(y = state, times = hver_list[['temp']]$time, func = ode_system, parms = parameters, method ='bdf')# 'daspk')
    model_out = as.data.frame(out)
    # Extract model predictions
    tt = unlist(hver_list[['obs_data']]$time)
    R_pred <- model_out[which(model_out$time %in% tt), "R"]
    C_pred <- model_out[which(model_out$time %in% tt), "C"]
    # cB_pred <- model_out[, "cB"]
    
    # conform the observations 
    obsDf = hver_list[['obs_data']][which(hver_list[['obs_data']]$time %in% model_out$time),]
    
    # Calculate loss
    loss <- sum((hver_list[['obs_data']]$R_obs - R_pred)^2, na.rm = TRUE) + 
      sum((hver_list[['obs_data']]$C_obs - C_pred)^2, na.rm = TRUE) #+ 
      # sum((obs_data$cB_obs - cB_pred)^2, na.rm = TRUE)
    
    return(loss)
  }

  # initial_guesses <- c(r = 15, Er = 0.32, a = 0.001,
  #                      K_base = 200, m_base = 0.001, h = 0.01,
  #                      Em = 0.65, e = 0.15)
  initial_guesses <- c(r = 15, a = 0.0051,
                           K_base = 500, K_scale = 0.5, alpha = 0.25,
                           m_base = 0.05, Em = 0.65, e = 0.16)
  hver_model <- ode(y=c(R = mean(hver_list[['obs_data']]$R_obs[2]), C = mean(hver_list[['obs_data']]$C_obs[2])), times = hver_list[['temp']]$time, func = ode_system, parms = initial_guesses, method = 'daspk')
  compare_results(hver_model, hver_list[['obs_data']])
  plot(hver_model)
  
plot_modObs = function(modODE = NULL, obsDf = NULL,...){
  mod = as.data.frame(modODE) |> dplyr::mutate(type = 'model') |> dplyr::filter(time > 730)
  obs = obsDf |> dplyr::mutate(type = 'observed') |> dplyr::filter(time > 730) |> dplyr::select(time, R = 'R_obs', C = 'C_obs', type)
  
  plotDf = bind_rows(mod, obs) |> group_by(time, type) |> pivot_longer(c(C,R), names_to = 'variable', values_to = 'value')
  
  plotDf |> ggplot()+
    geom_line(aes(x = time, y = log10(value), color = type))+
    facet_wrap(~variable)+
    theme_minimal()
  
  
}
plot_modObs(modODE = hver_model, obsDf = hver_list[['obs_data']])

     optim_result <- optim(par = initial_guesses, fn = loss_function, method = "L-BFGS-B", 
                           lower = rep(0.0001, length(initial_guesses)), 
                           upper = c(60,#r
                                     0.1,#a
                                     500,#K_base
                                     5,#K_scale
                                     2,#alpha
                                     1,#m_base
                                     2,#Em
                                     1),#e,
                           control = list(maxit = 3e5))
     optim_result$par
     
    hverCostFun = function(pars){
      out <- as.data.frame(ode(y=c(R = hver_list[['obs_data']]$R_obs[1], C = hver_list[['obs_data']]$C_obs[1]), times = hver_list[['temp']]$time, func = ode_system, parms = initial_guesses, method = 'daspk'))
      cost <- modCost(model = out, obs = hver_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs'), err = NULL)
      return(modCost(model = out, obs = hver_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs'), err = NULL, cost = cost))}
    
    hverCost = hverCostFun(pars = initial_guesses)
    hverCostFun(initial_guesses)$model
    
    hverCost2<- function(lpars){
      hverCostFun(c(exp(lpars)))# out <- ode(y=c(R = P[['R0']], C = P[['C0']]), times = c(0, hver_list[['temp']]$time), func = ode_system, parms = P, method = 'bdf')#daspk')
      # return(modCost(out, hver_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs')))
    }
    # debugonce(modFit)
  fit <- modFit(f = hverCost2, p = c(log(initial_guesses), R0 = hver_list[['obs_data']]$R_obs[1],
                                      C0 = hver_list[['obs_data']]$C_obs[1]))
  exp(coef(fit))
  var0 = fit$var_ms_unweighted
  cov0 = summary(fit)$cov.sc
  MCMC <- modMCMC(f = hverCost2, p = fit$par, niter = 5000, wvar0 = 0.1,  var0 = var0, updatecov = 50)
  MCMC$pars <- exp(MCMC$pars)
  plot(MCMC, Full = TRUE)
  # parameter fitting using levenberg marquart algorithm

  
  
  # Assuming `optim_result$par` contains your optimized parameters
  model_predictions <- ode(y = c(R = hver_list[['obs_data']]$R_obs[1],
                                 C = hver_list[['obs_data']]$C_obs[1]), #cB = round(obs_data$cB_obs[1],3)),
                                 times = hver_list[['obs_data']]$time, func = ode_system, parms = optim_result$par)
  
  # Convert model output to a data frame
  model_predictions_df <- as.data.frame(model_predictions)
  # model_predictions_df$cB <- model_predictions_df$V3 # Assuming the third state variable is cB in your model output

  
  # Add observed data to the data frame
  model_predictions_df$R_obs = hver_list[['obs_data']]$R_obs
  model_predictions_df$C_obs = hver_list[['obs_data']]$C_obs
  # model_predictions_df$cB_obs = obs_data$cB_obs
  
  # Melt the data frame for ggplot2
  plot_data <- melt(model_predictions_df, id.vars = "time", variable.name = "Variable", value.name = "Value")
  
  # Create a new column to differentiate between observed and predicted values
  plot_data$Type <- ifelse(grepl("_obs", plot_data$Variable), "Observed", "Predicted")
  
  # Adjust the Variable column to have cleaner names
  plot_data$Variable <- gsub("_obs", "", plot_data$Variable)
  plot_data$Variable <- gsub("V[1-3]", c("R", "C", "cB"), plot_data$Variable)

  ggplot(plot_data, aes(x = time, y = Value, color = Type)) + 
    geom_line(aes(linetype = Type), size = 1) + 
    geom_point(aes(shape = Type), size = 3) +
    facet_wrap(~Variable, scales = "free_y") + 
    theme_minimal() + 
    labs(title = "Model Predictions vs. Observations", x = "Time", y = "Value") + 
    scale_color_manual(values = c("Observed" = "red", "Predicted" = "blue")) + 
    scale_shape_manual(values = c("Observed" = 16, "Predicted" = 17))
  
#######
  #st6 ODE model
  # run it in the background after tweaking below
  rstudioapi::jobRunScript(
    path = here::here("ignore/st6_ODEfit.R"),
    name = "st6 ODE",
    workingDir = here::here(),
    importEnv = FALSE,
  )
  
  st6_list = obs_data[['st6']] 
  st6_list$obs_data = st6_list$obs_data |> dplyr::mutate(across(c(R_obs, C_obs, cB_obs), ~.x*10)) 
  
  T_func <- approxfun(st6_list[['temp']]$time, st6_list[['temp']]$tempC, method = "linear", yleft = min(st6_list[['temp']]$tempC), yright = max(st6_list[['temp']]$tempC))
  
  L_func <- approxfun(st6_list[['light']]$time, st6_list[['light']]$light_obs, method = "linear", yleft = min(st6_list[['light']]$light_obs), yright = max(st6_list[['light']]$light_obs))
  
  M_func = approxfun(st6_list[['obs_data']]$time, st6_list[['obs_data']]$M_obs, method = 'linear', yleft = min(st6_list[['obs_data']]$M_obs), yright = max(st6_list[['obs_data']]$M_obs))
  
  ode_system <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      T <- T_func(time)
      L <- L_func(time)
      M <- M_func(time)
      # 
      # r <- r_base * exp(-Er/C_to_overkt(T)) * min(c(1,(1-L/200)))
      K <- K_base + K_scale * L#min(c(1,(1-L/800)))#L
      # 
      # a <- a_base * L#exp(a_T_coef * (T))
      m <- m_base * exp(-Em/C_to_overkt(T))* (M^-alpha) #
      
      # cB <- a_T * R * C / (1 + a_T * h * R)
      # dRdt <- r_TL * R * (1 - R / K_TL) - cB
      # dCdt <- e * cB - m_T * C
      # dRdt <- r_TL * R * (1 - R / K_TL) - a_T * R * C / (1 + a_T * h * R)
      # dCdt <- e * a_T * R * C / (1 + a_T * h * R) - m_T * C
      dRdt <- r* R * (1 - R / K) - a * R * C #/ (1 + a * h * R)
      dCdt <- e * a * R * C - m*C#/ (1 + a * h * R) - m * C
      
      # return(list(c(dRdt, dCdt, cB)))
      return(list(c(dRdt, dCdt)))
    })
  }
  
  loss_function <- function(parameters) {
    # Initial state
    state <- c(R = st6_list[['obs_data']]$R_obs[1], C = st6_list[['obs_data']]$C_obs[1])#, cB = obs_data$cB_obs[1]) 
    # Solve ODE
    out <- ode(y = state, times = st6_list[['temp']]$time, func = ode_system, parms = parameters, method = 'daspk')
    model_out = as.data.frame(out)
    # Extract model predictions
    tt = unlist(st6_list[['obs_data']]$time)
    R_pred <- model_out[which(model_out$time %in% tt), "R"]
    C_pred <- model_out[which(model_out$time %in% tt), "C"]
    # cB_pred <- model_out[, "cB"]
    
    # conform the observations 
    obsDf = st6_list[['obs_data']][which(st6_list[['obs_data']]$time %in% model_out$time),]
    
    # Calculate loss
    loss <- sqrt((sum((obsDf$R_obs - R_pred)^2, na.rm = TRUE) + 
                    sum((obsDf$C_obs - C_pred)^2, na.rm = TRUE)/nrow(obsDf))) #+ 
    # sum((obs_data$cB_obs - cB_pred)^2, na.rm = TRUE)
    
    return(loss)
  }
  
  initial_guesses <- c(r = 15, a = 0.0041,
                       K_base = 100, K_scale = 1.7, alpha = 0.25,
                       m_base = 0.004, Em = 0.45, e = 0.12)
  
  st6_model <- ode(y=c(R = st6_list[['obs_data']]$R_obs[1], C = st6_list[['obs_data']]$C_obs[1]), times = st6_list[['temp']]$time, func = ode_system, parms = initial_guesses, method = 'daspk')
  compare_results(st6_model, st6_list[['obs_data']])
  
  plot( st6_model)
  # remove first two years and compare mean and CVs of data and prediction
  
  
  optim_result <- optim(par = initial_guesses, fn = loss_function, 
                        method = "L-BFGS-B", 
                        lower = rep(0.00001, length(initial_guesses)), 
                        upper = c(30,#r
                                  0.1,#a
                                  200,#K_base
                                  5,#K_scale
                                  2,#alpha
                                  1,#m_base
                                  2,#Em
                                  1),#e
                        control = list(maxit = 5e2))
  optim_result$par
  saveRDS(optim_result, here::here("data/st6ODE_optim.rds"))
  optim_result <- readRDS(here::here("data/st6ODE_optim.rds"))
  
  st6_model_op <- ode(y=c(R = st6_list[['obs_data']]$R_obs[1], C = st6_list[['obs_data']]$C_obs[1]), times = st6_list[['temp']]$time, func = ode_system, parms = optim_result$par, method = 'daspk')
  
  optim_result2 <- optim(par = optim_result$par, fn = loss_function, 
                         method = "L-BFGS-B", 
                         lower = rep(0.00001, length(initial_guesses)), 
                         upper = c(20,2,0.1,200,0.5,1,2,1),
                         control = list(maxit = 1e5))
  
  optim_result2$par
  saveRDS(optim_result2, here::here("data/st6ODE_optim2.rds"))
  optim_result2 <- readRDS(here::here("data/st6ODE_optim2.rds"))
  
  st6CostFun = function(pars){
    out <- as.data.frame( st6_model_op)
    cost <- modCost(model = out, obs = st6_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs'), err = NULL)
    return(modCost(model = out, obs = st6_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs'), err = NULL, cost = cost))}
  
  
  st6Cost = st6CostFun(pars = initial_guesses)
  st6CostFun(initial_guesses)$model
  
  modPars = optim_result$par*1.2
  
  st6Cost2<- function(lpars){
    st6CostFun(c(exp(lpars)))
  }
  # debugonce(modFit)
  fit <- modFit(f = st6Cost2, p = c(log(modPars),
                                    R0 = st6_list[['obs_data']]$R_obs[1],
                                    C0 = st6_list[['obs_data']]$C_obs[1]))
  
  var0 = fit$var_ms_unweighted
  cov0 = summary(fit)
  # parameter fitting using levenberg marquart algorithm
  # initial guess for parameters
  # fitting
  fitval=nls2::nls(par=initial_guesses,fn=loss_function)
  # Check optimized parameters
  
  
  # Assuming `optim_result$par` contains your optimized parameters
  model_predictions <- ode(y = c(R = st6_list[['obs_data']]$R_obs[1],
                                 C = st6_list[['obs_data']]$C_obs[1]), #cB = round(obs_data$cB_obs[1],3)),
                           times = st6_list[['obs_data']]$time, func = ode_system, parms = optim_result$par)
  
  # Convert model output to a data frame
  model_predictions_df <- as.data.frame(model_predictions)
  # model_predictions_df$cB <- model_predictions_df$V3 # Assuming the third state variable is cB in your model output
  
  
  # Add observed data to the data frame
  model_predictions_df$R_obs = st6_list[['obs_data']]$R_obs
  model_predictions_df$C_obs = st6_list[['obs_data']]$C_obs
  # model_predictions_df$cB_obs = obs_data$cB_obs
  
  # Melt the data frame for ggplot2
  plot_data <- melt(model_predictions_df, id.vars = "time", variable.name = "Variable", value.name = "Value")
  
  # Create a new column to differentiate between observed and predicted values
  plot_data$Type <- ifelse(grepl("_obs", plot_data$Variable), "Observed", "Predicted")
  
  # Adjust the Variable column to have cleaner names
  plot_data$Variable <- gsub("_obs", "", plot_data$Variable)
  plot_data$Variable <- gsub("V[1-3]", c("R", "C", "cB"), plot_data$Variable)
  
  ggplot(plot_data, aes(x = time, y = Value, color = Type)) + 
    geom_line(aes(linetype = Type), size = 1) + 
    geom_point(aes(shape = Type), size = 3) +
    facet_wrap(~Variable, scales = "free_y") + 
    theme_minimal() + 
    labs(title = "Model Predictions vs. Observations", x = "Time", y = "Value") + 
    scale_color_manual(values = c("Observed" = "red", "Predicted" = "blue")) + 
    scale_shape_manual(values = c("Observed" = 16, "Predicted" = 17))
  
  
  #######
  #st9 ODE model
  # run it in the background after tweaking below
  rstudioapi::jobRunScript(
    path = here::here("ignore/st9_ODEfit.R"),
    name = "st9 ODE",
    workingDir = here::here(),
    importEnv = FALSE,
  )
  
  st9_list = obs_data[['st9']] 
  st9_list$obs_data = st9_list$obs_data |> dplyr::mutate(across(c(R_obs, C_obs, cB_obs), ~.x*10)) 
  
  T_func <- approxfun(st9_list[['temp']]$time, st9_list[['temp']]$tempC, method = "linear", yleft = min(st9_list[['temp']]$tempC), yright = max(st9_list[['temp']]$tempC))
  
  L_func <- approxfun(st9_list[['light']]$time, st9_list[['light']]$light_obs, method = "linear", yleft = min(st9_list[['light']]$light_obs), yright = max(st9_list[['light']]$light_obs))
  
  M_func = approxfun(st9_list[['obs_data']]$time, st9_list[['obs_data']]$M_obs, method = 'linear', yleft = min(st9_list[['obs_data']]$M_obs), yright = max(st9_list[['obs_data']]$M_obs))
  
  ode_system <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      T <- T_func(time)
      L <- L_func(time)
      M <- M_func(time)
      # 
      # r <- r_base * exp(-Er/C_to_overkt(T)) * min(c(1,(1-L/200)))
      K <- K_base + K_scale * L#min(c(1,(1-L/800)))#L
      # 
      # a <- a_base * L#exp(a_T_coef * (T))
      m <- m_base * exp(-Em/C_to_overkt(T))* (M^-alpha) #
      
      # cB <- a_T * R * C / (1 + a_T * h * R)
      # dRdt <- r_TL * R * (1 - R / K_TL) - cB
      # dCdt <- e * cB - m_T * C
      # dRdt <- r_TL * R * (1 - R / K_TL) - a_T * R * C / (1 + a_T * h * R)
      # dCdt <- e * a_T * R * C / (1 + a_T * h * R) - m_T * C
      dRdt <- r* R * (1 - R / K) - a * R * C #/ (1 + a * h * R)
      dCdt <- e * a * R * C - m*C#/ (1 + a * h * R) - m * C
      
      # return(list(c(dRdt, dCdt, cB)))
      return(list(c(dRdt, dCdt)))
    })
  }
  
  loss_function <- function(parameters) {
    # Initial state
    state <- c(R = st9_list[['obs_data']]$R_obs[1], C = st9_list[['obs_data']]$C_obs[1])#, cB = obs_data$cB_obs[1]) 
    # Solve ODE
    out <- ode(y = state, times = st9_list[['temp']]$time, func = ode_system, parms = parameters, method = 'daspk')
    model_out = as.data.frame(out)
    # Extract model predictions
    tt = unlist(st9_list[['obs_data']]$time)
    R_pred <- model_out[which(model_out$time %in% tt), "R"]
    C_pred <- model_out[which(model_out$time %in% tt), "C"]
    # cB_pred <- model_out[, "cB"]
    
    # conform the observations 
    obsDf = st9_list[['obs_data']][which(st9_list[['obs_data']]$time %in% model_out$time),]
    
    # Calculate loss
    loss <- sqrt((sum((obsDf$R_obs - R_pred)^2, na.rm = TRUE) + 
      sum((obsDf$C_obs - C_pred)^2, na.rm = TRUE)/nrow(obsDf))) #+ 
    # sum((obs_data$cB_obs - cB_pred)^2, na.rm = TRUE)
    
    return(loss)
  }
  
  initial_guesses <- c(r = 10, a = 0.002,
                       K_base = 50, K_scale = 0.7, alpha = 0.25,
                       m_base = 0.012, Em = 0.65, e = 0.15)
  
  st9_model <- ode(y=c(R = st9_list[['obs_data']]$R_obs[1], C = st9_list[['obs_data']]$C_obs[1]), times = st9_list[['temp']]$time, func = ode_system, parms = st9ODE_optim$par, method = 'daspk')
  
  plot( st9_model)
  # remove first two years and compare mean and CVs of data and prediction
  compare_results(st9_model, st9_list[['obs_data']])
  
  
  optim_result <- optim(par = initial_guesses, fn = loss_function, 
                        method = "L-BFGS-B", 
                        lower = rep(0.00001, length(initial_guesses)), 
                        upper = c(20,0.1,200,1,2,0.5,2,1),
                        control = list(maxit = 1e5))
  optim_result$par
  saveRDS(optim_result, here::here("data/st9ODE_optim.rds"))
  optim_result <- readRDS(here::here("data/st9ODE_optim.rds"))
  
  st9_model_op <- ode(y=c(R = st9_list[['obs_data']]$R_obs[1], C = st9_list[['obs_data']]$C_obs[1]), times = st9_list[['temp']]$time, func = ode_system, parms = optim_result$par, method = 'daspk')
  
  optim_result2 <- optim(par = optim_result$par, fn = loss_function, 
                        method = "L-BFGS-B", 
                        lower = rep(0.00001, length(initial_guesses)), 
                        upper = c(20,2,0.1,200,0.5,1,2,1),
                        control = list(maxit = 1e5))
  
  optim_result2$par
  saveRDS(optim_result2, here::here("data/st9ODE_optim2.rds"))
  optim_result2 <- readRDS(here::here("data/st9ODE_optim2.rds"))
  
  st9CostFun = function(pars){
    out <- as.data.frame( st9_model_op)
    cost <- modCost(model = out, obs = st9_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs'), err = NULL)
    return(modCost(model = out, obs = st9_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs'), err = NULL, cost = cost))}
  
  
  st9Cost = st9CostFun(pars = initial_guesses)
  st9CostFun(initial_guesses)$model
  
  modPars = optim_result$par*1.2
  
  st9Cost2<- function(lpars){
    st9CostFun(c(exp(lpars)))
  }
  # debugonce(modFit)
  fit <- modFit(f = st9Cost2, p = c(log(modPars),
                                    R0 = st9_list[['obs_data']]$R_obs[1],
                                    C0 = st9_list[['obs_data']]$C_obs[1]))
  
  var0 = fit$var_ms_unweighted
  cov0 = summary(fit)
  # parameter fitting using levenberg marquart algorithm
  # initial guess for parameters
  # fitting
  fitval=nls2::nls(par=initial_guesses,fn=loss_function)
  # Check optimized parameters
  
  
  # Assuming `optim_result$par` contains your optimized parameters
  model_predictions <- ode(y = c(R = st9_list[['obs_data']]$R_obs[1],
                                 C = st9_list[['obs_data']]$C_obs[1]), #cB = round(obs_data$cB_obs[1],3)),
                           times = st9_list[['obs_data']]$time, func = ode_system, parms = optim_result$par)
  
  # Convert model output to a data frame
  model_predictions_df <- as.data.frame(model_predictions)
  # model_predictions_df$cB <- model_predictions_df$V3 # Assuming the third state variable is cB in your model output
  
  
  # Add observed data to the data frame
  model_predictions_df$R_obs = st9_list[['obs_data']]$R_obs
  model_predictions_df$C_obs = st9_list[['obs_data']]$C_obs
  # model_predictions_df$cB_obs = obs_data$cB_obs
  
  # Melt the data frame for ggplot2
  plot_data <- melt(model_predictions_df, id.vars = "time", variable.name = "Variable", value.name = "Value")
  
  # Create a new column to differentiate between observed and predicted values
  plot_data$Type <- ifelse(grepl("_obs", plot_data$Variable), "Observed", "Predicted")
  
  # Adjust the Variable column to have cleaner names
  plot_data$Variable <- gsub("_obs", "", plot_data$Variable)
  plot_data$Variable <- gsub("V[1-3]", c("R", "C", "cB"), plot_data$Variable)
  
  ggplot(plot_data, aes(x = time, y = Value, color = Type)) + 
    geom_line(aes(linetype = Type), size = 1) + 
    geom_point(aes(shape = Type), size = 3) +
    facet_wrap(~Variable, scales = "free_y") + 
    theme_minimal() + 
    labs(title = "Model Predictions vs. Observations", x = "Time", y = "Value") + 
    scale_color_manual(values = c("Observed" = "red", "Predicted" = "blue")) + 
    scale_shape_manual(values = c("Observed" = 16, "Predicted" = 17))
  
  
  #######
  #st7 ODE model
  # run it in the background after tweaking below
  rstudioapi::jobRunScript(
    path = here::here("ignore/st7_ODEfit.R"),
    name = "st7 ODE",
    workingDir = here::here(),
    importEnv = FALSE,
  )
  
  st7_list = obs_data[['st7']] 
  st7_list$obs_data = st7_list$obs_data |> dplyr::mutate(across(c(R_obs, C_obs, cB_obs), ~.x*10)) 
  
  T_func <- approxfun(st7_list[['temp']]$time, st7_list[['temp']]$tempC, method = "linear", yleft = min(st7_list[['temp']]$tempC), yright = max(st7_list[['temp']]$tempC))
  
  L_func <- approxfun(st7_list[['light']]$time, st7_list[['light']]$light_obs, method = "linear", yleft = min(st7_list[['light']]$light_obs), yright = max(st7_list[['light']]$light_obs))
  
  M_func = approxfun(st7_list[['obs_data']]$time, st7_list[['obs_data']]$M_obs, method = 'linear', yleft = min(st7_list[['obs_data']]$M_obs), yright = max(st7_list[['obs_data']]$M_obs))
  
  ode_system <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      T <- T_func(time)
      L <- L_func(time)
      M <- M_func(time)
      # 
      # r <- r_base * exp(-Er/C_to_overkt(T)) * min(c(1,(1-L/200)))
      K <- K_base + K_scale * L#min(c(1,(1-L/800)))#L
      # 
      # a <- a_base * L#exp(a_T_coef * (T))
      m <- m_base * exp(-Em/C_to_overkt(T))* (M^-alpha) #
      
      # cB <- a_T * R * C / (1 + a_T * h * R)
      # dRdt <- r_TL * R * (1 - R / K_TL) - cB
      # dCdt <- e * cB - m_T * C
      # dRdt <- r_TL * R * (1 - R / K_TL) - a_T * R * C / (1 + a_T * h * R)
      # dCdt <- e * a_T * R * C / (1 + a_T * h * R) - m_T * C
      dRdt <- r* R * (1 - R / K) - a * R * C #/ (1 + a * h * R)
      dCdt <- e * a * R * C - m*C#/ (1 + a * h * R) - m * C
      
      # return(list(c(dRdt, dCdt, cB)))
      return(list(c(dRdt, dCdt)))
    })
  }
  
  loss_function <- function(parameters) {
    # Initial state
    state <- c(R = st7_list[['obs_data']]$R_obs[1], C = st7_list[['obs_data']]$C_obs[1])#, cB = obs_data$cB_obs[1]) 
    # Solve ODE
    out <- ode(y = state, times = st7_list[['temp']]$time, func = ode_system, parms = parameters, method = 'daspk')
    model_out = as.data.frame(out)
    # Extract model predictions
    tt = unlist(st7_list[['obs_data']]$time)
    R_pred <- model_out[which(model_out$time %in% tt), "R"]
    C_pred <- model_out[which(model_out$time %in% tt), "C"]
    # cB_pred <- model_out[, "cB"]
    
    # conform the observations 
    obsDf = st7_list[['obs_data']][which(st7_list[['obs_data']]$time %in% model_out$time),]
    
    # Calculate loss
    loss <- sqrt((sum((obsDf$R_obs - R_pred)^2, na.rm = TRUE) + 
                    sum((obsDf$C_obs - C_pred)^2, na.rm = TRUE)/nrow(obsDf))) #+ 
    # sum((obs_data$cB_obs - cB_pred)^2, na.rm = TRUE)
    
    return(loss)
  }
  
  initial_guesses <- c(r = 40, a = 0.0025,
                       K_base = 200, K_scale = 1, alpha = 0.25,
                       m_base = 0.01, Em = 0.43, e = 0.20)
  
  st7_model <- ode(y=c(R = st7_list[['obs_data']]$R_obs[1], C = st7_list[['obs_data']]$C_obs[1]), times = st7_list[['temp']]$time, func = ode_system, parms = initial_guesses, method = 'daspk')
  #compare summary stats to observed
  compare_results(st7_model, st7_list[['obs_data']])
  
  plot( st7_model)
  # remove first two years and compare mean and CVs of data and prediction
  
  optim_result <- optim(par = initial_guesses, fn = loss_function, 
                        method = "L-BFGS-B", 
                        lower = rep(0.00001, length(initial_guesses)), 
                        upper = c(60,#r
                                  0.1,#a
                                  500,#K_base
                                  5,#K_scale
                                  2,#alpha
                                  1,#m_base
                                  2,#Em
                                  1),#e,
                        control = list(maxit = 5e2))
  optim_result$par
  saveRDS(optim_result, here::here("data/st7ODE_optim.rds"))
  optim_result <- readRDS(here::here("data/st7ODE_optim.rds"))
  
  st7_model_op <- ode(y=c(R = st7_list[['obs_data']]$R_obs[1], C = st7_list[['obs_data']]$C_obs[1]), times = st7_list[['temp']]$time, func = ode_system, parms = optim_result$par, method = 'daspk')
  
  optim_result2 <- optim(par = optim_result$par, fn = loss_function, 
                         method = "L-BFGS-B", 
                         lower = rep(0.00001, length(initial_guesses)), 
                         upper = c(60,#r
                                   0.1,#a
                                   500,#K_base
                                   5,#K_scale
                                   2,#alpha
                                   1,#m_base
                                   2,#Em
                                   1),#e,
                         control = list(maxit = 5e2))
  
  optim_result2$par
  saveRDS(optim_result2, here::here("data/st7ODE_optim2.rds"))
  optim_result2 <- readRDS(here::here("data/st7ODE_optim2.rds"))
  
  st7CostFun = function(pars){
    out <- as.data.frame(ode(y=c(R = st7_list[['obs_data']]$R_obs[1], C = st7_list[['obs_data']]$C_obs[1]), times = st7_list[['temp']]$time, func = ode_system, parms = pars, method = 'daspk'))
    cost <- modCost(model = out, obs = st7_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs'), err = NULL)
    return(modCost(model = out, obs = st7_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs'), err = NULL, cost = cost))}
  
  
  st7Cost = st7CostFun(pars = optim_result$par)
  st7CostFun(optim_result$par)$model
  
  modPars = initial_guesses*1.2
  
  st7Cost2<- function(lpars){
    st7CostFun(c(exp(lpars)))
  }
  # debugonce(modFit)
  fit <- modFit(f = st7Cost2, p = c(log(modPars),
                                    R0 = st7_list[['obs_data']]$R_obs[1],
                                    C0 = st7_list[['obs_data']]$C_obs[1]))
  
  exp(coef(fit))
  var0 = fit$var_ms_unweighted
  cov0 = summary(fit)$cov.scaled*2.4^2/5
  MCMC <- modMCMC(f = hverCost2, p = fit$par, niter = 5000, wvar0 = 0.1,  var0 = var0, updatecov = 50)
  MCMC$pars <- exp(MCMC$pars)
  saveRDS(list(fit = fit,MCMC = MCMC), here::here("data/st7_fmeFit.rds"))
  plot(MCMC, Full = TRUE)

  
  # Assuming `optim_result$par` contains your optimized parameters
  model_predictions <- ode(y = c(R = st7_list[['obs_data']]$R_obs[1],
                                 C = st7_list[['obs_data']]$C_obs[1]), #cB = round(obs_data$cB_obs[1],3)),
                           times = st7_list[['obs_data']]$time, func = ode_system, parms = optim_result$par)
  
  # Convert model output to a data frame
  model_predictions_df <- as.data.frame(model_predictions)
  # model_predictions_df$cB <- model_predictions_df$V3 # Assuming the third state variable is cB in your model output
  
  
  # Add observed data to the data frame
  model_predictions_df$R_obs = st7_list[['obs_data']]$R_obs
  model_predictions_df$C_obs = st7_list[['obs_data']]$C_obs
  # model_predictions_df$cB_obs = obs_data$cB_obs
  
  # Melt the data frame for ggplot2
  plot_data <- melt(model_predictions_df, id.vars = "time", variable.name = "Variable", value.name = "Value")
  
  # Create a new column to differentiate between observed and predicted values
  plot_data$Type <- ifelse(grepl("_obs", plot_data$Variable), "Observed", "Predicted")
  
  # Adjust the Variable column to have cleaner names
  plot_data$Variable <- gsub("_obs", "", plot_data$Variable)
  plot_data$Variable <- gsub("V[1-3]", c("R", "C", "cB"), plot_data$Variable)
  
  ggplot(plot_data, aes(x = time, y = Value, color = Type)) + 
    geom_line(aes(linetype = Type), size = 1) + 
    geom_point(aes(shape = Type), size = 3) +
    facet_wrap(~Variable, scales = "free_y") + 
    theme_minimal() + 
    labs(title = "Model Predictions vs. Observations", x = "Time", y = "Value") + 
    scale_color_manual(values = c("Observed" = "red", "Predicted" = "blue")) + 
    scale_shape_manual(values = c("Observed" = 16, "Predicted" = 17))
  #######
  #oh2 ODE model
  # run it in the background after tweaking below
  rstudioapi::jobRunScript(
    path = here::here("ignore/oh2_ODEfit.R"),
    name = "oh2 ODE",
    workingDir = here::here(),
    importEnv = FALSE,
  )
  
  oh2_list = obs_data[['oh2']] 
  oh2_list$obs_data = oh2_list$obs_data |> dplyr::mutate(across(c(R_obs, C_obs, cB_obs), ~.x*10)) 
  
  T_func <- approxfun(oh2_list[['temp']]$time, oh2_list[['temp']]$tempC, method = "linear", yleft = min(oh2_list[['temp']]$tempC), yright = max(oh2_list[['temp']]$tempC))
  
  L_func <- approxfun(oh2_list[['light']]$time, oh2_list[['light']]$light_obs, method = "linear", yleft = min(oh2_list[['light']]$light_obs), yright = max(oh2_list[['light']]$light_obs))
  
  M_func = approxfun(oh2_list[['obs_data']]$time, oh2_list[['obs_data']]$M_obs, method = 'linear', yleft = min(oh2_list[['obs_data']]$M_obs), yright = max(oh2_list[['obs_data']]$M_obs))
  
  ode_system <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      T <- T_func(time)
      L <- L_func(time)
      M <- M_func(time)
      # 
      # r <- r_base * exp(-Er/C_to_overkt(T)) * min(c(1,(1-L/200)))
      K <- K_base + K_scale * L#min(c(1,(1-L/800)))#L
      # 
      # a <- a_base * L#exp(a_T_coef * (T))
      m <- m_base * exp(-Em/C_to_overkt(T))* (M^-alpha) #
      
      # cB <- a_T * R * C / (1 + a_T * h * R)
      # dRdt <- r_TL * R * (1 - R / K_TL) - cB
      # dCdt <- e * cB - m_T * C
      # dRdt <- r_TL * R * (1 - R / K_TL) - a_T * R * C / (1 + a_T * h * R)
      # dCdt <- e * a_T * R * C / (1 + a_T * h * R) - m_T * C
      dRdt <- r* R * (1 - R / K) - a * R * C #/ (1 + a * h * R)
      dCdt <- e * a * R * C - m*C#/ (1 + a * h * R) - m * C
      
      # return(list(c(dRdt, dCdt, cB)))
      return(list(c(dRdt, dCdt)))
    })
  }
  
  loss_function <- function(parameters) {
    # Initial state
    state <- c(R = oh2_list[['obs_data']]$R_obs[1], C = oh2_list[['obs_data']]$C_obs[1])#, cB = obs_data$cB_obs[1]) 
    # Solve ODE
    out <- ode(y = state, times = oh2_list[['temp']]$time, func = ode_system, parms = parameters, method = 'daspk')
    model_out = as.data.frame(out)
    # Extract model predictions
    tt = unlist(oh2_list[['obs_data']]$time)
    R_pred <- model_out[which(model_out$time %in% tt), "R"]
    C_pred <- model_out[which(model_out$time %in% tt), "C"]
    # cB_pred <- model_out[, "cB"]
    
    # conform the observations 
    obsDf = oh2_list[['obs_data']][which(oh2_list[['obs_data']]$time %in% model_out$time),]
    
    # Calculate loss
    loss <- sqrt((sum((obsDf$R_obs - R_pred)^2, na.rm = TRUE) + 
                    sum((obsDf$C_obs - C_pred)^2, na.rm = TRUE)/nrow(obsDf))) #+ 
    # sum((obs_data$cB_obs - cB_pred)^2, na.rm = TRUE)
    
    return(loss)
  }
  
  initial_guesses <- c(r = 40, a = 0.0025,
                       K_base = 200, K_scale = 1, alpha = 0.25,
                       m_base = 0.01, Em = 0.43, e = 0.20)
  
  oh2_model <- ode(y=c(R = oh2_list[['obs_data']]$R_obs[1], C = oh2_list[['obs_data']]$C_obs[1]), times = oh2_list[['temp']]$time, func = ode_system, parms = initial_guesses, method = 'daspk')
  #compare summary stats to observed
  compare_results(oh2_model, oh2_list[['obs_data']])
  
  plot( oh2_model)
  # remove first two years and compare mean and CVs of data and prediction
  
  optim_result <- optim(par = initial_guesses, fn = loss_function, 
                        method = "L-BFGS-B", 
                        lower = rep(0.00001, length(initial_guesses)), 
                        upper = c(60,#r
                                  0.1,#a
                                  500,#K_base
                                  5,#K_scale
                                  2,#alpha
                                  1,#m_base
                                  2,#Em
                                  1),#e,
                        control = list(maxit = 5e2))
  optim_result$par
  saveRDS(optim_result, here::here("data/oh2ODE_optim.rds"))
  optim_result <- readRDS(here::here("data/oh2ODE_optim.rds"))
  
  oh2_model_op <- ode(y=c(R = oh2_list[['obs_data']]$R_obs[1], C = oh2_list[['obs_data']]$C_obs[1]), times = oh2_list[['temp']]$time, func = ode_system, parms = optim_result$par, method = 'daspk')
  
  optim_result2 <- optim(par = optim_result$par, fn = loss_function, 
                         method = "L-BFGS-B", 
                         lower = rep(0.00001, length(initial_guesses)), 
                         upper = c(60,#r
                                   0.1,#a
                                   500,#K_base
                                   5,#K_scale
                                   2,#alpha
                                   1,#m_base
                                   2,#Em
                                   1),#e,
                         control = list(maxit = 5e2))
  
  optim_result2$par
  saveRDS(optim_result2, here::here("data/oh2ODE_optim2.rds"))
  optim_result2 <- readRDS(here::here("data/oh2ODE_optim2.rds"))
  
  oh2CostFun = function(pars){
    out <- as.data.frame( oh2_model_op)
    cost <- modCost(model = out, obs = oh2_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs'), err = NULL)
    return(modCost(model = out, obs = oh2_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs'), err = NULL, cost = cost))}
  
  
  oh2Cost = oh2CostFun(pars = initial_guesses)
  oh2CostFun(initial_guesses)$model
  
  modPars = optim_result$par*1.2
  
  oh2Cost2<- function(lpars){
    oh2CostFun(c(exp(lpars)))
  }
  # debugonce(modFit)
  fit <- modFit(f = oh2Cost2, p = c(log(modPars),
                                    R0 = oh2_list[['obs_data']]$R_obs[1],
                                    C0 = oh2_list[['obs_data']]$C_obs[1]))
  
  var0 = fit$var_ms_unweighted
  cov0 = summary(fit)
  # parameter fitting using levenberg marquart algorithm
  # initial guess for parameters
  # fitting
  fitval=nls2::nls(par=initial_guesses,fn=loss_function)
  # Check optimized parameters
  
  
  # Assuming `optim_result$par` contains your optimized parameters
  model_predictions <- ode(y = c(R = oh2_list[['obs_data']]$R_obs[1],
                                 C = oh2_list[['obs_data']]$C_obs[1]), #cB = round(obs_data$cB_obs[1],3)),
                           times = oh2_list[['obs_data']]$time, func = ode_system, parms = optim_result$par)
  
  # Convert model output to a data frame
  model_predictions_df <- as.data.frame(model_predictions)
  # model_predictions_df$cB <- model_predictions_df$V3 # Assuming the third state variable is cB in your model output
  
  
  # Add observed data to the data frame
  model_predictions_df$R_obs = oh2_list[['obs_data']]$R_obs
  model_predictions_df$C_obs = oh2_list[['obs_data']]$C_obs
  # model_predictions_df$cB_obs = obs_data$cB_obs
  
  # Melt the data frame for ggplot2
  plot_data <- melt(model_predictions_df, id.vars = "time", variable.name = "Variable", value.name = "Value")
  
  # Create a new column to differentiate between observed and predicted values
  plot_data$Type <- ifelse(grepl("_obs", plot_data$Variable), "Observed", "Predicted")
  
  # Adjust the Variable column to have cleaner names
  plot_data$Variable <- gsub("_obs", "", plot_data$Variable)
  plot_data$Variable <- gsub("V[1-3]", c("R", "C", "cB"), plot_data$Variable)
  
  ggplot(plot_data, aes(x = time, y = Value, color = Type)) + 
    geom_line(aes(linetype = Type), size = 1) + 
    geom_point(aes(shape = Type), size = 3) +
    facet_wrap(~Variable, scales = "free_y") + 
    theme_minimal() + 
    labs(title = "Model Predictions vs. Observations", x = "Time", y = "Value") + 
    scale_color_manual(values = c("Observed" = "red", "Predicted" = "blue")) + 
    scale_shape_manual(values = c("Observed" = 16, "Predicted" = 17))
  
#########
  #st14 ODE model
  # run it in the background after tweaking below
  rstudioapi::jobRunScript(
    path = here::here("ignore/st14_ODEfit.R"),
    name = "st14 ODE",
    workingDir = here::here(),
    importEnv = FALSE,
  )
  
  st14_list = obs_data[['st14']] 
  st14_list$obs_data = st14_list$obs_data |> dplyr::mutate(across(c(R_obs, C_obs, cB_obs), ~.x*10)) 
  
  T_func <- approxfun(st14_list[['temp']]$time, st14_list[['temp']]$tempC, method = "linear", yleft = min(st14_list[['temp']]$tempC), yright = max(st14_list[['temp']]$tempC))
  
  L_func <- approxfun(st14_list[['light']]$time, st14_list[['light']]$light_obs, method = "linear", yleft = min(st14_list[['light']]$light_obs), yright = max(st14_list[['light']]$light_obs))
  
  M_func = approxfun(st14_list[['obs_data']]$time, st14_list[['obs_data']]$M_obs, method = 'linear', yleft = min(st14_list[['obs_data']]$M_obs), yright = max(st14_list[['obs_data']]$M_obs))
  
  ode_system <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      T <- T_func(time)
      L <- L_func(time)
      M <- M_func(time)
      # 
      # r <- r_base * exp(-Er/C_to_overkt(T)) * min(c(1,(1-L/200)))
      K <- K_base + K_scale * L#min(c(1,(1-L/800)))#L
      # 
      # a <- a_base * L#exp(a_T_coef * (T))
      m <- m_base * exp(-Em/C_to_overkt(T))* (M^-alpha) #
      
      # cB <- a_T * R * C / (1 + a_T * h * R)
      # dRdt <- r_TL * R * (1 - R / K_TL) - cB
      # dCdt <- e * cB - m_T * C
      # dRdt <- r_TL * R * (1 - R / K_TL) - a_T * R * C / (1 + a_T * h * R)
      # dCdt <- e * a_T * R * C / (1 + a_T * h * R) - m_T * C
      dRdt <- r* R * (1 - R / K) - a * R * C #/ (1 + a * h * R)
      dCdt <- e * a * R * C - m*C#/ (1 + a * h * R) - m * C
      
      # return(list(c(dRdt, dCdt, cB)))
      return(list(c(dRdt, dCdt)))
    })
  }
  
  loss_function <- function(parameters) {
    # Initial state
    state <- c(R = st14_list[['obs_data']]$R_obs[1], C = st14_list[['obs_data']]$C_obs[1])#, cB = obs_data$cB_obs[1]) 
    # Solve ODE
    out <- ode(y = state, times = st14_list[['temp']]$time, func = ode_system, parms = parameters, method = 'daspk')
    model_out = as.data.frame(out)
    # Extract model predictions
    tt = unlist(st14_list[['obs_data']]$time)
    R_pred <- model_out[which(model_out$time %in% tt), "R"]
    C_pred <- model_out[which(model_out$time %in% tt), "C"]
    # cB_pred <- model_out[, "cB"]
    
    # conform the observations 
    obsDf = st14_list[['obs_data']][which(st14_list[['obs_data']]$time %in% model_out$time),]
    
    # Calculate loss
    loss <- sqrt((sum((obsDf$R_obs - R_pred)^2, na.rm = TRUE) + 
                    sum((obsDf$C_obs - C_pred)^2, na.rm = TRUE)/nrow(obsDf))) #+ 
    # sum((obs_data$cB_obs - cB_pred)^2, na.rm = TRUE)
    
    return(loss)
  }
  
  initial_guesses <- c(r = 28, a = 0.0028,
                       K_base = 400, K_scale = 1, alpha = 0.25,
                       m_base = 0.0095, Em = 0.43, e = 0.20)
  
  st14_model <- ode(y=c(R = st14_list[['obs_data']]$R_obs[1], C = st14_list[['obs_data']]$C_obs[1]), times = st14_list[['temp']]$time, func = ode_system, parms = initial_guesses, method = 'daspk')
  #compare summary stats to observed
  compare_results(st14_model, st14_list[['obs_data']])
  
  plot( st14_model)
  # remove first two years and compare mean and CVs of data and prediction
  
  optim_result <- optim(par = initial_guesses, fn = loss_function, 
                        method = "L-BFGS-B", 
                        lower = rep(0.00001, length(initial_guesses)), 
                        upper = c(60,#r
                                  0.1,#a
                                  500,#K_base
                                  5,#K_scale
                                  2,#alpha
                                  1,#m_base
                                  2,#Em
                                  1),#e,
                        control = list(maxit = 5e2))
  optim_result$par
  saveRDS(optim_result, here::here("data/st14ODE_optim.rds"))
  optim_result <- readRDS(here::here("data/st14ODE_optim.rds"))
  
  st14_model_op <- ode(y=c(R = st14_list[['obs_data']]$R_obs[1], C = st14_list[['obs_data']]$C_obs[1]), times = st14_list[['temp']]$time, func = ode_system, parms = optim_result$par, method = 'daspk')
  
  optim_result2 <- optim(par = optim_result$par, fn = loss_function, 
                         method = "L-BFGS-B", 
                         lower = rep(0.00001, length(initial_guesses)), 
                         upper = c(60,#r
                                   0.1,#a
                                   500,#K_base
                                   5,#K_scale
                                   2,#alpha
                                   1,#m_base
                                   2,#Em
                                   1),#e,
                         control = list(maxit = 5e2))
  
  optim_result2$par
  saveRDS(optim_result2, here::here("data/st14ODE_optim2.rds"))
  optim_result2 <- readRDS(here::here("data/st14ODE_optim2.rds"))
  
  st14CostFun = function(pars){
    out <- as.data.frame( st14_model_op)
    cost <- modCost(model = out, obs = st14_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs'), err = NULL)
    return(modCost(model = out, obs = st14_list[['obs_data']] |> select(time, R = 'R_obs',C = 'C_obs'), err = NULL, cost = cost))}
  
  
  st14Cost = st14CostFun(pars = initial_guesses)
  st14CostFun(initial_guesses)$model
  
  modPars = optim_result$par*1.2
  
  st14Cost2<- function(lpars){
    st14CostFun(c(exp(lpars)))
  }
  # debugonce(modFit)
  fit <- modFit(f = st14Cost2, p = c(log(modPars),
                                    R0 = st14_list[['obs_data']]$R_obs[1],
                                    C0 = st14_list[['obs_data']]$C_obs[1]))
  
  var0 = fit$var_ms_unweighted
  cov0 = summary(fit)
  # parameter fitting using levenberg marquart algorithm
  # initial guess for parameters
  # fitting
  fitval=nls2::nls(par=initial_guesses,fn=loss_function)
  # Check optimized parameters
  
  
  # Assuming `optim_result$par` contains your optimized parameters
  model_predictions <- ode(y = c(R = st14_list[['obs_data']]$R_obs[1],
                                 C = st14_list[['obs_data']]$C_obs[1]), #cB = round(obs_data$cB_obs[1],3)),
                           times = st14_list[['obs_data']]$time, func = ode_system, parms = optim_result$par)
  
  # Convert model output to a data frame
  model_predictions_df <- as.data.frame(model_predictions)
  # model_predictions_df$cB <- model_predictions_df$V3 # Assuming the third state variable is cB in your model output
  
  
  # Add observed data to the data frame
  model_predictions_df$R_obs = st14_list[['obs_data']]$R_obs
  model_predictions_df$C_obs = st14_list[['obs_data']]$C_obs
  # model_predictions_df$cB_obs = obs_data$cB_obs
  
  # Melt the data frame for ggplot2
  plot_data <- melt(model_predictions_df, id.vars = "time", variable.name = "Variable", value.name = "Value")
  
  # Create a new column to differentiate between observed and predicted values
  plot_data$Type <- ifelse(grepl("_obs", plot_data$Variable), "Observed", "Predicted")
  
  # Adjust the Variable column to have cleaner names
  plot_data$Variable <- gsub("_obs", "", plot_data$Variable)
  plot_data$Variable <- gsub("V[1-3]", c("R", "C", "cB"), plot_data$Variable)
  
  ggplot(plot_data, aes(x = time, y = Value, color = Type)) + 
    geom_line(aes(linetype = Type), size = 1) + 
    geom_point(aes(shape = Type), size = 3) +
    facet_wrap(~Variable, scales = "free_y") + 
    theme_minimal() + 
    labs(title = "Model Predictions vs. Observations", x = "Time", y = "Value") + 
    scale_color_manual(values = c("Observed" = "red", "Predicted" = "blue")) + 
    scale_shape_manual(values = c("Observed" = 16, "Predicted" = 17))
  
  
  
  T_func <- approxfun(obs_data$time, obs_data$temp_obs, method = "linear", yleft = min(obs_data$temp_obs), yright = max(obs_data$temp_obs))
  L_func <- approxfun(obs_data$time, obs_data$light_obs, method = "linear", yleft = min(obs_data$light_obs), yright = max(obs_data$light_obs))
  
  
  ode_system <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      T <- T_func(t)
      L <- L_func(t)
      
      r_TL <- r_base * exp(r_T_coef * (T - T_opt))
      K_TL <- K_base * exp(K_T_coef * (T - T_opt) + K_L_coef * (L - L_opt))
      
      a_T <- a_base * exp(a_T_coef * (T - T_opt))
      m_T <- m_base * exp(m_T_coef * (T - T_opt))
      
      cB <- a_T * R * C / (1 + a_T * h * R)
      dRdt <- r_TL * R * (1 - R / K_TL) - cB
      dCdt <- e * cB - m_T * C
      
      return(list(c(dRdt, dCdt, cB)))
    })
  }
  
  loss_function <- function(parameters) {
    # Initial state
    state <- c(R = obs_data$R_obs[1], C = obs_data$C_obs[1], cB = obs_data$cB_obs[1]) 
    
    # Solve ODE
    model_out <- ode(y = state, times = obs_data$time, func = ode_system, parms = parameters)
    
    # Extract model predictions
    R_pred <- model_out[, "R"]
    C_pred <- model_out[, "C"]
    cB_pred <- model_out[, "cB"]
    
    # Calculate loss
    loss <- sum((obs_data$R_obs - R_pred)^2, na.rm = TRUE) + 
      sum((obs_data$C_obs - C_pred)^2, na.rm = TRUE) + 
      sum((obs_data$cB_obs - cB_pred)^2, na.rm = TRUE)
    
    return(loss)
  } 
  
  initial_guesses <- c(r_base = 0.5, r_T_coef = 0.01, r_L_coef = 0.01, 
                       L_opt = 200, T_opt = 25,
                       K_base = 500, K_T_coef = 0.01, K_L_coef = 0.1, 
                       a_base = 0.01, a_T_coef = 0.01, 
                       m_base = 0.3, m_T_coef = 0.01, 
                       e = 0.15, h = 0.1)
  # optim_result <- optim(par = initial_guesses, fn = loss_function, method = "BFGS", hessian = TRUE)
  obs_data = na.omit(obs_data)
  optim_result <- optim(par = initial_guesses, fn = loss_function, method = "L-BFGS-B", lower = rep(0, length(initial_guesses)), upper = c(1,100,100,30,Inf,1e4,Inf,Inf,1,10,1,10,1,Inf), control = list(maxit = 1e5))
  
  # optimization <- function(observed_data) {
  #   # Initial guesses for parameters
  #   initial_guesses <- c(r_base = 0.1, K_base = 10, a_base = 0.01, h = 0.1, e = 0.1, m_base = 0.01,
  #                        r_T_coef = 0.01, r_L_coef = 0.01, K_T_coef = 0.01, K_L_coef = 0.01, 
  #                        a_T_coef = 0.01, m_T_coef = 0.01, T_opt = 20, L_opt = 1000)
  #   
  #   # Optimizing
  #   opt_result <- optim(par = initial_guesses, fn = loss_function, observed_data = observed_data, method = "L-BFGS-B")
  #   
  #   return(opt_result$par)
  # }
  
  # x_R = x_R |>
  #   bind_rows(x_R |> dplyr::mutate(doy = doy + 365)) |>
  #   bind_rows(x_R |> dplyr::mutate(doy = doy + (365*2))) |>
  #   bind_rows(x_R |> dplyr::mutate(doy = doy + (365*3))) |>
  #   bind_rows(x_R |> dplyr::mutate(doy = doy + (365*4))) |>
  #   dplyr::mutate(doy = doy-207) |> dplyr::filter(doy > 1)
  
  # x_C = FW_inputs[[1]][1:1000] |> bind_rows(.id = 'boot_id') |>
  #   group_by(doy) |>
  #   dplyr::summarise(C = mean(bio_mg_m, na.rm = TRUE),
  #                    sigma_C = sd(bio_mg_m, na.rm = TRUE)) |>
  #   na.omit()
  
  # x_C =x_C |> bind_rows(x_C |> dplyr::mutate(doy = doy+365)) |>
  #   bind_rows(x_C |> dplyr::mutate(doy = doy + (365*2)))|>
  #   bind_rows(x_C |> dplyr::mutate(doy = doy + (365*3))) |>
  #   bind_rows(x_C |> dplyr::mutate(doy = doy + (365*4))) |>
  #   dplyr::mutate(doy = doy-207) |> dplyr::filter(doy > 1)
  # x_E = x_E |> bind_rows(x_E |> dplyr::mutate(doy = doy+365)) |>
  #   bind_rows(x_E |> dplyr::mutate(doy = doy + (365*2)))|>
  #   bind_rows(x_E |> dplyr::mutate(doy = doy + (365*3))) |>
  #   bind_rows(x_E |> dplyr::mutate(doy = doy + (365*4))) |>
  #   dplyr::mutate(doy = doy-207) |> dplyr::filter(doy > 1)

 #  # prepare the data
 #  stan_data <- list(T = nrow(stan_dataDf), ts = stan_dataDf$doy,
 #                    R_initial = stan_dataDf$R_obs[1]*10, 
 #                    C_initial = stan_dataDf$C_obs[1]*10,
 #                    cB_initial = stan_dataDf$cB_obs[1]*10,
 #                    R_obs = stan_dataDf$R_obs*10, 
 #                    sigma_R = stan_dataDf$sigma_R*10,
 #                    C_obs = stan_dataDf$C_obs*10, 
 #                    sigma_C = stan_dataDf$C_sigma*10,
 #                    cB_obs = stan_dataDf$cB_obs*10, 
 #                    sigma_cB = stan_dataDf$cB_sigma)#,
 #  # temp = x_E$tempC)
 # genCRmodel = cmdstanr::cmdstan_model(here::here("R/genCR_time.stan"))
 #  gen_CRmodel = stan_model(here::here("R/genCR_time.stan"))
 #  # gen_CRmodel_cmdstan = cmdstan_model(here::here("R/genCR_time.stan"))
 #  #run the model
 #  fitCR <- sampling(object = gen_CRmodel, data = stan_data, chains = 4, iter = 500, cores = 4)
 #  traceplot(fitCR)
 #  summary(fitCR)
 #  # Run Bayesian inference
 #    fit <- stan(
 #      # fit = fit2,
 #      # model_code = readLines(gen_CRmodel),
 #      file =here::here("R/genCR_time.stan"),
 #      data = stan_data,
 #      chains = 4,
 #      cores = 4,
 #      iter = 500
 #    )
  }
  
