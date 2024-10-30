#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param RM_models
#' @return
#' @author jimjunker1
#' @export
summarise_RMmodels <- function(obs_data, RM_models) {
  
  # create a sub function that calculates the relevant parameters for mean analysis from fit objects
  
  summarise_model = function(list = NULL, fit = NULL, lpar = TRUE,...){
    lpars = fit$par
    if(lpar == TRUE){
      pars = exp(lpars)
    } else{
      pars = lpars
    }
    T_func <- approxfun(list[['temp']]$time, list[['temp']]$tempC, method = "linear", yleft = min(list[['temp']]$tempC), yright = max(list[['temp']]$tempC))
    
    L_func <- approxfun(list[['light']]$time, list[['light']]$light_obs, method = "linear", yleft = min(list[['light']]$light_obs), yright = max(list[['light']]$light_obs))
    
    M_func = approxfun(list[['obs_data']]$time, list[['obs_data']]$M_obs, method = 'linear', yleft = min(list[['obs_data']]$M_obs), yright = max(list[['obs_data']]$M_obs))
    maxtime = max(list[['light']]$time, na.rm = TRUE)
    
    # ode system with m as output
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
        
        dRdt <- r* R * (1 - R / K) - a * R * C #/ (1 + a * h * R)
        dCdt <- e * a * R * C - m*C#/ (1 + a * h * R) - m * C
        
        # return(list(c(dRdt, dCdt, cB)))
        return(list(c(dRdt, dCdt)))
      })
    }
    # fit the ODE with the fit parameters
    model <- ode(y=c(R = list[['obs_data']]$R_obs[1], C = list[['obs_data']]$C_obs[1]), times = list[['temp']]$time, func = ode_system, parms = pars, method = 'daspk')
    
    Kvec = pars['K_base'] + pars['K_scale'] * L_func(1:maxtime)
    Kmean = mean(Kvec, na.rm = TRUE)
    
    mVec = pars['m_base'] * exp(-pars['Em']/C_to_overkt(T_func(1:maxtime))) * (M_func(1:maxtime)^-pars['alpha'])
    mMean = mean(mVec, na.rm = TRUE)
    
    # calculate the equilibrium at the mean values
    equilibrium <- function(state, params) {

      R <- state[1]
      C <- state[2]
      
      f1 <- params['r'] * R * (1 - R / params['K']) - params['a'] * R * C
      f2 <- params['e'] * params['a'] * R * C - params['m'] * C
      
      return(c(f1, f2))
    }
    
    # Find equilibrium numerically
    state <- c(R = mean(list[['obs_data']]$R_obs, na.rm = TRUE), C = mean(list[['obs_data']]$C_obs, na.rm = TRUE))
    
    params <- c(a = unname(pars['a']), r = unname(pars['r']), K = Kmean, e = unname(pars['e']), m = mMean)
    equilibriumPoint <- multiroot(f = equilibrium, start = state, parms = params)
    ## estimate the jacobian matrix
    jacobianMatrix <- function(R, C, params) {
      a <- unname(params['a'])
      J <- matrix(c(
        params['r'] * (1 - 2 * R / params['K']) - params['a'] * C, -params['a'] * R,
        params['e'] * params['a'] * C, params['e'] * params['a'] * R - params['m']
      ), nrow = 2, byrow = TRUE)
      
      return(J)
    }
    # calculated the Jacobian matrix at equlibrium
    J <- jacobianMatrix(equilibriumPoint$root[1], equilibriumPoint$root[2], params)
    
    eigenvalues <- eigen(J)$values
    realParts <- Re(eigenvalues)
    return(list(equilibriumPoint = equilibriumPoint,
                realParts = realParts,
                a = pars['a'],
                e = pars['e'],
                K = Kmean,
                m = mMean))
  }
  
  estimateHopfBifurcation <- function(r, K, e, m, a) {
    # Define the Jacobian matrix for the Rosenzweig-MacArthur model
    jacobianMatrix <- function(R, C, params) {
      a <- params$a
      J <- matrix(c(
        r * (1 - 2 * R / K) - a * C, -a * R,
        e * a * C, e * a * R - m
      ), nrow = 2, byrow = TRUE)
      
      return(J)
      
      # Function to find a critical point where the real part of the dominant eigenvalue is zero
      findCriticalPoint <- function(a, r, K, e, m) {
        # Guess initial values for R and C
        R_guess <- K / 2
        C_guess <- r / (a * e)
        
        # Find equilibrium numerically
        state <- c(R = R_guess, C = C_guess)
        params <- c(a = a, r = r, K = K, e = e, m = m)
        equilibriumPoint <- multiroot(f = equilibrium, start = state, parameters = params)
        
        # Calculate the Jacobian at the equilibrium
        J <- jacobianMatrix(equilibriumPoint$root[1], equilibriumPoint$root[2], params)
        
        # Find the eigenvalues
        eigenvalues <- eigen(J)$values
        
        # Check for Hopf bifurcation condition (dominant eigenvalue real part = 0)
        realParts <- Re(eigenvalues)
        if(any(realParts == 0)) {
          return(list(a = a, R = equilibriumPoint$root[1], C = equilibriumPoint$root[2], status = "Possible Hopf Bifurcation"))
        } else {
          return(list(status = "No Hopf Bifurcation"))
        }
      }
    }}
  
opListNames = list.files(path = here::here('data'), "._fmeFit_op.rds") %>% 
    gsub("(.)_.*_op.rds","\\1", x = .)

opFitList = list.files(path = here::here("data/"), "._fmeFit_op.rds", full.names = TRUE) |> map(\(x) readRDS(x)) |> setNames(opListNames) |> rlist::list.subset(names(stream_order_list)) %>% map(~.x %>% pluck('fit'))

initListNames = list.files(path = here::here('data'), "._fmeFit_init.rds") %>% 
  gsub("(.)_.*_init.rds","\\1", x = .)

initFitList = list.files(path = here::here("data/"), "._fmeFit_init.rds", full.names = TRUE) |> map(\(x) readRDS(x)) |> setNames(initListNames) |> rlist::list.subset(names(stream_order_list)) %>% map(~.x %>% pluck('fit'))

opFitSummList = map2(opFitList, obs_data, 
                     ~summarise_model(list = .y, fit = .x))
   
initFitSummList = map2(initFitList, obs_data, 
                       ~summarise_model(list = .y, fit = .x))

initAggDf = initFitSummList |> map(\(x) rho = unname(x[['e']])*unname(x[['a']]*unname(x[['K']])/unname(x[['m']]))) |> unlist() |> bind_cols(as.numeric(stream_temp_labels),y = _) |> setNames(c('temperature','rho'))

initAggTemp_gam <- mgcv::gam(rho ~ temperature, data = initAggDf)
summary(initAggTemp_gam)

png(here::here("data/Bcrplot.png"), res = 400, height = 5, width = 7, units = 'in')
initAggDf |> ggplot()+
  geom_point(aes(x = temperature, y = rho), size =5)+
  scale_y_continuous(name = expression(beta[CR]), limits = c(0,45), expand = c(0,0))+
  scale_x_continuous(name = expression('Temperature ('*degree*C*")"), limits = c(0,30), expand = c(0,0))+
  scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos])+
  theme_minimal()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20))+
  NULL
dev.off()

opVarDf = opFitSummList |> 
  map(\(x){
    e = unname(x[['e']])
    a = unname(x[['a']])
    k = unname(x[['K']])
    m = unname(x[['m']])
    return(data.frame(e = e,
                      a = a,
                      k = k,
                      m = m))
  }) |>
  bind_rows() |>
  bind_cols(temp_C = as.numeric(stream_temp_labels),y = _)|>
  pivot_longer(-temp_C, names_to = 'variable', values_to = 'value')

opVarDf %>% 
  ggplot()+
  geom_point(aes(x = temp_C, y = value, color = variable), size = 2)+
  geom_smooth(aes(x = temp_C , y = value, color = variable), se = FALSE,
              method = 'lm')+
  facet_wrap(~variable, ncol = 2, scales = 'free_y')


x = map2(list(st7_list, oh2_list), y = as.list(st7_fmeFit_init$fit, oh2_fmeFit_init$fit), ~summarise_model(list = .x, fit = .y))
   
}


   

   ###
   
   
   
   
   
   
   
   
#    # Function to find a critical point where the real part of the dominant eigenvalue is zero
#    findCriticalPoint <- function(a, r, K, e, m) {
#      # Guess initial values for R and C
#      R_guess <- K / 2
#      C_guess <- r / (a * e)
#      
#      # Solve for equilibrium
#      equilibrium <- function(state, params) {
#        R <- state[1]
#        C <- state[2]
#        
#        f1 <- r * R * (1 - R / K) - a * R * C
#        f2 <- e * a * R * C - m * C
#        
#        return(c(f1, f2))
#      }
#      
#      # Find equilibrium numerically
#      state <- c(R = R_guess, C = C_guess)
#      params <- c(a = a, r = r, K = K, e = e, m = m)
#      equilibriumPoint <- multiroot(f = equilibrium, start = state, parameters = params)
#      
#      # Calculate the Jacobian at the equilibrium
#      J <- jacobianMatrix(equilibriumPoint$root[1], equilibriumPoint$root[2], params)
#      
#      # Find the eigenvalues
#      eigenvalues <- eigen(J)$values
#      
#      # Check for Hopf bifurcation condition (dominant eigenvalue real part = 0)
#      realParts <- Re(eigenvalues)
#      if(any(realParts == 0)) {
#        return(list(a = a, R = equilibriumPoint$root[1], C = equilibriumPoint$root[2], status = "Possible Hopf Bifurcation"))
#      } else {
#        return(list(status = "No Hopf Bifurcation"))
#      }
#    }
#    
#    # Example usage
#    a_range <- seq(0.1, 1, length.out = 100)
#    results <- lapply(a_range, function(a) findCriticalPoint(a, r, K, e, m))
#    
#    return(results)
#  }
#  
#  # Example call to the function (adjust parameters as needed)
#  results <- estimateHopfBifurcation(r = 1, K = 100, e = 0.1, m = 0.1)
#  
# }
