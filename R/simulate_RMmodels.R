#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author jimjunker1
#' @export
simulate_RMmodels <- function() {

  # Define the Rosenzweig-MacArthur model function
  rosenzweig_macarthur <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      # consumption dynamics
      # dcB_dt <- 
      
      # Resource dynamics
      dR_dt <- r * R * (1 - (R / K)) - (a*R*C)/(1+a*h*R)
      
      # Consumer dynamics
      dC_dt <- e * (a*R*C)/(1+a*h*R) - m * C
      
      # Return the derivatives
      list(c(dR_dt, dC_dt))
    })
  }
  
  # Function to calculate equilibrium biomass
  equilibrium_biomass <- function(parameters) {
    with(parameters, {
      # Calculate equilibrium biomass
      R_eq <- (m * e * K) / (a * (1 - (m / r)))
      C_eq <- (r * R_eq) / (e * K)
      
      # Return equilibrium biomass values
      list(Resource = R_eq, Consumer = C_eq)
    })
  }
  
  # Example parameter values
  parameters <- list(r = 3, K = 100, a = 0.05, h = 0.5, e = 0.1, m = 0.01)
  
  # Example initial state
  state <- c(R = 10, C = 10)
  
  # Example time vector
  time <- seq(0, 10000, by = 0.1)
  
  # Solve the ODE system
  ode_result <- ode(y = state, times = time, func = rosenzweig_macarthur, parms = parameters)
  plot(ode_result)
  # Calculate equilibrium biomass
  equilibrium_result <- equilibrium_biomass(parameters)
  
  # Print equilibrium biomass
  print(equilibrium_result)
  
  # Function to estimate local stability
  # Function to estimate local stability and maximum eigenvalue
  estimate_stability <- function(R, C, parameters) {
    with(parameters, {
      # Update parameters with current values of R and C
      parameters <- c(parameters, R = R, C = C)
      
      # Calculate the Jacobian matrix
      J <- matrix(0, nrow = 2, ncol = 2)
      J[1, 1] <- r * (1 - (2 * R / K)) - (c * C)
      J[1, 2] <- -c * R
      J[2, 1] <- e * c * C
      J[2, 2] <- e * c * R - m
      
      # Calculate the eigenvalues
      eigenvalues <- eigen(J)$values
      
      # Maximum eigenvalue
      max_eigenvalue <- max(Re(eigenvalues))
      
      # Check the stability
      if(all(Re(eigenvalues) < 0)) {
        stability <- "locally stable"
      } else {
        stability <- "unstable"
      }
      
      # Return the stability result and maximum eigenvalue
      list(stability = stability, max_eigenvalue = max_eigenvalue)
    })
  }
  
  # Example parameter values
  parameters <- list(r = 0.5, K = 100, a = 0.01, h = 2, e = 0.1, m = 0.1)
  
  # Example values of R and C
  R <- 50
  C <- 5
  
  # Estimate local stability and maximum eigenvalue
  stability_result <- estimate_stability(R, C, parameters)
  
  # Print stability result and maximum eigenvalue
  print(stability_result)
  
  return(ode_result)

}
