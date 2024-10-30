# # Load the necessary package
# library(deSolve)
# 
# # Set the maximum time for the simulation
# tmax <- 20000
# 
# # Define the parameters as a named list
# param <- list(
#   a1 = 1.3,
#   r1 = 10.0,
#   K1 = 5.0,
#   d1 = 0.20,
#   b1 = 1.0,
#   e = 0.7,
#   Amp = 0.5,
#   Ps = 0.0,
#   p = 1
# )
# 
# # Define the forcing function
# Kforce <- function(t, Amp, p, Ps) {
#   return(Amp * sin((p * 2 * pi) * t + Ps * pi))
# }
# 
# # Define the differential equations for CR
# CR <- function(times, y, params) {
#   with(as.list(c(y, params)), {
#     # Calculate Kforce
#     Kforce_value <- Kforce(times, Amp, p, Ps)
#     
#     # Differential equations using f1 and f2
#     dR1 <- R1 * (r1 * (1 - R1 / (K1 + Kforce_value)) - a1 * C1 / (b1 + R1))
#     dC1 <- C1 * (e * a1 * R1 / (b1 + R1) - d1)
#     
#     # Return the rate of change
#     list(c(dR1, dC1))
#   })
# }
# 
# # Define the differential equations for CRu
# CRu <- function(times, y, params) {
#   with(as.list(c(y, params)), {
#     # Differential equations using f1u and f2u
#     dR1 <- R1 * (r1 * (1 - R1 / K1) - a1 * C1 / (b1 + R1))
#     dC1 <- C1 * (e * a1 * R1 / (b1 + R1) - d1)
#     
#     # Return the rate of change
#     list(c(dR1, dC1))
#   })
# }
# 
# # Function to solve CR equations
# sol <- function(R0, C0, Kval, rval, eval, freq, param) {
#   
#   # Update the parameters with the passed values
#   param <- modifyList(param, list(K1 = Kval, r1 = rval, e = eval, p = freq))
#   
#   # Initial conditions
#   y0 <- c(R1 = R0, C1 = C0)
#   
#   # Time vector
#   times <- seq(0, tmax, by = 0.01)
#   
#   # Solve the differential equations
#   out <- ode(y = y0, times = times, func = CR, parms = param,
#              atol = 1e-7, rtol = 1e-7, maxsteps = 100000000)
#   
#   return(out)
# }
# 
# # Function to solve CRu equations
# solu <- function(R0, C0, Kval, rval, eval, param) {
#   
#   # Update the parameters with the passed values
#   param <- modifyList(param, list(K1 = Kval, r1 = rval, e = eval))
#   
#   # Initial conditions
#   y0 <- c(R1 = R0, C1 = C0)
#   
#   # Time vector
#   times <- seq(0, tmax, by = 0.01)
#   
#   # Solve the differential equations
#   out <- ode(y = y0, times = times, func = CRu, parms = param,
#              atol = 1e-7, rtol = 1e-7, maxsteps = 100000000)
#   
#   return(out)
# }
# 
# # Example usage:
# # Solve CR equations with specific parameters
# result_sol <- sol(R0 = 1, C0 = 1, Kval = 5.0, rval = 10.0, eval = 0.7, freq = 1, param = param)
# 
# # Solve CRu equations with specific parameters
# result_solu <- solu(R0 = 1, C0 = 1, Kval = 5.0, rval = 10.0, eval = 0.7, param = param)
# 
# # The result can be printed or plotted
# print(result_sol)
# print(result_solu)
# 
# 
# ### second part----
# # Load necessary libraries
# library(deSolve)
# library(pracma)   # for Jacobian and Eigenvalues
# library(parallel) # for parallel processing
# library(ggplot2)  # for plotting
# 
# # Define the parameters as a named list
# param <- list(
#   a1 = 1.3,
#   r1 = 10.0,
#   K1 = 5.0,
#   d1 = 0.20,
#   b1 = 1.0,
#   e = 0.7,
#   Amp = 0.5,
#   Ps = 0.0,
#   p = 1
# )
# 
# # Define the functions f1u and f2u
# f1u <- function(R1, C1, params) {
#   with(as.list(params), {
#     return(R1 * (r1 * (1 - R1 / K1) - a1 * C1 / (b1 + R1)))
#   })
# }
# 
# f2u <- function(R1, C1, params) {
#   with(as.list(params), {
#     return(C1 * (e * a1 * R1 / (b1 + R1) - d1))
#   })
# }
# 
# # Define the Jacobian
# jacobian <- function(R1, C1, params) {
#   J <- jacobian(function(x) c(f1u(x[1], x[2], params), f2u(x[1], x[2], params)), c(R1, C1))
#   return(J)
# }
# 
# # Solve for equilibrium points
# solve_equilibrium <- function(params) {
#   eqs <- uniroot.all(function(R1) {
#     f1 <- f1u(R1, 1, params)
#     f2 <- f2u(R1, 1, params)
#     return(c(f1, f2))
#   }, interval = c(0, 10))
#   return(eqs)
# }
# 
# # Simplify and calculate Eigenvalues at equilibrium
# eigenvalues_at_eq <- function(eqs, params) {
#   eigenvals <- sapply(eqs, function(eq) {
#     J <- jacobian(eq, 1, params)
#     return(eigen(J)$values)
#   })
#   return(eigenvals)
# }
# 
# # Compute bifurcation parameters
# bifurcation_params <- function(params) {
#   with(as.list(params), {
#     b2 <- (a1 * d1 * e * (-b1 + K1) * r1 - d1^2 * (b1 + K1) * r1) / (2 * a1 * e * (-d1 + a1 * e) * K1)
#     return(solve(b2 == 0, K1))
#   })
# }
# 
# # Compute Hopf bifurcation value
# KHopf <- function(params) {
#   with(as.list(params), {
#     return(-((b1 * (d1 + a1 * e)) / (d1 - a1 * e)))
#   })
# }
# 
# # Compute eigenvalues over a range of K1 values
# compute_eigenvalues <- function(Kvals, params) {
#   eigenvalues <- mclapply(Kvals, function(Kval) {
#     params$K1 <- Kval
#     eqs <- solve_equilibrium(params)
#     max_eigen <- max(Re(eigenvalues_at_eq(eqs, params)))
#     return(c(Kval, max_eigen))
#   }, mc.cores = detectCores() - 1)
#   return(do.call(rbind, eigenvalues))
# }
# 
# # Plot the results
# plot_results <- function(eigenvalues, Keigs) {
#   p1 <- ggplot(data = as.data.frame(eigenvalues), aes(x = V1, y = V2)) +
#     geom_point(color = "black", size = 0.5) +
#     labs(x = "K", y = "Max Eigenvalue") +
#     theme_minimal()
#   
#   p2 <- ggplot(data = as.data.frame(Keigs), aes(x = V1, y = V2)) +
#     geom_line(color = "black") +
#     geom_ribbon(aes(ymin = 0, ymax = V2), fill = "green", alpha = 0.3) +
#     labs(x = "Kmean", y = "Max Eigenvalue") +
#     theme_minimal()
#   
#   return(gridExtra::grid.arrange(p1, p2, ncol = 2))
# }
# 
# # Calculate Keigs
# Keigs <- function(rval, eval, params) {
#   params$r1 <- rval
#   params$e <- eval
#   Kvals <- seq(KHopf(params) - params$Amp, KHopf(params) + params$Amp, by = 0.01)
#   return(compute_eigenvalues(Kvals, params))
# }
# 
# # Run the analysis
# params <- modifyList(param, list(r1 = 5, e = 0.5))
# eigenvalues <- compute_eigenvalues(seq(0.1, 3.5, by = 0.01), params)
# Keigs_result <- Keigs(5, 0.5, params)
# 
# # Plot the results
# plot_results(eigenvalues, Keigs_result)
# 
# # Calculate and plot Keigdiff
# Keigdiff <- function(rval, eval, params) {
#   params$r1 <- rval
#   params$e <- eval
#   Kmeans <- seq(KHopf(params), KHopf(params) + 2 * params$Amp, by = 0.05)
#   diffs <- sapply(Kmeans, function(Kmean) {
#     Kvals <- seq(Kmean - params$Amp, Kmean + params$Amp, by = 0.01)
#     sum_diff <- sum(sapply(Kvals, function(Kval) {
#       params$K1 <- Kval
#       eqs <- solve_equilibrium(params)
#       max_eig_val <- max(Re(eigenvalues_at_eq(eqs, params)))
#       max_eig_val - max(Re(eigenvalues_at_eq(eqs, params)))
#     }))
#     return(c(Kmean, sum_diff))
#   })
#   return(do.call(rbind, diffs))
# }
# 
# # Plot Keigdiff
# plot_keigdiff <- function(diffs) {
#   p <- ggplot(data = as.data.frame(diffs), aes(x = V1, y = V2)) +
#     geom_point(color = "black", size = 0.5) +
#     labs(x = "Kmean", y = "Integral of λ - λ mean") +
#     theme_minimal()
#   return(p)
# }
# 
# # Calculate Keigdiff and plot
# Keigdiff_results <- Keigdiff(1.5, 0.7, params)
# p_diff1 <- plot_keigdiff(Keigdiff_results)
# Keigdiff_results <- Keigdiff(5, 0.7, params)
# p_diff2 <- plot_keigdiff(Keigdiff_results)
# Keigdiff_results <- Keigdiff(10, 0.7, params)
# p_diff3 <- plot_keigdiff(Keigdiff_results)
# Keigdiff_results <- Keigdiff(20, 0.7, params)
# p_diff4 <- plot_keigdiff(Keigdiff_results)
# 
# # Display all plots in a grid
# gridExtra::grid.arrange(p_diff1, p_diff2, p_diff3, p_diff4, ncol = 2)
# 
