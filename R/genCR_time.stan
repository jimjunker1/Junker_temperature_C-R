functions {
  vector C_R(real t,           // Time
             vector y,         // State variables: y[1] = R, y[2] = C
             vector theta     // Parameters
                   // Real data (temperature, light)
             ) {     // Int data (empty here)
    // observed states
    real R = y[1];
    real C = y[2];
    real cB = y[3];
    // parameters
    real a = theta[1];
    real h = theta[2];
    real r = theta[3];
    real K = theta[4];
    real e = theta[5];
    real m = theta[6];
    // real data

    real dcB_dt = a*R*C/(1+a*h*R);    // Biomass consumption
    real dR_dt = r*R*(1 - R/K) - cB;  // Prey dynamics
    real dC_dt = (e * cB) - (m * C);  // Predator dynamics
    
    if(is_nan(dR_dt) || is_nan(dC_dt) || is_nan(dcB_dt)) {
    reject("NaN detected; t = ", t, "; derivs = ", [dR_dt, dC_dt, dcB_dt],
              "; y = ", y, "; theta = ", theta);}
    return {dR_dt, dC_dt, dcB_dt};
  }
}

data {
  int<lower=1> T;             // Number of time steps
  real t0;
  array[T] real<lower=0> ts;        // Time values
  vector[3] y_inits;
  // real<lower=0> R_initial;    // Initial resource population
  // real<lower=0> C_initial;    // Initial consumer population
  // real<lower=0> cB_initial;   // Initial consumption
  matrix[3,T] y_obs;
  // real<lower=0> C_obs[T];     // Observed mean consumer population data
  // real<lower=0> R_obs[T];     // Observed mean resource population data
  // real<lower=0> cB_obs[T];    // Observed Biomass specific consumption
  matrix[3,T] sigma;
  // real<lower=0> sigma_C[T];   // Observed sd consumer population data
  // real<lower=0> sigma_R[T];   // Observed sd resources population data
  // real<lower=0> sigma_cB[T];  // Observed sd consumption
  // real<lower=0> temp[T];      // observed temperature data
  // real<lower=0> light[T];     // observed light data
}
// transformed data{
//   real x_r[0];
//   int x_i[0];
// }
parameters {
  real<lower=0,upper=1> a;         // attack/clearance rate parameter
  real<lower=0> h;         // handling time
  real<lower=0> r;         // prey growth rate
  real<lower=0> K;         // Carrying capacity
  real<lower=0,upper=1> m;         // Consumer mortality rate
  real<lower=0,upper=1> e; // gross growth efficiency
  // real<lower=0> sigma;
}

model {
    real y_init[3] = {R_initial, C_initial, cB_initial};
    real theta[6] = {a, h, r, K, e, m};  // Parameters
    real y_sim[T,3];  // Simulated population dynamics

  // Priors
  r ~ lognormal(1, 1);
  K ~ lognormal(10, 3);
  a ~ beta(1, 1);
  h ~ lognormal(0, 1);
  e ~ beta(1, 4);
  m ~ beta(1, 1);

  // Solve ODE
  // y_sim = integrate_ode_bdf(C_R, y_init, 0, ts, theta, x_r, x_i);
  // y_sim = integrate_ode_rk45(C_R, y_init, 0, ts, theta, x_r, x_i);
   array[T] vector[3] y_sim = ode_rk45(C_R, y_init, 0, ts, theta, x_r, x_i);
  matrix[3,T] y_sim_mat; // define a matrix version
  print(y_sim);

// Likelihood
for (t in 1:T) {
        y_sim_mat[,t] = y_sim[t];
  // R_obs[t] ~ normal(y_sim[t, 1], sigma_R[t]);
  // C_obs[t] ~ normal(y_sim[t, 2], sigma_C[t]);
  // cB_obs[t] ~ normal(y_sim[t,3], sigma_cB[t]);    
  }
      to_vector(y_obs) ~ normal(to_vector(y_sim_mat), to_vector(sigma));

}
// model {
//     // define inits & such
//     real y_init[3] = {R_initial, C_initial, cB_initial};
//     real theta[6] = {a, h, r, K, e, m};  // Parameters
//     
//   // Priors
//   r ~ lognormal(1, 1);
//   K ~ lognormal(10, 3);
//   a ~ beta(1, 1);
//   h ~ lognormal(0, 1);
//   e ~ beta(1, 4);
//   m ~ beta(1, 1);
//   // sigma ~ exponential(0.5);
//   
//   array[T] vector[3] y_sim = integrate_ode_rk45(C_R, y_init, 0, ts, theta, x_r, x_i); 
//     matrix[3,T] y_sim_mat; // define a matrix version
// 
// // Likelihood
//     for(t in 1:T) {
//         y_sim_mat[,t] = y_sim[t];
//     }
//     // vectorized likelihood, where y_obs and sigma are matrices defined in data block
//     to_vector(y_obs) ~ normal(to_vector(y_sim_mat), to_vector(sigma));
// }

generated quantities {
  real y_init[3] = {R_initial, C_initial, cB_initial};
  real theta[6];
  real y_pred[T, 3];
  
  // Parameters
  theta[1] = a;            // Saturation parameter
  theta[2] = h;            // Half-saturation constant (for resource)
  theta[3] = r;            // Base prey growth rate
  theta[4] = K;            // Carrying capacity
  theta[5] = e;            // Error term for consumer population
  theta[6] = m;            // Consumer mortality rate

  // Solve ODE for generated quantities
    // y_pred = integrate_ode_bdf(C_R, y_init, 0, ts, theta, x_r, x_i);
    y_pred = integrate_ode_rk45(C_R, y_init, 0, ts, theta, x_r, x_i);

}
