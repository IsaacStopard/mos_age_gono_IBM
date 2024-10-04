
functions{
  // function to integrate the number of mosquitoes
  real ode_int(real o_, real mu_, real m_0, real t_){
    real result = exp(-mu_ * t_) * (m_0 * mu_ + o_ * (exp(mu_ * t_) - 1)) / mu_;
    return result;
  }
  
}

data {
  int<lower = 0> max_t; // max time
  
  int<lower = 0> N; // number of mosquito counts observations
  int t[N]; // mosquito count times
  int M[N]; // mosquito counts
  
  real<lower = 0> mu; // per mosquito mortality rate
}

parameters {
  real<lower = 0> m0;
  real<lower = 0> kappa; // overdispersion parameters
  
  real<lower = 0> c;
  real<lower = -1, upper = 1> h; // autocorrelation
  real<lower=0> sigma; //
  
  real<lower = 0> o[max_t - 1]; // daily emergence rate at the population level
}

transformed parameters{
  real m[max_t]; // predicted numbers of mosquitoes
  
  m[1] = m0;
  for(i in 2:max_t){
    m[i] = ode_int(o[i - 1], mu, m[(i - 1)], 1.0);
  }
}

model {
  //priors
  m0 ~ normal(15.0, 2.5);
  kappa ~ normal(1, 2.5);
  c ~ normal(0.132 * 15.0, 2.5);
  h ~ normal(0, 1);
  sigma ~ normal(0.5, 2.5);
  
  o[1] ~ normal(mu, sqrt(1 - h^2) * sigma);
  for(i in 2:(max_t - 1)){
    o[i] ~ normal(c + h * (o[i - 1] - c), sqrt(1 - h^2) * sigma);
  }
  
  // likelihood
  for(i in 1:N){
    target += neg_binomial_2_lpmf(M[i] | m[t[i]], kappa);
  }
}

generated quantities{
  real pred_m[max_t];
  
  for(i in 1:max_t){
    pred_m[i] = neg_binomial_2_rng(m[i], kappa);
  }
}

