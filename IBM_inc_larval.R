# discrete time model
library(malariasimulation); library(foresite)

orderly2::orderly_shared_resource("default_values.R", 
                                  "IBM_adult_only.R")

source("shared/default_values.R")

##################################
##### seasonality parameters #####
##################################

season <- foresite::BFA$seasonality |> filter(name_1 == "Cascades") |> select(g0, g1, g2, g3, h1, h2, h3) |> unlist()

# baseline mortality
mum <- mu
rm(list = c("mu"))

parameters <- list(
  "g0" = season["g0"],
  "g" = season[paste0("g", seq(1, 3))],
  "h" = season[paste0("h", seq(1, 3))],
  "rainfall_floor" = 0.0001,
  "mum" = mum,
  "foraging_time" = malariasimulation::gamb_params$foraging_time,
  "blood_meal_rates" = gono_rate,
  "me" = 0.0338,
  "ml" = 0.0348,
  "gamma" = 13.25,
  "mup" = 0.249,
  # transition probabilities
  "del" = 6.64,
  "dl" = 3.72,
  "dpl" = 0.643,
  "beta" = 21.2, # maximum number of eggs laid per day
  "blood_meal_rates" = 1/3,
  "m" = m_in, # initial numbers of mosquitoes
  "Q0" = malariasimulation::gamb_params$Q0
)

#####################
##### functions #####
#####################

### the functions for the fourier series
rainfall <- function(t, parameters){
  
  # so between 0 and 1
  ts <- t / 365
  tf <- ts - floor(ts)
  
  result <- parameters$g0
  for (i in 1:length(parameters$g)){
    
    result <- result + (parameters$g[i] * cos(2 * pi * tf * i)) +
      parameters$h[i] * sin(2 * pi * tf * i)
    
  }
  
  return(max(result, parameters$rainfall_floor))
}

### functions to calculate the carrying capacity
eggs_laid <- function(beta, mu, f){
  eov <- beta / mu * (exp(mu / f) - 1)
  return(eov * mu * exp(-mu / f) / (1 - exp(-mu / f)))
}



calculate_omega <- function(parameters,
                            mu,
                            f
                            ) {
  
  sub_omega <- parameters$gamma * parameters$ml / parameters$me - (
    parameters$del / parameters$dl) + ((parameters$gamma - 1) * parameters$ml * parameters$del)
  
  mum <- parameters$mum
  
  beta <- eggs_laid(parameters$beta,
    mum,
    parameters$blood_meal_rates
  )
  
  return(
    -.5 * sub_omega + sqrt(
    .25 * sub_omega^2 +
      .5 * parameters$gamma * beta * parameters$ml * parameters$del /
      (parameters$me * mum * parameters$dl * (
        1. + parameters$dpl * parameters$mup
      ))
  )
  )
}

calculate_K0 <- function(parameters, m) {
  omega <- calculate_omega(parameters)
  
  return(m * 2 * parameters$dl * parameters$mum * (
    1. + parameters$dpl * parameters$mup
  ) * parameters$gamma * (omega + 1) / (
    omega / (parameters$ml * parameters$del) - (
      1. / (parameters$ml * parameters$dl)
    ) - 1.
  )
  )
}

calculate_R_bar <- function(parameters) {
  mean(sapply(1:365, function(t){
    rainfall(t, 
             parameters)}
    )
    )
}


carrying_capacity <- function(timestep,
                              parameters
                              ){
  if (parameters$model_seasonality){
    r = rainfall(timestep, parameters);
    return(K0 * r / R_bar);
  } else{
    return(K0)
  }
}



### mortality rate calculations
get_gonotrophic_cycle <- function(parameters) {
  f <- parameters$blood_meal_rates
  gonotrophic_cycle <- 1 / f - parameters$foraging_time
}

# mortality rate varying as function of time
adult_death_rate <- function(t, parameters) {
  # initial mortality rate
  if(parameters$int_on & t > parameters$t_on){
    parameters$mum * parameters$step_percent
  } else{
    parameters$mum
  }
}


# calculates the probability of an event happening from a constant rate
m_rate_E <- function(me, E, L, K){
  me * (1 + (E + L) / K)
}

m_rate_L <- function(ml, E, L, K){
  ml * (1 + gamma * (E + L) / K)
}

# calculates the probability of an event happening from a constant rate
constant_rate_to_prob <- function(rate){
  1 - exp(-rate)
}

parameters <- append(parameters, c("m_prob_P" = constant_rate_to_prob(parameters$mup),
                                   "d_prob_E" = constant_rate_to_prob(1/parameters$del),
                                   "d_prob_L" = constant_rate_to_prob(1/parameters$dl),
                                   "d_prob_P" = constant_rate_to_prob(1/parameters$dp)))

parameters <- append(parameters, c("INT_on" = TRUE,
                                       "t_on" = 365,
                                       "step_percent" = 0.5))

# probability early larval stages are alive
parameters <- append(parameters, 
                     c("R_bar" = calculate_R_bar(parameters), 
                       "K0" = calculate_K0(parameters, m = m_in)))
run_IBM_larval <- function(parameters){
  
  # putting the parameters in the environment
  # states
  # A - adult female
  # E - eggs
  # L - larvae
  # P - pupae
  
  mos_states <- c("A", "E", "L", "P")
  
  mos_states_t0 <- rep("A", N0) # creating the initial states of each mosquito
  
  # initialising the numbers in each 
  mos_pop <- CategoricalVariable$new(categories = mos_states, 
                                     initial_values = mos_states_t0)
  
  # age variable to keep track of the ages
  # all mosquitoes start at age zero
  age <- DoubleVariable$new(initial_values = rep(0, N0)) 
  
  ### processes
  transition_process <- function(t){
    # getting index positions of different life stages
    i_E <- mos_pop$get_index_of("E")$to_vector()
    i_L <- mos_pop$get_index_of("L")$to_vector()
    i_P <- mos_pop$get_index_of("P")$to_vector()
    i_A <- mos_pop$get_index_of("A")$to_vector()
    
    pop_E <- mos_pop$get_size_of("E")
    pop_L <- mos_pop$get_size_of("L")
    pop_P <- mos_pop$get_size_of("P")
    pop_A <- mos_pop$get_size_of("A")
    
    Kt <- carrying_capacity(timestep = t, parameters = parameters)
    
    m_prob_Et <- constant_rate_to_prob(
      m_rate_E(E = pop_E, L = pop_L, K = Kt, me = parameters$me)
      )
    
    m_prob_Lt <- constant_rate_to_prob(
      m_rate_L(E = pop_E, L = pop_L, K = Kt, ml = parameters$ml, gamma = parameters$gamma)
    )
    
    m_prob_A <- constant_rate_to_prob(
      adult_death_rate(t = t, parameters = parameters)
    )
    
    # deaths
    deaths_E <- sample(i_E, rbinom(1, pop_E, m_prob_Et), replace = FALSE)
    deaths_L <- sample(i_L, rbinom(1, pop_L, m_prob_Lt), replace = FALSE)
    deaths_P <- sample(i_P, rbinom(1, pop_P, parameters$m_prob_P), replace = FALSE)
    deaths_A <- sample(i_A, rbinom(1, pop_A, m_prob_A), replace = FALSE)
    
    # transitions
    transitions_E <- sample(i_E[i_E %in% deaths_E == 0], rbinom(1, pop_E, parameters$d_prob_E), replace = FALSE)
    transitions_L <- sample(i_L[i_L %in% deaths_L == 0], rbinom(1, pop_L, parameters$d_prob_L), replace = FALSE)
    transitions_P <- sample(i_P[i_P %in% deaths_P == 0], rbinom(1, pop_P, parameters$d_prob_P), replace = FALSE)
    
    transitions_P_male <- sample(transitions_P, rbinom(1, length(transitions_P), 1/2), replace = FALSE)
    transitions_P_female <- transitions_P[transitions_P %in% transitions_P_male == 0]
    
    # deaths
    mos_pop$queue_shrink(index = c(deaths_A, deaths_E, deaths_L, deaths_P, transitions_P_male))
    age$queue_shrink(c(deaths_A, deaths_E, deaths_L, deaths_P, transitions_P_male))
    
    # transitions
    mos_pop$queue_update(value = "L", index = transitions_E)
    age$queue_update(value = 0, index = transitions_E)
    
    mos_pop$queue_update(value = "P", index = transitions_L)
    age$queue_update(value = 0, index = transitions_L)
    
    mos_pop$queue_update(value = "A", index = transitions_P_female)
    age$queue_update(value = 0, index = transitions_P_female)
    
    # age process
    age_inds <- c(i_E, i_L, i_P, i_A)[c(i_E, i_L, i_P, i_A) %in% c(deaths_E, transitions_E, deaths_L, transitions_L, deaths_P, transitions_P_female, deaths_A) == 0]
    mos_ages <- age$get_values() 
    mos_ages[age_inds] <- mos_ages[age_inds] + dt 
    age$queue_update(mos_ages, index = NULL)
    
    # birth process
    
    beta_larval * pop_A
  }
  
  ##### rendering
  mos_render <- Render$new(timesteps = steps)
  
  state_render <- categorical_count_renderer_process(
    renderer = mos_render,
    variable = mos_pop,
    categories = mos_states
  )
  
  integer_count_renderer <- function(renderer, variable, categories){
    stopifnot(inherits(variable, "IntegerVariable"))
    stopifnot(inherits(renderer, "Render"))
    function(t){
      for (c in categories){
        renderer$render(paste0(c), variable$get_size_of(c), t)
      }
    }
  }
  
  double_count_renderer <- function(renderer, variable, categories){
    stopifnot(inherits(variable, "DoubleVariable"))
    stopifnot(inherits(renderer, "Render"))
    function(t){
      for (c in categories){
        renderer$render(paste0(c), sum(round(variable$get_values(), digits = log10(dt)) == c), t)
      }
    }
  }
  
  double_renderer <- function(renderer, variable, category){
    stopifnot(inherits(variable, "DoubleVariable"))
    stopifnot(inherits(renderer, "Render"))
    function(t){
      renderer$render(paste0(category), variable$get_values(index = 1), t)
    }
  }
  
  age_render <- double_count_renderer(
    renderer = mos_render,
    variable = age,
    categories = seq(0, 150, dt)
  )
  
  # running the simulation
  simulation_loop(
    variables = list(tot_lag,
                     mos_pop,
                     age,
                     tot_deaths,
                     tot_births),
    events = list(
      parity_event
    ),
    processes = list(
      death_process,
      birth_process,
      parity_process,
      age_process,
      age_render,
      state_render,
      lag_render,
      death_render,
      birth_render
    ),
    timesteps = steps
  )
  
  states <- mos_render$to_dataframe() |> 
    mutate(state_tot = NP_count + P_count,
           age_tot = rowSums(pick(X0:X150)),
           age_10_plus = rowSums(pick(X10:X150)),
           prop_parous = P_count / state_tot,
           prop_10_plus = age_10_plus/age_tot,
           check = state_tot == age_tot,
           check_age = state_tot < age_tot)
  
  if(any(states$check == 0) & any(states$check_age == 1)){
    stop("state and age totals are unequal: age total is too high")
  } else{
    return(states)
  }
}

dead_prob_E <- 1 - exp(-)
dead_prob_L <- 1 - exp(-() * dt)
dead_prob_P <- 1 - exp(-mup * dt)
