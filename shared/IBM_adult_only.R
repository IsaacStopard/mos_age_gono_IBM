# function to run the IBM
# m_prob - daily mortality probability
run_IBM <- function(N0,
                    o_init,
                    gono_prob,
                    m_prob,
                    birth_model,
                    birth_params,
                    steps,
                    dt){
  
  # putting the parameters in the environment
  # states
  # NP - nulliparous
  # P - paraous
  # TL - total lagged by one time point
  mos_states <- c("NP", "P")
  mos_states_t0 <- rep("NP", N0) # creating the initial states of each mosquito
  
  # initialising the numbers in each 
  mos_pop <- CategoricalVariable$new(categories = mos_states, 
                                     initial_values = mos_states_t0)
  
  # age variable to keep track of the ages
  # all mosquitoes start at age zero
  age <- DoubleVariable$new(initial_values = rep(0, N0)) 
  
  tot_lag <- DoubleVariable$new(initial_values = o_init)
  tot_deaths <- DoubleVariable$new(initial_values = 0)
  tot_births <- DoubleVariable$new(initial_values = 0)
  
  # event to determine age when parous
  parity_event <- TargetedEvent$new(population_size = N0)
  
  parity_event$add_listener(
    function(t, target){mos_pop$queue_update("P", target)}
    )
  
  ### processes
  parity_process <- function(t){
    NP <- mos_pop$get_index_of("NP")
    already_scheduled <- parity_event$get_scheduled()
    NP$and(already_scheduled$not(inplace = TRUE))
    gono_times <- rgeom(n = NP$size(), prob = gono_prob) + 1
    parity_event$schedule(target = NP, delay = gono_times)
  }
  
  # converting mortality rate to the daily survival probability
  death_process <- function(t){
    pop_size <- mos_pop$size()
    
    deaths <- sample.int(n = pop_size, size = rbinom(1, pop_size, m_prob), replace = FALSE)
    
    n_deaths <- length(deaths)
    tot_deaths$queue_update(values = n_deaths, index = 1)
    
    # updating the values
    mos_pop$queue_shrink(deaths)
    age$queue_shrink(deaths)
    parity_event$queue_shrink(deaths)
  }
  
  age_process <- function(t){
    mos_ages <- age$get_values() + dt
    age$queue_update(mos_ages, index = NULL)
  }
  
 if(birth_model == "auto"){
    
    # params
    # o_init - initial number emerged
    # o_bar - mean number emerging
    # h - autocorrelation
    # o_sd - standard deviation
    o_bar <- birth_params[["o_bar"]]
    h <- birth_params[["h"]]
    o_sd <- birth_params[["o_sd"]]
    
    birth_process <- function(t){
      
      o_lag <- tot_lag$get_values(index = 1)
      
      mu_births <- rnorm(1,
                         mean = o_bar + h * (o_lag - o_bar),
                         sd = o_sd * sqrt(1 - h^2)
                         )
      
      tot_lag$queue_update(values = mu_births, index = 1)
      
      n_births <- rpois(1, lambda = mu_births)
      
      tot_births$queue_update(values = n_births, index = 1)
      
      if(n_births > 0){
        mos_pop$queue_extend(rep('NP', n_births))
        age$queue_extend(rep(0, n_births))
        parity_event$queue_extend(n_births)
      }
    }
    
  } else if(birth_model == "freq"){
    # params
    f <- birth_params[["f"]]
    a <- birth_params[["a"]]
    o <- birth_params[["o"]]
    
    birth_process <- function(t){
      
      m_b <- o * dt + a * dt * sin(t * f)
      
      n_births <- rpois(1, lambda = m_b)
      
      tot_births$queue_update(values = n_births, index = 1)
      
      if(n_births > 0){
        mos_pop$queue_extend(rep('NP', n_births))
        age$queue_extend(rep(0, n_births))
        parity_event$queue_extend(n_births)
      }
      
      tot_lag$queue_update(values = m_b, index = 1)
      
    }
    
  } else if(birth_model == "user"){
    
    all_i <- seq(0, steps, dt)
    
    birth_process <- function(t){
      i <- which(all_i == t)
      n_births <- rpois(1, lambda = birth_params[[i]])
      tot_births$queue_update(values = n_births, index = 1)
      
      if(n_births > 0){
        mos_pop$queue_extend(rep('NP', n_births))
        age$queue_extend(rep(0, n_births))
        parity_event$queue_extend(n_births)
      }
      
      tot_lag$queue_update(values = birth_params[[i]], index = 1)
      
    }
    
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
  
  lag_render <- double_renderer(
    renderer = mos_render,
    variable = tot_lag,
    category = "tot_lag"
  )
  
  death_render <- double_renderer(
    renderer = mos_render,
    variable = tot_deaths,
    category = "tot_deaths"
  )
  
  birth_render <- double_renderer(
    renderer = mos_render,
    variable = tot_births,
    category = "tot_births"
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

run_IBM_w_repeat <- function(N0,
                             o_init,
                             gono_prob,
                             m_prob,
                             birth_model,
                             birth_params,
                             steps,
                             dt,
                             max_attempts){
  attempts <- 1
  states <- tryCatch({run_IBM(N0 = N0,
                    o_init = o_init,
                    gono_prob = gono_prob,
                    m_prob = m_prob,
                    birth_model = birth_model,
                    birth_params = birth_params,
                    steps = steps,
                    dt = dt)},
                    error = function(cond){
                      return(NULL)
                    })

  while(attempts <= max_attempts & is.null(states)){
    attempts <- attempts + 1
    states <- tryCatch({run_IBM(N0 = N0,
                      o_init = o_init,
                      gono_prob = gono_prob,
                      m_prob = m_prob,
                      birth_model = birth_model,
                      birth_params = birth_params,
                      steps = steps,
                      dt = dt)},
                      error = function(cond){
                        return(NULL)
                      })
  }
  
  if(!is.null(states)){
    states <- states |> mutate(attempts = attempts)
  }
  
  return(states )
}
