library(odin.dust); library(mcstate); library(tidyverse); library(malariasimulation); library(coda); library(foresite)

orderly2::orderly_shared_resource("default_values.R", 
                                  "IBM_adult_only.R")

source("default_values.R")

# compile model
age_model_dust <- odin.dust::odin_dust("mos_age_odin.R")

theme_set(theme_bw() + 
            theme(text = element_text(size = 14), 
                  legend.background = element_rect(color = NA, fill = NA),
                  legend.text = element_text(size = 14)))

#####################
##### functions #####
#####################

run_model <- function(num_int, 
                      itn_cov, 
                      irs_cov, 
                      ITN_IRS_on, 
                      n_p = 50, 
                      timesteps = steps,
                      season){
  
  # Catch all: Not defined the correct number of interventions
  if (itn_cov > 0 & num_int == 1){
    stop(message("Incorrect number of interventions for definied ITN coverage. Please ensure you have correctly
                 specified the number of interventions."))
  }
  if (irs_cov > 0 & num_int < 3){
    stop(message("Incorrect number of interventions for definied IRS coverage. Please ensure you have correctly
                 specified the number of interventions."))
  }
  
  # Sets start time of coverage
  
  age_model <- age_model_dust$new(pars = list(num_int = num_int,
                                              m0 = 200,
                                              g0 = season[["g0"]],
                                              g1 = season[["g1"]],
                                              g2 = season[["g2"]],
                                              g3 = season[["g3"]],
                                              h1 = season[["h1"]],
                                              h2 = season[["h2"]],
                                              h3 = season[["h3"]],
                                              eff_itn_cov = itn_cov,
                                              irs_cov = irs_cov,
                                              num_int = num_int,
                                              ITN_IRS_on = ITN_IRS_on,
                                              irs_loss = log(2)/(0.5 * 365), # need to check
                                              itn_loss = log(2)/(2.64 * 365),
                                              rainfall_floor = 0.1),
                                  time = 1,
                                  n_particles = n_p,
                                  n_threads = 4L,
                                  seed = 1L)
  
  
  x <- array(NA, dim = c(age_model$info()$len, n_p, timesteps))
  
  for (t in seq_len(timesteps)) {
    x[ , , t] <- age_model$run(t)
  }
  
  return(x)
}

# m refers to the total counts of mosquitoes
extract_m <- function(x, i){
  if(i == 5){
    times <- x[1, 1, ] - 1
  } else{
    times <- x[1, 1, ]
  }
  t(x[i,,]) %>% as.data.frame() %>% mutate(t = times) %>% 
    pivot_longer(-t)
}

extract_age <- function(x, t){
  
  times <- x[1, 1, ]
  
  out <- x[8:(8 + 99),,t] 
  
  # checking the total numbers match
  if(t < (steps)){
    if(any(colSums(out) != x[5, , (t+1)])){
    warning("total mosquito numbers do not match")
    }
    }
  
  return(out |> as.data.frame() |> mutate(age = row_number() - 1) |> 
    pivot_longer(-age) |> 
    mutate(t = times[t]))
}

# to discard the initial years of the simulations
s_times <- function(.df){
  .df |> filter(t > 4 * 365) |> 
    mutate(t_plot = t - 4 * 365)
}

calc_mean_age <- function(x){
  all_ages <- lapply(1:(steps),
                     extract_age,
                     x = x) |> bind_rows()
  
  out <- unique(all_ages[,c("name", "t")])
  
  all_age_long <- all_ages |> 
    tidyr::uncount(value)
  
  all_age_mean <- all_age_long |> group_by(name, t) |> summarise(mean_age = mean(age))
  
  out <- left_join(out, all_age_mean, by = c("name", "t"))
  
  return(out)
}

##################################
##### seasonality parameters #####
##################################

season_bf <- foresite::BFA$seasonality |> filter(name_1 == "Cascades") |> select(g0, g1, g2, g3, h1, h2, h3) |> unlist()
season_perennial <- c("g0" = 0, "g1" = 0, "g2" = 0, "g3" = 0, "h1" = 0, "h2" = 0, "h3" = 0)

######################
##### model runs #####
######################

params <- expand.grid("itn_cov" = c(0, 1), 
            "ITN_IRS_on" = c((180 + 4 * 365), (280 + 4 * 365)),
            "season" = c("season_bf", "season_perennial")
            ) |> mutate(ITN_time = ifelse(ITN_IRS_on == (180 + 4 * 365), "ITN time: population increasing", "ITN time: population decreasing"),
                        seasonality = ifelse(season == "season_bf", "seasonal", "perennial"))

results <- vector(mode = "list", length = nrow(params))

for(i in 1:nrow(params)){
  results[[i]] <- run_model(num_int = 2, itn_cov = params[i, "itn_cov"], 
                            irs_cov = 0, ITN_IRS_on = params[i, "ITN_IRS_on"], 
                            season = get(paste0(params[i, "season"])))
}

# checks
for(i in 1:length(results)){
  if(any(dim(results[[1]]) != c(107, 50, 1825))){
    warning("dims not correct")
  }
}

# processing
total_mos <- vector(mode = "list", length = nrow(params))

for(i in 1:nrow(params)){
  total_mos[[i]] <- extract_m(x = results[[i]], i = 5) |> 
    mutate(ITN_time = params[i, "ITN_time"],
           seasonality = params[i, "seasonality"],
           itn_cov = params[i, "itn_cov"])
}

for(i in 1:length(results)){
  if(nrow(total_mos[[i]]) != 50 * 1825){
    warning("dims not correct")
  }
}

total_mos <- total_mos |> bind_rows() |> s_times()

mean_ages <- vector(mode = "list", length = nrow(params))

for(i in 1:nrow(params)){
  mean_ages[[i]] <- calc_mean_age(x = results[[i]]) |> 
    mutate(ITN_time = params[i, "ITN_time"],
           seasonality = params[i, "seasonality"],
           itn_cov = params[i, "itn_cov"])
}

mean_ages <- mean_ages |> bind_rows() |> s_times()

plot_df <- left_join(total_mos, mean_ages)

if(any(plot_df[is.na(plot_df$mean_age), "value"] |> unique() != 0)){
  warning("missing mean age values")
}

saveRDS(plot_df, file = "ssm_plot_df.rds")
saveRDS(params, file = "ssm_params.rds")
