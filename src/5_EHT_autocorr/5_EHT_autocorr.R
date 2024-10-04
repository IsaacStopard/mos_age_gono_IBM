library(readxl); library(rstan); library(pammtools); 
library(haven); library(lubridate); library(MASS); library(RColorBrewer)

orderly2::orderly_shared_resource("default_values.R")

source("default_values.R")

options(mc.cores = 4)
rstan_options(auto_write = TRUE)

#####################
##### functions #####
#####################

# fit each treatment separately
data_in_fun <- function(df){
  return(list(max_t = as.integer(max(df$day_continuous)),
              unq_t = as.integer(seq(1, max(df$day_continuous))),
              N = as.integer(nrow(df)),
              t = as.integer(df$day_continuous),
              M = as.integer(df$total),
              mu = mu))
}

extract_posteriors <- function(fit, df){
  
  pred_m <- rstan::extract(fit, "pred_m")$pred_m
  
  pred_m_df <- cbind(data.frame(t = seq(1, max(df$day_continuous))),
                     apply(pred_m, 2, quantile, probs = c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame())
  
  colnames(pred_m_df) <- c("day_continuous", "lower", "median", "upper")
  
  o <- rstan::extract(fit, "o")$o
  
  o_df <- cbind(data.frame(t = seq(1, max(df$day_continuous)-1)),
                apply(o, 2, quantile, probs = c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame())
  
  colnames(o_df) <- c("day_continuous", "lower", "median", "upper")
  
  return(list(pred_m_df = pred_m_df,
              o_df = o_df))
  
}

run_model <- function(df, iterations = 7500, control_in = list(adapt_delta = 0.99, stepsize = 1)){
  
  data_in <- data_in_fun(df)
  
  warmup <- round(iterations/2, digits = 0)
  chains <- 4
  
  #fit <- stan(file = "ode_model.stan", data = data_in, iter=iterations, chains = chains, seed=12345,
  #   warmup = warmup)
  
  fit_nc <- stan(file = "ode_model.stan", 
                 data = data_in, iter=iterations, chains = chains, seed=12345,
                 warmup = warmup, control = control_in)
  
  p <- extract_posteriors(fit = fit_nc, df = df)
  
  return(list(df = df,
              data_in = data_in,
              fit = fit_nc,
              p = p))
}

# Toe et al

T_df <- read_excel("EHT_data/Toe_et_al_2018.xlsx") %>% mutate(Week = ceiling(day/6),
                                                                   day_week = day - (Week - 1) * 6,
                                                                   day_continuous = (Week - 1) + day)

set.seed(123)
out_T <- run_model(df = subset(T_df, locality == "Tengrela")[,c("day_continuous", "hut", "treat", "total")])
out_V <- run_model(df = subset(T_df, locality == "VK5")[,c("day_continuous", "hut", "treat", "total")])

out_T$df <- out_T$df %>% mutate(hut = as.numeric(gsub("TEN-H", "", hut)))
out_V$df <- out_V$df %>% mutate(hut = as.numeric(gsub("VK5-H", "", hut)))

saveRDS(list("out_T" = out_T, "out_V" = out_V), file = "EHT_autocorr.rds")
