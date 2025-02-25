orderly2::orderly_shared_resource("default_values.R", 
                                  "IBM_adult_only.R")

orderly2::orderly_dependency("1_perennial", "latest", c("p_all.rds", "perennial_sims.rds"))
orderly2::orderly_dependency("2_frequency", "latest", c("f_all.rds", "freq_sims.rds"))
orderly2::orderly_dependency("3_autocorrelation", "latest", c("a_all.rds", "auto_sims.rds"))

source("default_values.R")
source("IBM_adult_only.R")

p_all <- readRDS("p_all.rds")
perennial_sims <- readRDS("perennial_sims.rds")

f_all <- readRDS("f_all.rds")
freq_sims <- readRDS("freq_sims.rds")

a_all <- readRDS("a_all.rds")
auto_sims <- readRDS("auto_sims.rds")

##################################
##### processing the results #####
##################################

process_sim <- function(i, all_sims, info, sim_type){
  out <- if(sim_type == "freq"){
    all_sims[[i]] |> 
      mutate(sim = info[i, "sim"],
           f = info[i, "f"],
           a = info[i, "a"])
  } else{
    all_sims[[i]] |> 
      mutate(sim = info[i, "sim"],
             h = info[i, "h"],
             o_sd = info[i, "o_sd"])
  }
  
  out <- out |> filter(timestep <= steps & timestep > steps - (365/dt)) |> 
    mutate(timestep = timestep - (steps - 365))
  
  if(any((out$state_tot > out$age_tot) == 1)){
    return(NULL)
  } else{
    return(out)
    }
  
}

perennial_sims <- lapply(1:length(perennial_sims), 
                      process_sim, 
                      all_sims = perennial_sims, 
                      info = p_all, 
                      sim_type = "freq") |> bind_rows()

freq_sims <- lapply(1:length(freq_sims), 
                    process_sim, 
                    all_sims = freq_sims, 
                    info = f_all, 
                    sim_type = "freq") |> bind_rows()

auto_sims <- lapply(1:length(auto_sims), process_sim, 
                    all_sims = auto_sims, info = a_all, sim_type = "auto") |> bind_rows()

#######################################################
##### quantifying changes in the age distribution #####
#######################################################

### estimating the mortality rate

estimate_mu <- function(i, 
                        df){
  
  ages <- df[i,]
  
  if(nrow(ages) > 1){
    stop("too many rows selected")
  }
  
  ages <- ages |> dplyr::select(X0:X150) |> tidyr::pivot_longer(X0:X150) |> 
    dplyr::mutate(name = as.numeric(gsub("X", "", name))) |> tidyr::uncount(value) |> 
    mutate(status = 1) |> 
    mutate(s_value = name + 1e-5) |> 
    as.data.frame() 
  
  exp_model <- survival::survreg(survival::Surv(time = s_value, event = status) ~ 1, 
                                 data = ages, dist = "exponential")
  
  mu_rate <- 1 / exp(coef(exp_model))
  
  mean <- mean(ages[,"name"])
  sd <- sd(ages[,"name"])
  cv <- sd/mean
  lower <- quantile(ages[,"name"], probs =  0.25)[[1]]
  upper <- quantile(ages[,"name"], probs =  0.75)[[1]]
  median <- median(ages[,"name"])
  
  return(data.frame("mean" = mean,
                    "sd_age" = sd,
                    "cv_age" = cv,
                    "mu_rate" = mu_rate[[1]], 
                    "median" = median, 
                    "lower" = lower, 
                    "upper" = upper)
         )
}

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("perennial_sims", "estimate_mu",
                                   "freq_sims", "auto_sims"))
p_mu <- foreach(i=1:nrow(perennial_sims),
                .packages = (.packages())) %dopar% {
                  tryCatch({estimate_mu(
                    i = i,
                    df = perennial_sims)
                  }, error = function(cond){
                    return(NA)
                  }
                  )
                }

a_mu <- foreach(i=1:nrow(auto_sims),
                .packages = (.packages())) %dopar% {
                  tryCatch({estimate_mu(
                    i = i,
                    df = auto_sims)
                  }, error = function(cond){
                    return(NA)
                  }
                  )
                }

f_mu <- foreach(i=1:nrow(freq_sims),
                .packages = (.packages())) %dopar% {
                  tryCatch({estimate_mu(
                    i = i,
                    df = freq_sims)
                  }, error = function(cond){
                    return(NA)
                  }
                  )
                }
stopCluster(cl)

perennial_sims <- perennial_sims |> cbind(bind_rows(p_mu))
freq_sims <- freq_sims |> cbind(bind_rows(f_mu))
auto_sims <- auto_sims |> cbind(bind_rows(a_mu))

saveRDS(perennial_sims, file = "perennial_sims_p.rds")
saveRDS(freq_sims, file = "freq_sims_p.rds")
saveRDS(auto_sims, file = "auto_sims_p.rds")
