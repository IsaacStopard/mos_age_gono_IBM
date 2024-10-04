orderly2::orderly_shared_resource("default_values.R", 
                                  "IBM_adult_only.R")

orderly2::orderly_dependency("5_EHT_autocorr", "latest", c("EHT_autocorr.rds"))

source("default_values.R")
source("IBM_adult_only.R")

EHT_fits <- readRDS("EHT_autocorr.rds")

library(rstan); library(ggridges)

# extracting the emergence rates
n_iter <- rstan::extract(EHT_fits$out_T$fit, "m0")[[1]] |> nrow()

n_runs <- 750
n_warmup <- 50

set.seed(12345)
inds <- sample(1:n_iter, n_runs, replace = FALSE)

saveRDS(inds, file = "inds.rds")

run_model_MCMC_iter <- function(fit, 
                                i,
                                n_warmup){
  
  m0_d <- rstan::extract(fit, "m0")[[1]][i]
  m0 <- m0_d |> round(digits = 0)
  
  N_EHT <- m0 * 6 # because there are multiple experimental huts
  
  sc <- N_EHT/m0
  
  o_t <- rstan::extract(fit, "o")[[1]][i,]*sc
  
  c_e <- m_prob * N_EHT
  
  out <- run_IBM_w_repeat(N0 = N_EHT,
                          o_init = c_e,
                          gono_prob = gono_prob,
                          m_prob = m_prob,
                          birth_model = "user",
                          birth_params = c(rep(c_e, n_warmup), o_t, c_e),
                          steps = n_warmup + length(o_t),
                          dt = dt,
                          max_attempts = max_attempts)
  return(out)
}

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("n_warmup", "inds", "EHT_fits", 
                                   "gono_prob", "m_prob", "dt", 
                                   "run_IBM", "run_IBM_w_repeat",
                                   "run_model_MCMC_iter"))

EHT_sims_T <- foreach(i = 1:length(inds),
                      .packages = (.packages())
                      ) %dopar% {
                        run_model_MCMC_iter(fit = EHT_fits$out_T$fit,
                                            i = inds[i],
                                            n_warmup = n_warmup)
                      }

EHT_sims_V <- foreach(i = 1:length(inds),
                      .packages = (.packages())
) %dopar% {
  run_model_MCMC_iter(fit = EHT_fits$out_V$fit,
                      i = inds[i],
                      n_warmup = n_warmup)
}

stopCluster(cl)

saveRDS(list("EHT_sims_T" = EHT_sims_T, 
             "EHT_sims_V" = EHT_sims_V,
             "n_warmup" = n_warmup,
             "n_iter" = n_iter),
        file = "EHT_sims.rds")
