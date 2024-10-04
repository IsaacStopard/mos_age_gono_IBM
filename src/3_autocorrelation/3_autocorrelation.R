orderly2::orderly_shared_resource("default_values.R", 
                                  "IBM_adult_only.R")

source("default_values.R")
source("IBM_adult_only.R")

###########################
##### autocorrelation #####
###########################

# yearly fluctuations
h <- seq(-0.9, 0.9, 0.1) # select f values that give an index that is an integer
o_bar <- o_init
o_sd <- c(1, 3, 5)

a_all <- expand.grid("h" = h, "o_bar" = o_bar, "o_sd" = o_sd)

a_all <- a_all[rep(1:nrow(a_all), n_reps),] |> mutate(sim = rep(1:n_reps, each = nrow(a_all)))

saveRDS(a_all, "a_all.rds")

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("N0", "o_init", "gono_prob", "m_prob", "a_all", "steps", "dt", 
                                   "run_IBM"))

auto_sims <- foreach(i=1:nrow(a_all),
                     .packages = (.packages())
) %dopar% {
  run_IBM_w_repeat(N0 = N0,
                     o_init = o_init,
                     gono_prob = gono_prob,
                     m_prob = m_prob,
                     birth_model = "auto",
                     birth_params = a_all[i,],
                     steps = steps,
                     dt = dt,
                     max_attempts = max_attempts)
}

stopCluster(cl)

saveRDS(auto_sims, 
        file = "auto_sims.rds")

