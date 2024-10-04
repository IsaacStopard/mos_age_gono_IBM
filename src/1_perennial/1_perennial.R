orderly2::orderly_shared_resource("default_values.R", 
                                  "IBM_adult_only.R")

source("default_values.R")
source("IBM_adult_only.R")

#####################
##### perennial #####
#####################

p_all <- data.frame("f" = 0, "a" = 0, "o" = o_init)

p_all <- p_all[rep(1:nrow(p_all), n_reps),] |> 
  mutate(sim = rep(1:n_reps, each = nrow(p_all)))

saveRDS(p_all, "p_all.rds")

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("N0", "o_init", "gono_prob", "m_prob", "p_all", "steps", "dt", 
                                   "run_IBM"))

perennial_sims <- foreach(i=1:nrow(p_all),
                     .packages = (.packages())
) %dopar% {
  run_IBM_w_repeat(N0 = N0,
                    o_init = o_init,
                    gono_prob = gono_prob,
                    m_prob = m_prob,
                    birth_model = "freq",
                    birth_params = p_all[i,],
                    steps = steps,
                    dt = dt,
                    max_attempts = max_attempts)
}
stopCluster(cl)

saveRDS(perennial_sims, 
        file = "perennial_sims.rds")


