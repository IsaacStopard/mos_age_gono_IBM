orderly2::orderly_shared_resource("default_values.R", 
                                  "IBM_adult_only.R")

source("default_values.R")
source("IBM_adult_only.R")

#####################
##### frequency #####
#####################

### differences in amplitude
f <- c(2*pi/365) # yearly fluctuations

a <- seq(1, floor(o_init * 0.825), 1) |> round(digits = 0)

f_all_a <- expand.grid("f" = f, "a" = a, "o" = o_init)

f_all_a <- f_all_a[rep(1:nrow(f_all_a), n_reps),] |> mutate(sim = rep(1:n_reps, each = nrow(f_all_a)))

### differences in frequency
f_f <- (seq(2, 50, 4) * 2 * pi)/365

a_f <- c(min(a), a[floor(length(a)/2)])#floor(o)/2

f_all_f <- expand.grid("f" = f_f, "a" = a_f, "o" = o_init)

f_all_f <- f_all_f[rep(1:nrow(f_all_f), n_reps),] |> mutate(sim = rep(1:n_reps, each = nrow(f_all_f)))

f_all <- rbind(f_all_a, f_all_f)

saveRDS(f_all, "f_all.rds")

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("N0", "o_init", "gono_prob", "m_prob", "f_all", "steps", "dt", 
                                   "run_IBM"))

freq_sims <- foreach(i=1:nrow(f_all),
                       .packages = (.packages())
) %dopar% {
  run_IBM_w_repeat(N0 = N0,
                    o_init = o_init,
                    gono_prob = gono_prob,
                    m_prob = m_prob,
                    birth_model = "freq",
                    birth_params = f_all[i,],
                    steps = steps,
                    dt = dt,
                    max_attempts = max_attempts)
}

stopCluster(cl)

saveRDS(freq_sims, 
        file = "freq_sims.rds")
