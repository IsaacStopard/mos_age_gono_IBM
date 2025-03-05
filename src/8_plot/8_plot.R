orderly2::orderly_shared_resource("default_values.R", 
                                  "IBM_adult_only.R")

orderly2::orderly_dependency("4_process_f", "latest", c("freq_sims_p.rds", "auto_sims_p.rds", "perennial_sims_p.rds"))
orderly2::orderly_dependency("5_EHT_autocorr", "latest", c("EHT_autocorr.rds"))
orderly2::orderly_dependency("6_EHT_runs", "latest", c("EHT_sims.rds"))
orderly2::orderly_dependency("7_ssm_inc_larval", "latest", c("ssm_plot_df.rds", "ssm_params.rds"))

##################
##### set up #####
##################

source("default_values.R")
source("IBM_adult_only.R")

freq_sims <- readRDS("freq_sims_p.rds")
auto_sims <- readRDS("auto_sims_p.rds")
perennial_sims <- readRDS("perennial_sims_p.rds")

EHT_autocorr <- readRDS("EHT_autocorr.rds")

EHT_sims <- readRDS("EHT_sims.rds")

o <- m_in * m_prob

theme_set(theme_bw() + 
            theme(text = element_text(size = 14), 
                  legend.background = element_rect(color = NA, fill = NA),
                  legend.text = element_text(size = 14)))

#####################
##### figure 1 #####
####################

################### perennial dynamics

perennial <- perennial_sims |> mutate(season = "perennial")

one_peak <- freq_sims |> filter(round(f, digits = 3) == 0.017 & a == 20) |> mutate(season = "one peak")

limits_df <- rbind(perennial, one_peak) |> group_by(timestep, sim, season) |> dplyr::select(X0:X150) |> pivot_longer(X0:X150) |> 
  mutate(name = as.numeric(gsub("X", "", name))) |> filter(value != 0)

max_count <- max(limits_df$value)
max_age <- max(limits_df$name)

time_plot <- function(df){
  
  mean_abundance <- df |> group_by(timestep, season) |> summarise(m = mean(state_tot))
  
  return(
    ggplot() + 
      geom_line(data = df |> filter(sim != 1), aes(x = timestep, y = state_tot, group = sim, col = "Other simulations"), linewidth = 0.3) +
      geom_line(data = df |> filter(sim == 1), aes(x = timestep, y = state_tot, group = sim, col = "One simulation"), linewidth = 1) +
      geom_line(data = mean_abundance, aes(x = timestep, y = m, col = "Mean"), linewidth = 0.75) +
      ylim(0, 450) +
      theme(legend.position = "inside", legend.position.inside = c(0.1, 0.9)) +
      ylab("Mosquito abundance") + xlab("Day of the year") +
      facet_wrap(~ season) +
      scale_colour_manual(values = c("Other simulations" = "grey30", "One simulation" = "#56B4E9", "Mean" = "black"), name = "")
  )
}

hist_plot <- function(df, sim_rep, time){
  ages <- df |> filter(sim %in% sim_rep & timestep == time)
  
  ages <- ages |> dplyr::select(X0:X150) |> pivot_longer(X0:X150) |> 
    mutate(name = as.numeric(gsub("X", "", name))) |> uncount(value)
  
  tot <- nrow(ages)
  
  return(
    list("tot" = tot,
         "plot" = ggplot(data = ages, aes(x = name)) + 
           geom_histogram(binwidth = 1, fill = "grey70", col = "black") +
           geom_vline(xintercept = mean(ages$name), col = "black", linetype = 2, linewidth = 1) +
           theme_classic() +
           theme(text = element_text(size = 10)) +
           xlab("Age") +
           ylab("Frequency") +
           coord_cartesian(xlim = c(0, 70), ylim = c(0, max_count)) +
           scale_x_continuous(breaks = seq(0, max_age, 20)) +
           annotate("text", x = 50, y = 50, label = paste0("Mean: ", round(mean(ages$name), digits = 1)),
                    col = "black", size = 14/.pt)
    )
  )
}

comb_plot <- function(x1 = 100, 
                      x2 = 300,
                      data = perennial,
                      pos_a1 = c(363.5, 75),
                      pos_a2 = c(1 * 365, 0.647 * 515),
                      pos_b1 = c(0.95, 1.35, 0.01, 0.35),
                      pos_b2 = c(0.8, 1.2, 0.65, 0.99)
){
  
  h_plot_1 <- hist_plot(df = data, sim_rep = 1, time = x1)
  h_plot_2 <- hist_plot(df = data, sim_rep = 1, time = x2)
  
  plot <- list(
    geom_segment(aes(x = pos_a2[[1]], xend = x2,
                     y = pos_a2[[2]], yend = h_plot_2$tot),
                 arrow = arrow(length = unit(0.25, "cm")), col = "#56B4E9", linewidth = 1),
    
    geom_segment(aes(x = pos_a1[[1]], xend = x1, 
                     y = pos_a1[[2]], yend = h_plot_1$tot),
                 arrow = arrow(length = unit(0.25, "cm")), col = "#56B4E9", linewidth = 1),
    
    inset_element(h_plot_2$plot +
                    theme(plot.background = element_rect(colour = "#56B4E9", fill=NA, linewidth = 1.25)),
                  left = pos_b2[[1]], bottom = pos_b2[[3]], right = pos_b2[[2]], top = pos_b2[[4]]),
    
    inset_element(h_plot_1$plot + theme(plot.background = element_rect(colour = "#56B4E9", fill=NA, linewidth = 1.25)),
                  left = pos_b1[[1]], bottom = pos_b1[[3]], right = pos_b1[[2]], top = pos_b1[[4]])
  )
  
  return(plot)
}

hp <- 0.325
wp <- 0.325

exp_plot_p <- time_plot(df = perennial) +
  labs(title = "A") + comb_plot(data = perennial, x1 = 70, x2 = 280,
                                                    pos_b1 = c(0.08, 0.08 + wp, 0.01, 0.01 + hp),
                                                    pos_a1 = c(85, 142.5),
                                                    pos_b2 = c(0.54, 0.54 + wp, 0.61, 0.61 + hp),
                                                    pos_a2 = c(265, 280))

##### age distribution population size 
age_dist_time_plot <- function(df, sim_rep_max, time, time_diff){
  
  ages <- lapply(sim_rep_max, 
         function(m, df){
           ages <- df |> filter(sim %in% seq(1, m) & timestep %in% time)
           ages <- ages |> dplyr::select(X0:X150, timestep) |> pivot_longer(X0:X150) |> 
             mutate(name = as.numeric(gsub("X", "", name)),
                    n_simulations = paste0("Number of pooled simulations: ", m)) |> group_by(timestep, n_simulations) |> uncount(value) |> as.data.frame()
           return(ages)
           },
     df = df
     ) |> bind_rows()
  
  ages$n_simulations <- factor(ages$n_simulations, levels = paste0("Number of pooled simulations: ", sort(sim_rep_max)))
  
  mean_age <- ages |> group_by(n_simulations, timestep) |> summarise(m = mean(name))
  
  plot <- ggplot() + 
    ggridges::geom_density_ridges(data = ages, 
                       aes(x = name, y = factor(timestep), fill = timestep),
      stat = "binline", binwidth = 1, scale = 1.9, 
      draw_baseline = FALSE, alpha = 0.75) +
    coord_cartesian(xlim = c(0, 40)) +
    scale_x_continuous(breaks = seq(0, max_age, 5)) +
    facet_wrap(~n_simulations) +
    theme(legend.position = "none") +
    xlab("Age (days)") + ylab("Day of year") +
    scale_fill_distiller() +
    geom_segment(data = mean_age, 
                 aes(x = m, y = factor(timestep), xend = m, yend = factor(timestep + time_diff)),
                 col = "black", linewidth = 0.7, linetype = 1) +
    scale_y_discrete(breaks = factor(time[round(seq(1, length(time), 3), digits = 0)]))
  
  return(plot)
}

age_dist_plot_p <- age_dist_time_plot(df = perennial, 
                                           sim_rep_max = c(1, 5, 20), 
                                           time = seq(25, 360, 25), time_diff = 25)

ggsave("age_dist_plot_p.pdf",
       age_dist_plot_p,
       device = "pdf",
       #width = 20, height = 25,
       width = 20, height = 12.5,
       units = "cm")

######### seasonal dynamics #########
##### one peak per year ####
exp_plot_s <- time_plot(df = one_peak) + comb_plot(data = one_peak, x1 = 25, x2 = 230,
                                                    pos_b1 = c(0.08, 0.08 + wp, 0.01, 0.01 + hp),
                                                    pos_a1 = c(85, 142.5),
                                                    pos_b2 = c(0.58, 0.58 + wp, 0.6, 0.6 + hp),
                                                    pos_a2 = c(283, 275))

p_plot <- (exp_plot_p + exp_plot_s) + plot_layout(guides = "collect")

season_plot_df <- rbind(perennial, one_peak)

season_plot_df$season <- factor(season_plot_df$season, levels = c("perennial", "one peak"))

sum_age_df <- season_plot_df |> 
  group_by(season, timestep) |> summarise(m = mean(mean),
                                          m_sd = mean(sd_age),
                                          m_cv = mean(cv_age),
                                          m_p = mean(prop_parous),
                                          m_p_10_plus = mean(prop_10_plus),
                                          m_mu = mean(mu_rate))

sum_age_df[sum_age_df$season == "perennial", "m"] |> min()

age_plot <- ggplot() +
  geom_line(data = season_plot_df,
            aes(x = timestep, y = mean, col = season,
                group = interaction(sim, season)), linewidth = 0.15, alpha = 0.25) +
  geom_line(data = sum_age_df, 
                   aes(x = timestep, y = m,
                       col = season), linewidth = 0.9) +
  theme(legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.8)) +
  ylab("Mean mosquito age (days)") +
  xlab("Day of the year") +
  scale_y_continuous(limits = c(4, 11)) +
  scale_colour_manual(values = c("#E69F00", "#CC79A7")) + labs(title = "B")

parity_plot <- ggplot() +
  geom_line(data = season_plot_df,
            aes(x = timestep, y = prop_parous, 
                col = season,
                group = interaction(sim, season)), linewidth = 0.15, alpha = 0.25) +
  geom_line(data = sum_age_df, 
            aes(x = timestep, y = m_p,
                col = season), linewidth = 0.9) +
  theme(legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.8)) +
  ylab("Percent parous") +
  xlab("Day of the year") +
  scale_y_continuous(limits = c(0.25, 0.85), labels = scales::percent) +
  scale_colour_manual(values = c("#E69F00", "#CC79A7"))

age_10_plot <- ggplot() +
  geom_line(data = season_plot_df,
            aes(x = timestep, y = prop_10_plus, 
                col = season,
                group = interaction(sim, season)), linewidth = 0.15, alpha = 0.25) +
  geom_line(data = sum_age_df, 
            aes(x = timestep, y = m_p_10_plus,
                col = season), linewidth = 0.9) +
  theme(legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.8)) +
  ylab("Percent of the population\nat least 10 days old") +
  xlab("Day of the year") +
  scale_y_continuous(limits = c(0.05, 0.55), labels = scales::percent) +
  scale_colour_manual(values = c("#E69F00", "#CC79A7"))

time_metric_plots <- (age_plot + parity_plot + age_10_plot) + 
  plot_layout(guides = "collect",
              axes = "collect")

calc_sum_t_df <- function(.df){
  .df <- .df |> 
  group_by(season, sim) |> summarise(m_age = mean(mean),
                                     sd_m_age = sd(mean),
                                     m_m = mean(state_tot),
                                     sd_m = sd(state_tot),
                                     m_prop_p = mean(prop_parous),
                                     sd_prop_p = sd(prop_parous),
                                     m_prop_10_plus = mean(prop_10_plus),
                                     sd_prop_10_plus = sd(prop_10_plus)) |> 
  mutate(cv_age = sd_m_age / m_age,
         cv_m = sd_m / m_m,
         cv_prop_p = sd_prop_p / m_prop_p,
         cv_prop_10_plus = sd_prop_10_plus / m_prop_10_plus)
  
  .df$season <- factor(.df$season, levels = c("one peak", "perennial"))
  
  return(.df)
}

sum_t_df <- calc_sum_t_df(season_plot_df)

# two week period during transmission season
sum_t_df_p_ts <- calc_sum_t_df(season_plot_df |> filter(round(timestep, digits = 0) >= round(365/4 - 7, digits = 0) &
                                                        round(timestep, digits = 0) <= round(365/4 + 7, digits = 0))
                               )

# transmission season
sum_t_df_ts <- calc_sum_t_df(season_plot_df |> filter(round(timestep, digits = 0) >= round(0, digits = 0) &
                                                          round(timestep, digits = 0) <= round(180, digits = 0))
                             )

cv_plot_fun <- function(.df, l_lim, u_lim, title){
  ggplot(data = .df, aes(y = season)) +
  geom_boxplot(aes(x = cv_m, col = "Mosquito count", fill = "Mosquito count"), alpha = 0.25) +
  geom_boxplot(aes(x = cv_age, col = "Mean age", fill = "Mean age"), alpha = 0.25) +
  geom_boxplot(aes(x = cv_prop_p, col = "Parity", fill = "Parity"), alpha = 0.25) +
  geom_boxplot(aes(x = cv_prop_10_plus, col = "Proportion of mosquitoes\nat least 10 days old", fill = "Proportion of mosquitoes\nat least 10 days old"), alpha = 0.25) +
  scale_x_continuous(limits = c(l_lim, u_lim), breaks = seq(l_lim, u_lim, 0.05), labels = scales::percent) +
  ylab("Seasonality in the\nmean emergence rate") +
  xlab("Coefficient of variation in the daily values") +
  scale_colour_manual(values = c("Proportion of mosquitoes\nat least 10 days old" = "#009E73", 
                                 "Parity" = "#D55E00",
                                 "Mosquito count" = "#0072B2",
                                 "Mean age" = "#000000"), 
                      name = "Entomological metric") +
    scale_fill_manual(values = c("Proportion of mosquitoes\nat least 10 days old" = "#009E73", 
                                   "Parity" = "#D55E00",
                                   "Mosquito count" = "#0072B2",
                                   "Mean age" = "#000000"), 
                        name = "Entomological metric") +
  labs(title = title)
}
  
cv_plot <- cv_plot_fun(sum_t_df_ts, l_lim = 0, u_lim = 0.21, title = "C Primary transmission season (half year)") +
  cv_plot_fun(sum_t_df_p_ts, l_lim = 0, u_lim = 0.21, title = "Primary transmission season (peak two-weeks)") + 
  plot_layout(axes = "collect", guides = "collect")

one_peak_plots <- (p_plot / time_metric_plots / cv_plot) + plot_layout(heights = c(1, 1, 0.875))

ggsave("op_plots.pdf", 
       one_peak_plots,
       device = "pdf",
       width = 50, height = 37.5,
       units = "cm")

freq_dens_plot <- ggplot(data = season_plot_df |> filter(sim %in% seq(1, 5)),
       aes(x = mean, 
           group = interaction(sim, season), 
           fill = factor(season), 
           col = factor(season))) +
  geom_density(alpha = 0.15) +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8)) +
  xlab("Mean daily ages over a single year (days)") +
  ylab("Density") +
  scale_fill_manual(values = c("#E69F00", "#CC79A7"),
                    name = "seasonality") +
  scale_colour_manual(values = c("#E69F00", "#CC79A7"),
                    name = "seasonality") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(2, 14), breaks = seq(2, 14, 2)) + labs(title = "C")


# age_dist_plot_p + 
#   plot_layout(nrow = 1, 
#               heights = c(2.5, 3))

# one_peak_plots <-  + ((age_plot + freq_dens_plot) / (parity_plot + age_10_plot) + 
#                 plot_layout(nrow = 2)) +
#   plot_layout(nrow = 2, heights = c(2, 3))

# values
round(min(perennial$mean), digits = 1)
round(max(perennial$mean), digits = 1)
round(min(one_peak$mean), digits = 1)
round(max(one_peak$mean), digits = 1)

round(min(perennial$mu_rate), digits = 2)
round(max(perennial$mu_rate), digits = 2)
round(min(one_peak$mu_rate), digits = 2)
round(max(one_peak$mu_rate), digits = 2)

round(min(perennial$prop_10_plus), digits = 2) * 100
round(max(perennial$prop_10_plus), digits = 2) * 100
round(min(one_peak$prop_10_plus), digits = 2) * 100
round(max(one_peak$prop_10_plus), digits = 2) * 100

round(min(perennial$prop_parous), digits = 2) * 100
round(max(perennial$prop_parous), digits = 2) * 100
round(min(one_peak$prop_parous), digits = 2) * 100
round(max(one_peak$prop_parous), digits = 2) * 100

round(median(subset(sum_t_df, season == "perennial")$cv_m), digits = 3)*100
round(median(subset(sum_t_df, season == "perennial")$cv_age), digits = 3)*100

round(median(subset(sum_t_df, season == "one peak")$cv_m), digits = 3)*100
round(median(subset(sum_t_df, season == "one peak")$cv_age), digits = 3)*100

round(median(subset(sum_t_df_ts, season == "perennial")$cv_m), digits = 3)*100
round(median(subset(sum_t_df_ts, season == "perennial")$cv_age), digits = 3)*100

round(median(subset(sum_t_df_ts, season == "one peak")$cv_m), digits = 3)*100
round(median(subset(sum_t_df_ts, season == "one peak")$cv_age), digits = 3)*100

round(median(subset(sum_t_df_p_ts, season == "perennial")$cv_m), digits = 3)*100
round(median(subset(sum_t_df_p_ts, season == "perennial")$cv_age), digits = 3)*100

round(median(subset(sum_t_df_p_ts, season == "one peak")$cv_m), digits = 3)*100
round(median(subset(sum_t_df_p_ts, season == "one peak")$cv_age), digits = 3)*100

# mortality
mu_plot <- ggplot() +
  geom_line(data = season_plot_df,
            aes(x = timestep, y = mu_rate, col = season,
                group = interaction(sim, season)), linewidth = 0.1, alpha = 0.25) +
  geom_line(data = sum_age_df, 
            aes(x = timestep, y = m_mu,
                col = season), linewidth = 0.9) +
  theme(legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.175, 0.1)) +
  ylab("Per-capita mosquito mortality rate\nestimated from the daily age distributions (per day)") +
  xlab("Day of the year") +
  scale_y_continuous(limits = c(0, 0.2)) +
  scale_colour_manual(values = c("#E69F00", "#CC79A7"))

ggsave("mu_plot.pdf", 
       mu_plot,
       device = "pdf",
       #width = 20, height = 25, 
       width = 20, height = 15,
       units = "cm")

####################
##### Figure 2 #####
####################

###############################
##### supplementary plots #####
###############################

### sensitivity analysis
freq_sims <- freq_sims |> mutate(w = 365 * f / (2*pi))

# amplitude
a_sens_sims <- subset(freq_sims, round(f, digits = 3) == 0.017) |> 
  mutate(a_p = a/o_init,
         a_p_lab = round(a_p * 100, digits = 1))

m_a_sens_sims <- a_sens_sims |> group_by(f, w, a, timestep, a_p, a_p_lab) |> 
  summarise(m_abundance = mean(state_tot),
            m_0 = mean(X0))

cv_a_sens_sims <- a_sens_sims |> group_by(f, w, a, a_p, a_p_lab, sim) |> 
  summarise(m_age = mean(mean),
            sd_age = sd(mean)) |> mutate(cv = sd_age / m_age)

supp_a_plot <- ggplot(data = a_sens_sims, 
                      aes(x = timestep)) +
  geom_line(aes(y = state_tot, group = sim, col = "total"), linewidth = 0.1, alpha = 0.5) +
  geom_line(aes(y = X0, group = sim, col = "aged 0"), linewidth = 0.1, alpha = 0.5) +
  geom_line(data = m_a_sens_sims, 
            aes(y = m_abundance, col = "total"), linewidth = 0.75) +
  geom_line(data = m_a_sens_sims, 
            aes(y = m_0, col = "aged 0"), linewidth = 0.75) +
  facet_wrap(~a_p_lab) +
  theme_bw() +
  ylab("Mosquito abundance") + xlab("Day of year") +
  scale_colour_manual(values = c("aged 0" = "skyblue", "total" = "grey50"),
                      name = "") + 
  guides(colour = guide_legend(override.aes = list(linewidth = 2)))

a_sens_plot <- ggplot(data = cv_a_sens_sims, 
                      aes(x = a_p, y = cv, group = factor(a_p))) +
  geom_boxplot() +
  ylim(0, 0.25) +
  xlab(expression(paste("Amplitude (",eta,") as a % of the mean emergence rate (o)"))) +
  ylab("Coefficient of variation in the\ndaily mean ages over a single year") +
  scale_x_continuous(labels = scales::percent) +
  theme(legend.position = "none")

# frequency
f_sens_sims <- subset(freq_sims, a == 10) |> 
  mutate(a_p = a/o_init)

m_f_sens_sims <- f_sens_sims |> group_by(w, f, a, timestep) |> 
  summarise(m_abundance = mean(state_tot),
            m_0 = mean(X0))

# mean values over the
sum_freq_f <- f_sens_sims |> group_by(a, f, sim, w) |>
  summarise(m = mean(mean),
            sd = sd(mean),
            cv = sd/m)

supp_f_plot <- ggplot(data = f_sens_sims, 
                      aes(x = timestep)) +
  geom_line(aes(y = state_tot, group = sim, col = "total"), linewidth = 0.1, alpha = 0.5) +
  geom_line(aes(y = X0, group = sim, col = "aged 0"), linewidth = 0.1, alpha = 0.5) +
  geom_line(data = m_f_sens_sims, aes(y = m_abundance, col = "total"), linewidth = 0.75) +
  geom_line(data = m_f_sens_sims, aes(y = m_0, col = "aged 0"), linewidth = 0.75) +
  facet_wrap(~w) +
  theme_bw() +
  ylab("Mosquito abundance") + xlab("Day of year") +
  scale_colour_manual(values = c("aged 0" = "skyblue", "total" = "grey50"),
                      name = "") + 
  guides(colour = guide_legend(override.aes = list(linewidth = 2)))

freq_sens_plot <- ggplot(data = sum_freq_f,  
                    aes(x = w, y = cv, group = factor(w))) +
  geom_boxplot() +
  xlab(expression(paste("Frequency of seasonality in ", lambda,  " (365", kappa, ") (rad per year)"))) +
  ylab("Coefficient of variation in the\ndaily mean ages over a single year") +
  theme_bw() +
  scale_x_continuous(breaks = c(seq(0, 50, 10))) +
  theme(text = element_text(size = 14)) +
  scale_y_continuous(limits = c(0, 0.25))

ggsave("supp_a_plot.pdf", 
       supp_a_plot,
       device = "pdf",
       width = 20, height = 15, 
       units = "cm")

ggsave("supp_f_plot.pdf", 
       supp_f_plot,
       device = "pdf",
       width = 20, height = 15, 
       units = "cm")

sensitivity_f_plot <- freq_sens_plot + a_sens_plot +
  plot_layout(nrow = 1) + plot_annotation(tag_levels = 'A') +
  plot_layout(axis = "collect")

ggsave("sensitivity_f_plot.pdf", 
       sensitivity_f_plot,
       device = "pdf",
       #width = 20, height = 25, 
       width = 30, height = 15,
       units = "cm")

####################
##### Figure 3 #####
####################

####### autocorrelation #######

auto_sims_plot <- auto_sims |> subset(o_sd == 5)

auto_sum <- auto_sims |> group_by(h, o_sd, sim) |> summarise(m = mean(mean),
                                                             sd = sd(mean),
                                                             cv = sd/m)

er_plot <- ggplot(data = auto_sims_plot |> subset(sim == 1 & h %in% c(-0.9, 0, 0.9)) |> mutate(h = paste0("h: ", h)), 
       aes(x = timestep+1, y = tot_lag, col = factor(h))) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = c("#0072B2", "black", "#D55E00"), name = "") + 
  facet_wrap(~h) +
  theme(legend.position = "none") +
  ylab(expression(paste("Mean mosquito emergence rate (",lambda,") (per day)"))) + xlab("Day of year")

abundance_plot <- ggplot(data = auto_sims_plot |> subset(sim == 1 & h %in% c(-0.9, 0, 0.9)) |> mutate(h = paste0("h: ", h)),
         aes(x = timestep)) +
  geom_line(aes(y = X0, col = "Mosquitoes aged 0"), linewidth = 0.5) +
  geom_line(aes(y = state_tot, col = "All mosquitoes"), linewidth = 0.5) +
  scale_colour_manual(values = c("All mosquitoes" = "black", "Mosquitoes aged 0" = "skyblue"), name = "") + 
  facet_wrap(~h) +
  theme_bw() + theme(text = element_text(size = 14), legend.position = "inside", legend.position.inside = c(0.1, 0.5)) +
  ylab("Mosquito abundance") + 
  xlab("Day of year") + 
  guides(colour = guide_legend(override.aes = list(linewidth = 2)))

auto_dens_plot <- ggplot(data = auto_sims_plot |> filter(h %in% c(-0.9, 0, 0.9)) |> mutate(h = paste0("h: ", h)),
                      aes(x = mean, group = interaction(sim, h), fill = factor(h))) +
  geom_density(alpha = 0.1, col = "black", linewidth = 0.1) +
  theme(legend.position = "inside", legend.position.inside = c(0.95, 0.875), 
        legend.background = element_rect(color = NA, fill = NA)) +
  xlab("Mean daily ages over a single year (days)") +
  ylab("Density") +
  scale_fill_manual(values = c("#0072B2", "grey", "#D55E00"), name = "") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(2, 14), breaks = seq(2, 14, 2)) +
  facet_wrap(~h) + 
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

auto_sens_plot <- ggplot(data = auto_sum, 
                         aes(x = h, y = cv, 
                             group = interaction(h, o_sd), col = factor(o_sd))) +
  geom_boxplot() +
  ylim(0, 0.25) +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8), text = element_text(size = 14)) +
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#56B4E9"), name = expression(paste(sigma))) +
  scale_fill_manual(values = c("#CC79A7", "#E69F00", "#56B4E9"), name = expression(paste(sigma))) +
  ylab("Coefficient of variation in the\nmean daily ages over a single year") +
  xlab("Autocorrelation (h)") +
  scale_x_continuous(breaks = seq(-0.9, 0.9, 0.2))

ggsave("auto_sens_plot.pdf", 
       auto_sens_plot,
       device = "pdf",
       #width = 20, height = 25, 
       width = 20, height = 15,
       units = "cm")
  
autocorr_plots <- er_plot / abundance_plot / auto_dens_plot + plot_annotation(tag_levels = 'A')

ggsave("autocorr_plots.pdf", 
       autocorr_plots,
       device = "pdf",
       #width = 20, height = 25, 
       width = 30, height = 35,
       units = "cm")

supp_auto_plots <- ggplot(data = auto_sims |> subset(sim == 1 & h %in% seq(-0.9, 0.9, 0.2)),
                          aes(x = timestep)) +
  geom_line(aes(y = X0, col = "Mosquitoes aged 0"), linewidth = 0.5) +
  geom_line(aes(y = state_tot, col = "All mosquitoes"), linewidth = 0.5) +
  scale_colour_manual(values = c("All mosquitoes" = "black", "Mosquitoes aged 0" = "skyblue"), name = "") + 
  facet_grid(vars(h), vars(o_sd)) +
  theme_bw() + theme(text = element_text(size = 14)) +
  ylab("Mosquito abundance") + 
  xlab("Day of year") + 
  guides(colour = guide_legend(override.aes = list(linewidth = 2)))

ggsave("supp_auto_plots.pdf", 
       supp_auto_plots,
       device = "pdf",
       #width = 20, height = 25, 
       width = 30, height = 50,
       units = "cm")

# mu_plot <- ggplot() +
#   geom_line(data = rbind(perennial, one_peak),
#             aes(x = timestep, y = mu_rate, 
#                 col = season,
#                 group = interaction(sim, season)), linewidth = 0.05, alpha = 0.25) +
#   geom_line(data = sum_age_df, 
#             aes(x = timestep, y = m_mu,
#                 col = season), linewidth = 0.9) +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         legend.position = "inside",
#         legend.position.inside = c(0.1, 0.9),
#         text = element_text(size = 14)) +
#   ylab("Estimated mortality rate\nfrom the age distribution") +
#   xlab("Day of the year") +
#   #scale_y_continuous(limits = c(0.4, 0.75), labels = scales::percent) +
#   scale_colour_manual(values = c("#D55E00", "#0072B2"))
# 
# ggsave("exp_plot.pdf", 
#        exp_plot,
#        device = "pdf",
#        #width = 20, height = 25, 
#        width = 35, height = 35,
#        units = "cm")

##########################
##### EHT model fits #####
##########################

EHT_fit_plot <- function(out){
  ggplot(data = out$df) +
    geom_ribbon(data = out$p$pred_m_df, aes(x = day_continuous, ymin = lower, ymax = upper), fill = "grey50", alpha = 0.225) +
    geom_line(aes(x = day_continuous, y = total, col = treat), linewidth = 0.4, alpha = 0.5) +
    geom_point(aes(x = day_continuous, y = total, shape = factor(hut), col = treat, fill = treat), size = 2.75, alpha = 0.7) +
    geom_line(data = out$p$pred_m_df, aes(x = day_continuous, y = median), linewidth = 1.5) +
    theme_bw() + ylab("Total mosquito count\nper experimental hut") + xlab("Days since start of EHT") +
    scale_shape_manual(values = c(21, 22, 23, 24, 25, 7), name = "Experimental\nhut") +
    theme(text = element_text(size = 14), 
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 8),
          legend.box = "horizontal") +
    scale_colour_viridis_d(name = "ITN") +
    scale_fill_viridis_d(name = "ITN") +
    scale_y_continuous(limits = c(0, 225), breaks = seq(0, 200, 50))
}

T_plot <-  EHT_fit_plot(out = EHT_autocorr$out_T) +
  ggtitle("Tengrela")

V_plot <- EHT_fit_plot(out = EHT_autocorr$out_V) +
  ggtitle("Vallée du Kou 5")

median(c(EHT_autocorr$out_T$df$total, EHT_autocorr$out_V$df$total)) * 6 * 3

dens_plot <- function(par, 
                      f_T = EHT_autocorr$out_T, 
                      f_V = EHT_autocorr$out_V){
  df <- rbind(rstan::extract(f_T$fit, par) %>% as.data.frame() %>% mutate(study = "Tengrela"),
              rstan::extract(f_V$fit, par) %>% as.data.frame() %>% mutate(study = "Vallée du Kou 5"))
  
  colnames(df)[1] <- "par"
  
  ggplot(data = df) +
    geom_density(aes(par, col = study, fill = study), alpha = 0.5) +
    theme_bw()  +
    theme(text = element_text(size = 14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylab("Posterior density") +
    scale_fill_manual(values = c("#E69F00", "#CC79A7"), name = "") +
    scale_colour_manual(values = c("#E69F00", "#CC79A7"), name = "")
}

h_plot <- dens_plot(par = "h") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
  xlab("Estimated lag-1 autocorrelation (h)") +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 1) +
  theme(legend.position = "inside", legend.position.inside = c(0.1, 0.9))

sd_plot <- dens_plot(par = "sigma") +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, 1)) +
  xlab(expression(paste("Lag-1 regression model standard deviation (", sigma, ")")))

c_plot <- dens_plot("c") +
  scale_x_continuous(limits = c(0, 12.5), breaks = seq(0, 12.5, 2.5)) +
  xlab("Lag-1 regression model mean (c)")

kappa_plot <- dens_plot("kappa") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 2.5, 0.5)) +
  xlab(expression(paste("Overdispersion parameter (", phi,")")))

ggsave("EHT_par_plots.pdf",
       c_plot + sd_plot + kappa_plot + 
         plot_layout(nrow = 1, guides = "collect") +
         plot_annotation(tag_levels = c("A")),
       width = 25,
       height = 5.75
)
  
####################
##### EHT runs #####
####################

list2env(EHT_sims, envir = .GlobalEnv)

process_EHT_sims <- function(EHT_sims,
                             n_warmup){
  EHT_sims <- lapply(1:length(EHT_sims), function(i){
  EHT_sims[[i]] |> 
    filter(timestep > n_warmup & timestep < nrow(EHT_sims[[i]])) |> 
    mutate(timestep = row_number(),
           sim = i)
  }) |> bind_rows()
  
  EHT_sims_mean <- EHT_sims |> group_by(timestep) |> summarise(m = mean(state_tot))
  
  return(list("sims" = EHT_sims, "mean" = EHT_sims_mean))
}

EHT_sims_T <- process_EHT_sims(EHT_sims_T, n_warmup = n_warmup)

EHT_sims_V <- process_EHT_sims(EHT_sims_V, n_warmup = n_warmup)

EHT_sim_plot <- ggplot() +
  geom_line(data = rbind(EHT_sims_T$sims |> mutate(EHT = "Tengrela"),
                         EHT_sims_V$sims |> mutate(EHT = "Vallée du Kou 5")), 
            aes(x = timestep, y = state_tot, group = interaction(factor(sim), EHT),
                col = EHT),
            linewidth = 0.075, alpha = 0.1) +
  geom_line(data = rbind(EHT_sims_T$mean |> mutate(EHT = "Tengrela"),
                         EHT_sims_V$mean |> mutate(EHT = "Vallée du Kou 5")),
            aes(x = timestep, y = m, col = EHT), linewidth = 1) +
  ylim(0, 500) +
  xlab("Day of EHT") + ylab("Simulated mosquito abundance") +
  scale_colour_manual(values = c("#E69F00", "#CC79A7"), name = "") +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.9))


ages_T <- EHT_sims_T$sims |> filter(sim %in% seq(1)) |> group_by(timestep) |> 
  dplyr::select(X0:X150) |> tidyr::pivot_longer(X0:X150) |> 
  dplyr::mutate(name = as.numeric(gsub("X", "", name))) |> tidyr::uncount(value) |> 
  mutate(EHT = "Tengrela")

ages_V <- EHT_sims_V$sims |> filter(sim %in% seq(1)) |> group_by(timestep) |> 
  dplyr::select(X0:X150) |> tidyr::pivot_longer(X0:X150) |> 
  dplyr::mutate(name = as.numeric(gsub("X", "", name))) |> tidyr::uncount(value) |> 
  mutate(EHT = "Vallée du Kou 5")

EHT_ages <- rbind(ages_T, ages_V)

mean_age_EHT <- EHT_ages |> group_by(EHT, timestep) |> summarise(m = mean(name))

EHT_ages_plot <- ggplot() +
  ggridges::geom_density_ridges(data = EHT_ages, 
                                aes(x = name, y = factor(timestep), fill = timestep),
                                stat = "binline", binwidth = 1, scale = 1.9, 
                                draw_baseline = FALSE, alpha = 0.75) +
  geom_segment(data = mean_age_EHT, 
               aes(x = m, y = factor(timestep), xend = m, yend = factor(timestep + 1)),
               col = "black", linewidth = 0.7, linetype = 1) +
  coord_cartesian(xlim = c(0, 50)) +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_y_discrete(breaks = seq(39, 1, -1)) +
  theme(legend.position = "none") +
  xlab("Age (days)") + ylab("Day of EHT") +
  scale_fill_distiller() +
  facet_wrap(~EHT) +
  ggtitle("Simulated mosquito age distributions")

ggsave("EHT_plot.pdf",
       (T_plot + V_plot + plot_layout(guides = "collect"))/
         (h_plot + EHT_sim_plot) /
         EHT_ages_plot +
         plot_layout(nrow = 3, heights = c(1, 1, 2)) +
         plot_annotation(tag_levels = c("A")),
       width = 17.5,
       height = 21.5
)

#####################
##### ITN model #####
#####################

ssm_plot_df <- readRDS("ssm_plot_df.rds")
ssm_params <- readRDS("ssm_params.rds")

get_vals_ssm <- function(time_in = 210, ITN_time_in, seasonality_in){
  return(left_join(
    subset(ssm_plot_df, ITN_time == ITN_time_in & 
             seasonality == seasonality_in &
                   t_plot == time_in & 
             itn_cov == 1) |> select(-itn_cov),
    
          subset(ssm_plot_df, ITN_time == ITN_time_in & seasonality == seasonality_in &
                   t_plot == time_in & itn_cov == 0) |> select(-itn_cov) |> 
            rename(value_0 = value, mean_age_0 = mean_age)
          ) |> mutate(diff_m = value_0 - value,
                      diff_mp = diff_m / value_0 * 100,
                      diff_age = mean_age_0 - mean_age,
                      diff_agep = diff_age / mean_age_0 * 100) |>
    group_by(t_plot, seasonality, ITN_time) |> 
  summarise(m_diff_m = round(mean(diff_m), digits = 1),
            m_diff_mp = round(mean(diff_mp), digits = 1),
            m_diff_age = round(mean(diff_age), digits = 1),
            m_diff_agep = round(mean(diff_agep), digits = 1),
            l_diff_m = round(min(diff_m), digits = 1),
            l_diff_mp = round(min(diff_mp), digits = 1),
            l_diff_age = round(min(diff_age), digits = 1),
            l_diff_agep = round(min(diff_agep), digits = 1),
            u_diff_m = round(max(diff_m), digits = 1),
            u_diff_mp = round(max(diff_mp), digits = 1),
            u_diff_age = round(max(diff_age), digits = 1),
            u_diff_agep = round(max(diff_agep), digits = 1)))
}

get_vals_ssm(time_in = 210, ITN_time_in = "ITN time: population increasing", seasonality_in = "perennial")
get_vals_ssm(time_in = 210, ITN_time_in = "ITN time: population increasing", seasonality_in = "seasonal")

get_vals_ssm(time_in = 310, ITN_time_in = "ITN time: population decreasing", seasonality_in = "perennial")
get_vals_ssm(time_in = 310, ITN_time_in = "ITN time: population decreasing", seasonality_in = "seasonal")

get_age_ssm <- function(time_in = 210, ITN_time_in, seasonality_in, itn_cov_in){
  return(subset(ssm_plot_df, ITN_time == ITN_time_in & 
         seasonality == seasonality_in &
         t_plot == time_in & 
         itn_cov == itn_cov_in) |> group_by(ITN_time, t_plot, seasonality, itn_cov) |> 
  summarise(m_mean_age = round(mean(mean_age), digits = 1),
            l_mean_age = round(min(mean_age), digits = 1),
            u_mean_age = round(max(mean_age), digits = 1)))
}

get_age_ssm(time_in = 210, ITN_time_in = "ITN time: population increasing", seasonality_in = "perennial", itn_cov_in = 0)
get_age_ssm(time_in = 210, ITN_time_in = "ITN time: population increasing", seasonality_in = "perennial", itn_cov_in = 1)

get_age_ssm(time_in = 210, ITN_time_in = "ITN time: population increasing", seasonality_in = "seasonal", itn_cov_in = 0)
get_age_ssm(time_in = 210, ITN_time_in = "ITN time: population increasing", seasonality_in = "seasonal", itn_cov_in = 1)

get_age_ssm(time_in = 310, ITN_time_in = "ITN time: population decreasing", seasonality_in = "seasonal", itn_cov_in = 0)
get_age_ssm(time_in = 310, ITN_time_in = "ITN time: population decreasing", seasonality_in = "seasonal", itn_cov_in = 1)

mos_plot <- ggplot() +
  geom_line(data = ssm_plot_df,
            aes(x = t_plot, y = value, group = interaction(seasonality, itn_cov, ITN_time, name), col = factor(itn_cov)),
            linewidth = 0.05, alpha = 0.1) +
  geom_line(data = ssm_plot_df |> group_by(seasonality, itn_cov, ITN_time, t_plot) |> summarise(m = mean(value)),
            aes(x = t_plot, y = m, group = interaction(seasonality, itn_cov, ITN_time), col = factor(itn_cov))) +
  geom_vline(data = ssm_params |> mutate(ITN_IRS_on = ITN_IRS_on - 4 * 365), 
             aes(xintercept = ITN_IRS_on),
             col = "black", linetype = 2) +
  facet_grid(vars(seasonality), vars(factor(ITN_time, levels = c("ITN time: population increasing", "ITN time: population decreasing")))) +
  theme_bw() +
  scale_colour_manual(values = c("#E69F00", "#56B4E9"), name = "ITN coverage") +
  xlab("Day") +
  ylab("Total mosquito abundance") 

age_plot <- 
  ggplot() + 
  geom_line(data = ssm_plot_df |> na.omit(),
            aes(x = t_plot, y = mean_age, 
                group = interaction(seasonality, itn_cov, ITN_time, name), 
                col = factor(itn_cov)),
            linewidth = 0.05, alpha = 0.1) +
  geom_line(data = ssm_plot_df |> na.omit() |> 
              group_by(seasonality, itn_cov, ITN_time, t_plot) |> summarise(m = mean(mean_age)),
            aes(x = t_plot, y = m, group = interaction(seasonality, itn_cov, ITN_time), col = factor(itn_cov))) +
  geom_vline(data = ssm_params |> mutate(ITN_IRS_on = ITN_IRS_on - 4 * 365), 
             aes(xintercept = ITN_IRS_on),
             col = "black", linetype = 2) +
  facet_grid(vars(seasonality),
             vars(factor(ITN_time, levels = c("ITN time: population increasing", "ITN time: population decreasing")))) +
  theme_bw() +
  scale_colour_manual(values = c("#E69F00", "#56B4E9"), name = "ITN coverage") +
  xlab("Day") +
  ylab("Mean mosquito age (days)") +
  coord_cartesian(ylim = c(0, 30))

mean_age_diff <- ssm_plot_df |> filter(itn_cov == 0) |> 
  left_join(ssm_plot_df |> filter(itn_cov == 1), 
            by = c("name", "t", "ITN_time", "t_plot", "seasonality")) |> 
  mutate(diff_age = mean_age.x - mean_age.y)

diff_plot <- ggplot() +
  geom_point(data = mean_age_diff |> na.omit(), 
             aes(x = t_plot, y = diff_age, col = ITN_time, group = interaction(name)),
             size = 0.25, alpha = 0.1) +
  geom_line(data = mean_age_diff |> na.omit() |> group_by(seasonality, t_plot, ITN_time) |> summarise(m = mean(diff_age)),
            aes(x = t_plot, y = m, col = ITN_time),
            linewidth = 1) +
  
  facet_wrap(~seasonality) +
  coord_cartesian(ylim = c(-5, 30))+
  scale_colour_manual(values = c("skyblue", "#009E73"), name = "") +
  xlab("Day") +
  ylab("Increase in mean mosquito age without ITNs (days)") +
  theme_bw() +
  geom_vline(data = ssm_params |> mutate(ITN_IRS_on = ITN_IRS_on - 4 * 365), 
             aes(xintercept = ITN_IRS_on, col = ITN_time),
             linetype = 2)

int_plots <- mos_plot + age_plot + diff_plot + 
  plot_layout(nrow = 3) + plot_annotation(tag_levels = 'A')

ggsave("int_plots.pdf", 
       int_plots,
       device = "pdf",
       #width = 20, height = 25, 
       width = 30, height = 40,
       units = "cm")
