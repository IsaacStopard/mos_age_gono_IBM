# State space model

## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

#################################
##### larval dynamics model #####
#################################

ts <- (time / 365) 
tf <- ts - floor(ts)

pi <- user(3.141593) 

# the parameters for the fourier series
g0 <- user()
g1 <- user()
g2 <- user()
g3 <- user()
h1 <- user()
h2 <- user()
h3 <- user()

rainfall_floor <- user(0.0001)

rainfall <- if(g0 == 0 && g1  == 0 && g2  == 0 && g3  == 0 && h1  == 0 && h2  == 0 && h3  == 0) 1 else max(g0 +
                                                                                                             g1 * cos(2 * pi * tf * 1) + g2 * cos(2 * pi * tf * 2) + g3 * cos(2 * pi * tf * 3) +
                                                                                                             h1 * sin(2 * pi * tf * 1) + h2 * sin(2 * pi * tf * 2) + h3 * sin(2 * pi * tf * 3), rainfall_floor)

r_f_pred[1:365] <- if(g0 == 0 && g1  == 0 && g2  == 0 && g3  == 0 && h1  == 0 && h2  == 0 && h3  == 0) 1 else max(g0 +
                                                                                                                    g1 * cos(2 * pi * i/365 * 1) + g2 * cos(2 * pi * i/365 * 2) + g3 * cos(2 * pi * i/365 * 3) +
                                                                                                                    h1 * sin(2 * pi * i/365 * 1) + h2 * sin(2 * pi * i/365 * 2) + h3 * sin(2 * pi * i/365 * 3), rainfall_floor)

dim(r_f_pred) <- 365

me <- user(0.0338)
ml <- user(0.0348)
gamma <- user(13.25)
mup <- user(0.249)

# transition probabilities
del <- user(6.64)
dl <- user(3.72)
dpl <- user(0.643)

mum <- user(0.132) # adult mosquito mortality rate
beta <- user(21.2)
blood_meal_rates <- user(1/3)


eov <- beta / mum * (exp(mum / blood_meal_rates) - 1)
beta_larval <- eov * mum * exp(-mum / blood_meal_rates) / (1 - exp(-mum / blood_meal_rates))

sub_omega <- gamma * ml / me - (del / dl) + ((gamma - 1) * ml * del)
omega <- -0.5 * sub_omega + sqrt(0.25 * sub_omega^2 + 0.5 * gamma * beta * ml * del / (me * mum * dl * (1. + dpl * mup)))

m <- user()

K0 <- m * 2 * dl * mum * (1. + dpl * mup) * gamma * (omega + 1) / (omega / (ml * del) - (1. / (ml * dl)) - 1) 
r_bar <- sum(r_f_pred)/365 # mean rainfall throughout the year
K <- K0 * rainfall / r_bar

# calculates the probability of an event happening from a constant rate

# probability early larval stages are alive
dead_prob_E <- 1 - exp(-(me * (1 + (E + L) / K)) * dt)
dead_prob_L <- 1 - exp(-(ml * gamma * (1 + (E + L) / K)) * dt)
dead_prob_P <- 1 - exp(-mup * dt)

develop_prob_E <- 1 - exp((-1/del) * dt)
develop_prob_L <- 1 - exp((-1/dl) * dt)
develop_prob_P <- 1 - exp((-1/dpl) * dt)

initial(E) <- init_E
init_E <- user(0)

initial(L) <- init_L
init_L <- user(0)

initial(P) <- init_P
init_P <- user(0)

n_births <- rpois(beta_larval * N)
# keep track of the number of births as another variable.

# # (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
# deriv(E) <- beta_larval * N - me*(1+(E+L)/K)*E - E/del
# # egg hatching - den. dep. mortality - maturing larvae
# deriv(L) <- E/del - ml*(1+gamma*(E + L)/K)*L - L/dl
# # pupae - mortality - fully developed pupae
# deriv(P) <- L/dl - mup*P - P/dpl

alive_E <- rbinom(E, (1 - dead_prob_E))
alive_undevelop_E <- rbinom(alive_E, (1 - develop_prob_E))
alive_develop_E <- alive_E - alive_undevelop_E

alive_L <- rbinom(L, (1 - dead_prob_L))
alive_undevelop_L <- rbinom(alive_L, (1 - develop_prob_L))
alive_develop_L <- alive_L - alive_undevelop_L

alive_P <- rbinom(P, (1 - dead_prob_P))
alive_undevelop_P <- rbinom(alive_P, (1 - develop_prob_P))
alive_develop_P <- rbinom(alive_P - alive_undevelop_P, 0.5) # divided by two to give just the female mosquitoes

# d_P_f <- rbinom(E, d_prob_P/2)
# d_P_m <- rbinom(E, d_prob_P/2)
# 
update(E) <- n_births + alive_undevelop_E
update(L) <- alive_develop_E + alive_undevelop_L
update(P) <- alive_develop_L + alive_undevelop_P

#######################################################
##### state space adult mosquito population model #####
#######################################################

n_age <- user(200)
init_A[1] <- m
init_A[2:n_age] <- 0
dim(init_A) <- n_age

## Model parameters (default in parenthesis)
mu[1:n_age] <- mum # constant mortality rate
dim(mu) <- n_age

## Initial conditions
initial(A[]) <- init_A[i]
dim(A) <- n_age

# Mortality process
m_prob[1:n_age] <- 1 - exp(-mu[i] * dt) # Individual probabilities of transition
dim(m_prob) <- n_age

n_mu_A[1:n_age] <- rbinom(A[i], m_prob[i]) # Draws from binomial distributions for numbers changing between compartments
dim(n_mu_A) <- n_age

update(A[1]) <- alive_develop_P
update(A[2:n_age]) <- (A[(i-1)] - n_mu_A[(i-1)])

initial(N) <- sum(init_A)
update(N) <- sum(A)
