# default values and packages

library(individual); library(tidyverse); library(foreach); library(doParallel); library(patchwork)
library(survival); library(flexsurv)

# setting up the initial compartments
max_t <- 365 * 5
dt <- 1 # must be 1 for mosquito aging process
steps <- max_t/dt # only run with daily time steps

# parameters
N0 <- m_in <- 200
mu <- 0.132 # constant mortality rate
m_prob <- 1 - exp(-mu * dt) # probability of adult mortality
o_init <- N0 * mu #N0 * m_prob

gono_rate <- 1/3 # biting rate
gono_prob <- 1 - exp(-gono_rate * dt)

# number of simulations
n_reps <- 20

max_attempts <- 5

