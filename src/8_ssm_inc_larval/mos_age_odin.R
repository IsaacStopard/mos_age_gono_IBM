# State space model

## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

######################
##### parameters #####
######################

m0 <- user()

# the parameters for the fourier series
g0 <- user()
g1 <- user()
g2 <- user()
g3 <- user()
h1 <- user()
h2 <- user()
h3 <- user()

### larval model

# mortality rate
me <- user(0.0338)
ml <- user(0.0348)
mup <- user(0.249)

# development times
del <- user(6.64)
dl <- user(3.72)
dpl <- user(0.643)

gamma <- user(13.25)
beta <- user(21.2)

### adult mosquito model

mu0 <- user(0.132) # adult mosquito mortality rate

pi <- user(3.141593) 
rainfall_floor <- user(0.001)

# setting up variable sizes
n_age <- user(100) 

tau1 <- user(0.69) # duration of host-seeking behaviour
tau2 <- user(2.31) # duration of resting behaviour
p10 <- exp(-mu0 * tau1) # prob of surviving 1 feeding cycle
p2 <- exp(-mu0 * tau2) #prob of surviving one resting cycle

###########################
##### initial numbers #####
###########################

initial(E) <- init_E
init_E <- user(0)

initial(L) <- init_L
init_L <- user(0)

initial(P) <- init_P
init_P <- user(0)

#######################
##### seasonality #####
#######################

### rainfall seasonality 
ts <- (time / 365) 
tf <- ts - floor(ts)

rainfall <- if(g0 == 0 && g1  == 0 && g2  == 0 && g3  == 0 && h1  == 0 && h2  == 0 && h3  == 0) 1 else max(g0 +
  g1 * cos(2 * pi * tf * 1) + g2 * cos(2 * pi * tf * 2) + g3 * cos(2 * pi * tf * 3) +
  h1 * sin(2 * pi * tf * 1) + h2 * sin(2 * pi * tf * 2) + h3 * sin(2 * pi * tf * 3), rainfall_floor)

r_f_pred[1:365] <- if(g0 == 0 && g1  == 0 && g2  == 0 && g3  == 0 && h1  == 0 && h2  == 0 && h3  == 0) 1 else max(g0 +
                         g1 * cos(2 * pi * i/365 * 1) + g2 * cos(2 * pi * i/365 * 2) + g3 * cos(2 * pi * i/365 * 3) +
                         h1 * sin(2 * pi * i/365 * 1) + h2 * sin(2 * pi * i/365 * 2) + h3 * sin(2 * pi * i/365 * 3), rainfall_floor)

dim(r_f_pred) <- 365

r_bar <- sum(r_f_pred)/365 # mean rainfall throughout the year

### seasonality in the carrying capacity
eov <- beta / mu * (exp(mu / fv) - 1)
beta_larval <- eov * mu * exp(-mu / fv) / (1 - exp(-mu / fv))

sub_omega <- gamma * ml / me - (del / dl) + ((gamma - 1) * ml * del)
omega <- -0.5 * sub_omega + sqrt(0.25 * sub_omega^2 + 0.5 * gamma * beta_larval * ml * del / (me * mu0 * dl * (1. + dpl * mup)))
K0 <- m0 * 2 * dl * mu0 * (1. + dpl * mup) * gamma * (omega + 1) / (omega / (ml * del) - (1. / (ml * dl)) - 1) 
K <- K0 * rainfall / r_bar

#################################
##### larval dynamics model #####
#################################

# # (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
# deriv(E) <- beta_larval * N - me*(1+(E+L)/K)*E - E/del
# # egg hatching - den. dep. mortality - maturing larvae
# deriv(L) <- E/del - ml*(1+gamma*(E + L)/K)*L - L/dl
# # pupae - mortality - fully developed pupae
# deriv(P) <- L/dl - mup*P - P/dpl

# calculates the probability of an event happening from a constant rate
# probability early larval stages are alive

dead_rate_E <- me * (1 + (E + L) / K)
develop_rate_E <- 1 / del
tot_rate_E <- dead_rate_E + develop_rate_E
move_E <- rbinom(E, 1 - exp(-tot_rate_E * dt))
dead_E <- rbinom(move_E, dead_rate_E / tot_rate_E)
develop_E <- move_E - dead_E

dead_rate_L <- ml * gamma * (1 + (E + L) / K)
develop_rate_L <- 1 / dl
tot_rate_L <- dead_rate_L + develop_rate_L

move_L <- rbinom(L, 1 - exp(-tot_rate_L * dt))
dead_L <- rbinom(move_L, dead_rate_L / tot_rate_L)
develop_L <- move_L - dead_L

dead_rate_P <- mup
develop_rate_P <- 1 / dpl
tot_rate_P <- dead_rate_P + develop_rate_P
move_P <- rbinom(P, 1 - exp(-tot_rate_P * dt))
dead_P <- rbinom(move_P, dead_rate_P / tot_rate_P)
develop_P <- move_P - dead_P

n_births <- rpois(beta_larval * N)

update(E) <- E + n_births - dead_E - develop_E
update(L) <- L + develop_E - dead_L - develop_L
update(P) <- P + develop_L - develop_P - dead_P

develop_P_to_A <- rbinom(develop_P, 0.5)

#######################################################
##### state space adult mosquito population model #####
#######################################################

N <- sum(A)

initial(Nout) <- m0
update(Nout) <- sum(A)

init_A[1] <- m0
init_A[2:n_age] <- 0
dim(init_A) <- n_age

## Model parameters (default in parenthesis)
mu_A[1:n_age] <- mu # constant mortality rate
dim(mu_A) <- n_age

## Initial conditions
initial(A[]) <- init_A[i]
dim(A) <- n_age

# Mortality process
m_prob[1:n_age] <- 1 - exp(-mu_A[i] * dt) # Individual probabilities of transition
dim(m_prob) <- n_age

n_mu_A[1:n_age] <- rbinom(A[i], m_prob[i]) # Draws from binomial distributions for numbers changing between compartments
dim(n_mu_A) <- n_age

update(A[1]) <- develop_P_to_A
update(A[2:n_age]) <- (A[(i-1)] - n_mu_A[(i-1)])


##############################
##### intervention model #####
##############################

# intervention parameters

irs_loss <- user()
itn_loss <- user()

ITN_IRS_on <- user() # days after which interventions begin
num_int <- user() # number of intervention categorys, ITN only, IRS only, neither, both
irs_cov <- user(0) # proportion of population covered by IRS

# cov is a vector of coverages for each intervention category:
#dim(itn_vector) <- user()
#itn_vector[] <- user()
#dim(t_vector) <- length(itn_vector)
#t_vector[] <- user()
eff_itn_cov <- user() #interpolate(t_vector, itn_vector, "constant")
#int_itn_irs_on <- interpolate(t_vector, t_vector, "constant")
#eff_ITN_IRS_on <- if (time < ITN_IRS_on) ITN_IRS_on else time + ITN_IRS_on #int_itn_irs_on

dim(cov_) <- 4
cov_[1] <- (1-eff_itn_cov)*(1-irs_cov)  # {No intervention}
cov_[2] <- eff_itn_cov*(1-irs_cov) # 	   {ITN only}
cov_[3] <- (1-eff_itn_cov)*irs_cov	#      {IRS only}
cov_[4] <- eff_itn_cov*irs_cov #	   {Both ITN and IRS}
cov[] <- cov_[i]
dim(cov) <- num_int

IRS_interval <- user(1 * 365) # how long IRS lasts
ITN_interval <- user(3 * 365) # how long ITN lasts

chi <- user(0.86) # proportion of vector endophily
Q0 <- user(0.92) # proportion of anthropophagy
bites_Bed <- user(0.89) # endophagy in bed
bites_Indoors <- user(0.97) # endophagy indoors

# General intervention model terminology:
# r - probability of trying to repeat feed after hitting ITN/IRS
# d - probability of dying after hitting ITN/IRS
# s - probability of successful feed after hitting ITN/IRS

# The maximum (and then minimum) r and d values for ITN/IRS on day 0 before they decay
r_ITN0 <- user(0.56)
d_ITN0 <- user(0.41)
d_IRS0 <- user(1)
r_IRS0 <- user(0.6)
r_ITN1 <- user(0.24)

# Calculates decay for ITN/IRS
ITN_decay = if(time < ITN_IRS_on) 0 else exp(-((time - ITN_IRS_on)%%ITN_interval) * itn_loss)
IRS_decay = if(time < ITN_IRS_on) 0 else exp(-((time - ITN_IRS_on)%%IRS_interval) * irs_loss)

# The r,d and s values turn on after ITN_IRS_on and decay accordingly
d_ITN <- if(time < ITN_IRS_on) 0 else d_ITN0*ITN_decay
r_ITN <- if(time < ITN_IRS_on) 0 else r_ITN1 + (r_ITN0 - r_ITN1)*ITN_decay
s_ITN <- if(time < ITN_IRS_on) 1 else 1 - d_ITN - r_ITN

r_IRS <- if(time < ITN_IRS_on) 0 else r_IRS0*IRS_decay
d_IRS <- if(time < ITN_IRS_on) 0 else chi*d_IRS0*IRS_decay
s_IRS <- if(time < ITN_IRS_on) 1 else 1 - d_IRS

# probability that mosquito bites and survives for each intervention category
dim(w_) <- 4
w_[1] <- 1
w_[2] <- 1 - bites_Bed + bites_Bed*s_ITN
w_[3] <- 1 - bites_Indoors + bites_Indoors*(1-r_IRS)*s_IRS
w_[4] <- 1 - bites_Indoors + bites_Bed*(1-r_IRS)*s_ITN*s_IRS + (bites_Indoors - bites_Bed)*(1-r_IRS)*s_IRS
w[] <- w_[i]
dim(w) <- num_int

# probability that mosquito is repelled during a single attempt for each int. cat.
dim(z_) <- 4
z_[1] <- 0
z_[2] <- bites_Bed*r_ITN
z_[3] <- bites_Indoors*r_IRS
z_[4] <- bites_Bed*(r_IRS+ (1-r_IRS)*r_ITN) + (bites_Indoors - bites_Bed)*r_IRS
z[] <- z_[i]
dim(z) <- num_int

# Calculating Z (zbar) and W (wbar) - see Supplementary materials 2 for details
dim(zhi) <- num_int
dim(whi) <- num_int
zhi[1:num_int] <- cov[i]*z[i]
whi[1:num_int] <- cov[i]*w[i]
zh <- if(time < ITN_IRS_on) 0 else sum(zhi)
wh <- if(time < ITN_IRS_on) 1 else sum(whi)

# Z (zbar) - average probability of mosquito trying again during single feeding attempt
zbar <- Q0*zh
# W (wbar) - average probability of mosquito successfully feeding during single attempt
wbar <- 1 - Q0 + Q0*wh

# p1 is the updated p10 given that interventions are now in place:
p1 <- wbar*p10/(1-zbar*p10)

# adult mosquito mortality and feeding rates after ITN use is accounted for
initial(fv) <- 1/( tau1/(1-zbar) + tau2 )
update(fv) <- 1/( tau1/(1-zbar) + tau2 ) # mosquito feeding rate (zbar from intervention param.)

initial(mu) <- -fv*log(p1*p2) 
update(mu) <- -fv*log(p1*p2) # mosquito death rate

