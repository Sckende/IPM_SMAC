# IPMs scripts 
# Bayesian population analysis using WinBUGS (book)
# Chapter 11 - Example
# Real data example: Hoopoe population dynamics
# 3 types of demographic data
# (1) a population count = max number of simultaneously occupied nest boxes in each year
# (2) capture-recapture data of adults and young
# (3) data on reproductive output
# OBJECTIVES : Estimate the demographic parameters and their temporal variability to understand the dynamics of the hoopoe populations
# Annual variation of pop resulting of survival, productivity, immigration & emigration
# short-lived species, so 2 stage-class are considered
# pop surveyed just before the repro = prebreeding census (1 year old ind. are differentiated from older ones)
# Stochasticity is included with appropriate distribution
rm(list = ls())
setwd("C:/Users/Etudiant/Desktop/SMAC/Projet_publi/3-Modelisation_translocation_petrels/SCRIPTS/Bayesian_pop_analysis_WinBUGS_real_data_example")

# Load packages
library(R2WinBUGS)
library(coda)
library(lme4)

# Localization of WinBUGS.exe
bugs.dir <- "C:/Users/Etudiant/Documents/WinBUGS14/"

# Load data
nyears <- 9

# Capture-recapture data: m-array of juveniles and adults (these are males and females together)
marray.j <- matrix (c(15, 3, 0, 0, 0, 0, 0, 0, 198, 0, 34, 9, 1, 0, 0, 0, 0, 287, 0, 0, 56, 8, 1, 0, 0, 0, 455, 0, 0, 0, 48, 3, 1, 0, 0, 518, 0, 0, 0, 0, 45, 13, 2, 0, 463, 0, 0, 0, 0, 0, 27, 7, 0, 493, 0, 0, 0, 0, 0, 0, 37, 3, 434, 0, 0, 0, 0, 0, 0, 0, 39, 405),
                    nrow = 8,
                    ncol = 9,
                    byrow = TRUE)

marray.a <- matrix(c(14, 2, 0, 0, 0, 0, 0, 0, 43, 0, 22, 4, 0, 0, 0, 0, 0, 44, 0, 0, 34, 2, 0, 0, 0, 0, 79, 0, 0, 0, 51, 3, 0, 0, 0, 94, 0, 0, 0, 0, 45, 3, 0, 0, 118, 0, 0, 0, 0, 0, 44, 3, 0, 113, 0, 0, 0, 0, 0, 0, 48, 2, 99, 0, 0, 0, 0, 0, 0, 0, 51, 90),
                   nrow = 8,
                   ncol = 9,
                   byrow = TRUE)

# Population count data
popcount <- c(32, 42, 64, 85, 82, 78, 73, 69, 79)

# Productivity data
J <- c(189, 274, 398, 538, 520, 476, 463, 438) # Number of offspring
R <- c(28, 36, 57, 77, 81, 83, 77, 72) # Number of surveyed broods

# Specify model in BUGS language
sink("ipm.hoopoe.bug")
cat("
model {
# -------------------------------------------------------
# Integrated population model
# - Age structured model with 2 age classes:
# 1-year old and adults (at least 2-years old)
# - Age at first breeding = 1 year
# - Prebreeding census, female-based
# - All vital rates are assumed to be time-dependent (random)
# - Explicit estimation of immigration
# -------------------------------------------------------

# ---------------------------------------
# 1. Define the priors for the parameters
# ---------------------------------------
# Initial population sizes
N1[1] ~ dnorm(100, 0.0001)I(0,) # 1-year old individuals
NadSurv[1] ~ dnorm(100, 0.0001)I(0,) # Adults >= 2 years
Nadimm[1] ~ dnorm(100, 0.0001)I(0,) # Immigrants

# Mean demographic parameters (on appropriate scale)
l.mphij ~ dnorm(0, 0.0001)I(-10,10) # Bounded to help with convergence
l.mphia ~ dnorm(0, 0.0001)I(-10,10)
l.mfec ~ dnorm(0, 0.0001)I(-10,10)
l.mim ~ dnorm(0, 0.0001)I(-10,10)
l.p ~ dnorm(0, 0.0001)I(-10,10)

# Precision of standard deviations of temporal variability
sig.phij ~ dunif(0, 10)
tau.phij <- pow(sig.phij, -2)
sig.phia ~ dunif(0, 10)
tau.phia <- pow(sig.phia, -2)
sig.fec ~ dunif(0, 10)
tau.fec <- pow(sig.fec, -2)
sig.im ~ dunif(0, 10)
tau.im <- pow(sig.im, -2)

# Distribution of error terms (bounded to help with convergence)
for (t in 1:(nyears-1)){
epsilon.phij[t] ~ dnorm(0, tau.phij)I(-15,15)
epsilon.phia[t] ~ dnorm(0, tau.phia)I(-15,15)
epsilon.fec[t] ~ dnorm(0, tau.fec)I(-15,15)
epsilon.im[t] ~ dnorm(0, tau.im)I(-15,15)
}

# -----------------------
# 2. Constrain parameters
# -----------------------
for (t in 1:(nyears-1)){
logit(phij[t]) <- l.mphij + epsilon.phij[t] # Juv. apparent survival
logit(phia[t]) <- l.mphia + epsilon.phia[t] # Adult apparent survival
log(f[t]) <- l.mfec + epsilon.fec[t] # Productivity
log(omega[t]) <- l.mim + epsilon.im[t] # Immigration
logit(p[t]) <- l.p # Recapture probability
}

# -----------------------
# 3. Derived parameters
# -----------------------
mphij <- exp(l.mphij)/(1+exp(l.mphij)) # Mean juvenile survival probability
mphia <- exp(l.mphia)/(1+exp(l.mphia)) # Mean adult survival probability
mfec <- exp(l.mfec) # Mean productivity
mim <- exp(l.mim) # Mean immigration rate

# Population growth rate
for (t in 1:(nyears-1)){
lambda[t] <- Ntot[t+1] / Ntot[t]
logla[t] <- log(lambda[t])
}
mlam <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)])) # Geometric mean

# -----------------------------------------
# 4. The likelihoods of the single data sets
# -----------------------------------------
# 4.1. Likelihood for population count data (state-space model)
# 4.1.1 System process
for (t in 2:nyears){
mean1[t] <- 0.5 * f[t-1] * phij[t-1] * Ntot[t-1]
N1[t] ~ dpois(mean1[t])
NadSurv[t] ~ dbin(phia[t-1], Ntot[t-1])
mpo[t] <- Ntot[t-1] * omega[t-1]
Nadimm[t] ~ dpois(mpo[t])
}

# 4.1.2 Observation process
for (t in 1:nyears){
Ntot[t] <- NadSurv[t] + Nadimm[t] + N1[t]
y[t] ~ dpois(Ntot[t])
}

# 4.2 Likelihood for capture-recapture data: CJS model (2 age classes)
# Multinomial likelihood
for (t in 1:(nyears-1)){
marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t])
marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])
}

# Calculate number of released individuals
for (t in 1:(nyears-1)){
r.j[t] <- sum(marray.j[t,])
r.a[t] <- sum(marray.a[t,])
}

# m-array cell probabilities for juveniles
for (t in 1:(nyears-1)){
q[t] <- 1-p[t]

# Main diagonal
pr.j[t,t] <- phij[t]*p[t]

# Above main diagonal
for (j in (t+1):(nyears-1)){
pr.j[t,j] <- phij[t]*prod(phia[(t+1):j])*prod(q[t:(j-1)])*p[j]
}#j

# Below main diagonal
for (j in 1:(t-1)){
pr.j[t,j] <- 0
}#j

# Last column
pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
}#t

# m-array cell probabilities for adults
for (t in 1:(nyears-1)){

# Main diagonal
pr.a[t,t] <- phia[t]*p[t]

# above main diagonal
for (j in (t+1):(nyears-1)){
pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
}#j

# Below main diagonal
for (j in 1:(t-1)){
pr.a[t,j] <- 0
}#j

# Last column
pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
}#t

# 4.3. Likelihood for productivity data: Poisson regression
for (t in 1:(nyears-1)){
J[t] ~ dpois(rho[t])
rho[t] <- R[t] * f[t]
}
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(nyears = nyears,
                  marray.j = marray.j,
                  marray.a = marray.a,
                  y = popcount,
                  J = J,
                  R = R)
# Initial values
inits <- function(){list(l.mphij = rnorm(1, 0.2, 0.5),
                         l.mphia = rnorm(1, 0.2, 0.5),
                         l.mfec = rnorm(1, 0.2, 0.5),
                         l.mim = rnorm(1, 0.2, 0.5),
                         l.p = rnorm(1, 0.2, 1),
                         sig.phij = runif(1, 0.1, 10),
                         sig.phia = runif(1, 0.1, 10),
                         sig.fec = runif(1, 0.1, 10),
                         sig.im = runif(1, 0.1, 10),
                         N1 = round(runif(nyears, 1, 50), 0),
                         NadSurv = round(runif(nyears, 5, 50), 0),
                         Nadimm = round(runif(nyears, 1, 50), 0))}

# Parameters monitored
parameters <- c("phij",
                "phia",
                "f",
                "omega",
                "p",
                "lambda",
                "mphij",
                "mphia",
                "mfec",
                "mim",
                "mlam",
                "sig.phij",
                "sig.phia",
                "sig.fec",
                "sig.im",
                "N1",
                "NadSurv",
                "Nadimm",
                "Ntot")
# MCMC settings
ni <- 20000
nt <- 6
nb <- 10000
nc <- 3
# Call WinBUGS from R (BRT 5 min)
ipm.hoopoe <- bugs(bugs.data,
                   inits,
                   parameters, "ipm.hoopoe.bug",
                   n.chains = nc,
                   n.thin = nt,
                   n.iter = ni,
                   n.burnin = nb,
                   bugs.directory = bugs.dir,
                   #working.directory = getwd(),
                   debug = TRUE)
# Close the WinbUGS window & run the following command
print(ipm.hoopoe)
