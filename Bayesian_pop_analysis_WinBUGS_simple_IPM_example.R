# IPMs scripts 
# Bayesian population analysis using WinBUGS (book)
# Chapter 11 - Example of simple IPM (counts, CMR, reproduction)

rm(list = ls())
setwd("C:/Users/Etudiant/Desktop/SMAC/Projet_publi/3-Modelisation_translocation_petrels/SCRIPTS/Bayesian_pop_analysis_WinBUGS_simple_IPM_example")

# Load packages
library(R2WinBUGS)
library(coda)
library(lme4)
#library(rjags)

# Localization of WinBUGS.exe
bugs.dir <- "C:/Users/Etudiant/Documents/WinBUGS14/"

#### ----------- Load data ---------- ####
# Population counts (from years 1 to 10)
y <- c(45, 48, 44, 59, 62, 62, 55, 51, 46, 42)

# Capture-recapture data (in m-array format, from years 1 to 10)
m <- matrix(c(11, 0, 0, 0, 0, 0, 0, 0, 0, 70,
              0, 12, 0, 1, 0, 0, 0, 0, 0, 52,
              0, 0, 15, 5, 1, 0, 0, 0, 0, 42,
              0, 0, 0, 8, 3, 0, 0, 0, 0, 51,
              0, 0, 0, 0, 4, 3, 0, 0, 0, 61,
              0, 0, 0, 0, 0, 12, 2, 3, 0, 66,
              0, 0, 0, 0, 0, 0, 16, 5, 0, 44,
              0, 0, 0, 0, 0, 0, 0, 12, 0, 46,
              0, 0, 0, 0, 0, 0, 0, 0, 11, 71,
              10, 2, 0, 0, 0, 0, 0, 0, 0, 13,
              0, 7, 0, 1, 0, 0, 0, 0, 0, 27,
              0, 0, 13, 2, 1, 1, 0, 0, 0, 14,
              0, 0, 0, 12, 2, 0, 0, 0, 0, 20,
              0, 0, 0, 0, 10, 2, 0, 0, 0, 21,
              0, 0, 0, 0, 0, 11, 2, 1, 1, 14,
              0, 0, 0, 0, 0, 0, 12, 0, 0, 18,
              0, 0, 0, 0, 0, 0, 0, 11, 1, 21,
              0, 0, 0, 0, 0, 0, 0, 0, 10, 26),
            ncol = 10,
            byrow = TRUE) # The lastcolumn in matrix m contains the number of released individuals that are never recaptured. The top half of the array contains the data on birds marked as juveniles and the bottom half on those marked as adults.

# Productivity data (from years 1 to 9) - The numbers recorded in the last year are not considered here because they are not needed in the population model

J <- c(64, 132, 86, 154, 156, 134, 116, 106, 110) # Total number of nestlings recorded
R <- c(21, 28, 26, 38, 35, 33, 31, 30, 33) # Annual number of surveyed broods


#### ---------- Analysis of the model ---------- ####

# Specify model in BUGS language
sink("ipm.bug")
cat("
model {
# -----------------------------------------------
# Integrated population model
# - Age structured model with 2 age classes:
# 1-year old and adults (at least 2 years old)
# - Age at first breeding = 1 year
# - Prebreeding census, female-based
# - All vital rates assumed to be constant
# -----------------------------------------------

# -----------------------------------------------
# 1. Define the priors for the parameters
# -----------------------------------------------

# Observation error
tauy <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 50)
sigma2.y <- pow(sigma.y, 2)

# Initial population sizes
N1[1] ~ dnorm(100, 0.0001)I(0,) # 1-year
Nad[1] ~ dnorm(100, 0.0001)I(0,) # Adults

# Survival and recapture probabilities, as well as productivity
for (t in 1:(nyears-1)){
sjuv[t] <- mean.sjuv
sad[t] <- mean.sad
p[t] <- mean.p
f[t] <- mean.fec
}
mean.sjuv ~ dunif(0, 1)
mean.sad ~ dunif(0, 1)
mean.p ~ dunif(0, 1)
mean.fec ~ dunif(0, 20)

# -----------------------------------------------
# 2. Derived parameters
# -----------------------------------------------

# Population growth rate
for (t in 1:(nyears-1)){
lambda[t] <- Ntot[t+1] / Ntot[t]
}

# -----------------------------------------------
# 3. The likelihoods of the single data sets
# -----------------------------------------------

# 3.1. Likelihood for population population count data (state-space model)

# 3.1.1 System process
for (t in 2:nyears){
mean1[t] <- f[t-1] / 2 * sjuv[t-1] * Ntot[t-1]
N1[t] ~ dpois(mean1[t])
Nad[t] ~ dbin(sad[t-1], Ntot[t-1])
}
for (t in 1:nyears){
Ntot[t] <- Nad[t] + N1[t]
}

# 3.1.2 Observation process
for (t in 1:nyears){
y[t] ~ dnorm(Ntot[t], tauy)
}

# 3.2 Likelihood for capture-recapture data: CJS model (2 age classes)

# Multinomial likelihood
for (t in 1:2*(nyears-1)){
m[t,1:nyears] ~ dmulti(pr[t,], r[t])
}

# Calculate the number of released individuals
for (t in 1:2*(nyears-1)){
r[t] <- sum(m[t,])
}

# m-array cell probabilities for juveniles
for (t in 1:(nyears-1)){

  # Main diagonal
q[t] <- 1-p[t]
pr[t,t] <- sjuv[t] * p[t]

  # Above main diagonal
for (j in (t+1):(nyears-1)){
pr[t,j] <- sjuv[t]*prod(sad[(t+1):j])*prod(q[t:(j-1)])*p[j]
}#j

  # Below main diagonal
for (j in 1:(t-1)){
pr[t,j] <- 0
}#j

  # Last column: probability of non-recapture
pr[t,nyears] <- 1-sum(pr[t,1:(nyears-1)])
}#t

  # m-array cell probabilities for adults
for (t in 1:(nyears-1)){

  # Main diagonal
pr[t+nyears-1,t] <- sad[t] * p[t]

  # Above main diagonal
for (j in (t+1):(nyears-1)){
pr[t+nyears-1,j] <- prod(sad[t:j])*prod(q[t:(j-1)])*p[j]
}#j

  # Below main diagonal
for (j in 1:(t-1)){
pr[t+nyears-1,j] <- 0
}#j

  # Last column
pr[t+nyears-1,nyears] <- 1 - sum(pr[t+nyears-1,1:(nyears-1)])
}#t

# 3.3. Likelihood for productivity data: Poisson regression

for (t in 1:(nyears-1)){
J[t] ~ dpois(rho[t])
rho[t] <- R[t]*f[t]
}
}
",
fill = TRUE)

sink()

# Bundle data
bugs.data <- list(m = m,
                  y = y,
                  J = J,
                  R = R,
                  nyears = dim(m)[2])
# Initial values
inits <- function(){
  list(mean.sjuv = runif(1, 0, 1),
       mean.sad = runif(1, 0, 1),
       mean.p = runif(1, 0, 1),
       mean.fec = runif(1, 0, 10),
       N1 = rpois(dim(m)[2], 30),
       Nad = rpois(dim(m)[2], 30),
       sigma.y = runif(1, 0, 10))
}

# Parameters monitored
parameters <- c("mean.sjuv",
                "mean.sad",
                "mean.p",
                "mean.fec",
                "N1",
                "Nad",
                "Ntot",
                "lambda",
                "sigma2.y")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 2 min)
ipm <- bugs(bugs.data,
            inits,
            parameters,
            "ipm.bug",
            n.chains = nc,
            n.thin = nt,
            n.iter = ni,
            n.burnin = nb,
            bugs.directory = bugs.dir,
            #working.directory = NULL,
            debug = TRUE)

#### ---------- Close the WinbUGS window & run the following command ---------- ####
print(ipm, digits = 3)

#### ---------- Producing figure ---------- ####
par(cex = 1.2)
lower <- upper <- numeric()
for (i in 1:10){
  lower[i] <- quantile(ipm$sims.list$Ntot[,i], 0.025)
  upper[i] <- quantile(ipm$sims.list$Ntot[,i], 0.975)
}

plot(ipm$mean$Ntot,
     type = "b",
     ylim = c(35, 65),
     ylab = "Population size",
     xlab = "Year",
     las = 1,
     pch = 16,
     col = "blue",
     frame = F,
     cex = 1.5)
segments(1:10,
         lower,
         1:10,
         upper,
         col = "blue")
points(y,
       type = "b",
       col = "black",
       pch = 16,
       lty = 2,
       cex = 1.5)
legend(x = 1,
       y = 65,
       legend = c("Counts", "Estimates"),
       pch = c(16, 16),
       col = c("black", "blue"),
       lty = c(2, 1),
       bty = "n")
