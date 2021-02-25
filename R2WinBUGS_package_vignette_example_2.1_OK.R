#R2WinBUGS: A package for running WinBUGS from R
# Sturtz et al 2005
rm(list = ls())
setwd("C:/Users/Etudiant/Desktop/SMAC/Projet_publi/3-Modelisation_translocation_petrels/SCRIPTS/R2WinBUGS_package_vignette_example_2.1_OK")
# Example 2.1. School data
library(R2WinBUGS)
library(coda)
library(lme4)

data(schools)
# Specify model in BUGS language
sink("schools.bug") # store the model in separated file
cat("
model{
  for(j in 1:J){
    y[j] ~ dnorm(theta[j], tau.y[j])
    theta[j] ~ dnorm(mu.theta, tau.theta)
    tau.y[j] <- pow(sigma.y[j], -2)
  }
  mu.theta ~ dnorm(0.0, 1.0E-6)
  tau.theta <- pow(sigma.theta, -2)
  sigma.theta ~ dunif(0, 1000)
}
",fill = TRUE)
sink()

# Preparation of initial data
J <- nrow(schools)
y <- schools$estimate
sigma.y <- schools$sd
data <- list("J", "y", "sigma.y")

# How many chain to be run
n.chain <- 3
n.iter <- 1000

# Initial values
inits <- function(){
  list(theta = rnorm(J, 0, 100),
       mu.theta = rnorm(1, 0, 100),
       sigma.theta = runif(1, 0, 10))
}

# Start of the MCMC simulation

schools.sim <- bugs(data,
                    inits,
                    model.file = "C:/Users/Etudiant/Documents/schools.bug",
                    parameters = c("theta", "mu.theta", "sigma.theta"),
                    n.chains = n.chain,
                    n.iter = n.iter,
                    bugs.directory = "C:/Users/Etudiant/Documents/WinBUGS14/",
                    #working.directory = getwd(),
                    debug = T)

# Print results
print(schools.sim)
plot(schools.sim)
