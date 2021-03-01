rm(list = ls())
getwd()
setwd("C:/Users/Etudiant/Desktop/SMAC/Projet_publi/3-Modelisation_translocation_petrels/DATA/2017_Saunders_et_al_DATA")

# Localization of WinBUGS.exe
bugs.dir <- "C:/Users/Etudiant/Documents/WinBUGS14/"

sink("ipm.pva.jags")
cat("
model {
    #----------------------------------------------------------------------------
    #  Integrated population model BPVA (10 yr predictions)
    #  - Stage structured model with 2 stages: juvenile and adult
    #  - Age at first breeding = 1 year
    #  - Prebreeding census, female-based
    #  - All vital rates assumed to be time-dependent (random env. stochasticity)
    #  - Includes env. stochasticity thru random time effects for all params
    #  - Explicit estimation of immigration as expected number of individuals
    #  - Merlin effect (latent abundance) on adult survival estimated by state-
    # space model
    #----------------------------------------------------------------------------
    
    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes
    n1 ~ dnorm(100, 0.001)I(0,)           # HY individuals
    nadSurv ~ dnorm(100, 0.001)I(0,)      # Adults >= 2 years
    nadimm ~ dnorm(100, 0.001)I(0,)       # Immigrants
    
    N1[1] <- round(n1)
    NadSurv[1] <- round(nadSurv)
    Nadimm[1] <- round(nadimm)
    Ntot[1] <- N1[1] + NadSurv[1] + Nadimm[1]
    
    # Mean demographic parameters (on appropriate scale)
    l.mphij ~ dnorm(0, 0.001)              
    l.mphia ~ dnorm(0, 0.001)
    l.mfec ~ dnorm(0, 0.001)
    b0.omm ~ dunif(0, 20)  #expected number of immigrants              
    l.p ~ dnorm(0, 0.001)
    beta.phia ~ dnorm(0, 0.1)               
    
    #back transformation
    log.b0.omm <- log(b0.omm)
    
    
    # Precision of standard deviations of temporal variability
    sig.phij ~ dunif(0, 10)
    tau.phij <- pow(sig.phij, -2)
    sig.phia ~ dunif(0, 10)                  
    tau.phia <- pow(sig.phia, -2)
    sig.fec ~ dunif(0, 10)
    tau.fec <- pow(sig.fec, -2)
    sig.im ~ dunif(0, 10)
    tau.im <- pow(sig.im, -2)
    
    sig.obs ~ dunif(0.5, 50)
    tau.obs <- pow(sig.obs, -2)
    
    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1+K)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)	
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.fec[t] ~ dnorm(0, tau.fec)T(-5,5)
    epsilon.im[t] ~ dnorm(0, tau.im)T(-5,5)   
    }
    
    #-----------------------------------------------
    # 2. Constrain parameters (for temp variability)
    #-----------------------------------------------
    
    # Juv. apparent survival
    for (t in 1:(nyears-1+K)){
    logit(phij[t]) <- l.mphij + epsilon.phij[t]          
    
    # Adult apparent survival with merlin effect
    logit(phia[t]) <- l.mphia + beta.phia*N.cor[t] + epsilon.phia[t]    
    
    log(f[t]) <- l.mfec + epsilon.fec[t] # Productivity
    log(omega[t]) <- log.b0.omm + epsilon.im[t] # Immigration
    logit(p[t]) <- l.p # Recapture probability
    }
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------

    mphij <- exp(l.mphij)/(1+exp(l.mphij))   # Mean juvenile survival probability
    mphia <- exp(l.mphia)/(1+exp(l.mphia))   # Mean adult survival probability
    mfec <- exp(l.mfec)                      # Mean productivity
    
    # Population growth rate (1993 to 2016)
    for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / (Ntot[t] + 0.0001)                             
    logla[t] <- log(lambda[t])
    }
    mlam <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)])) # Geometric mean
    
    
    # Population growth rate (merlins)
    for (t in 1:(nyears-1+K)){
    lambda.mer[t] <- N.est[t+1] / (N.est[t] + 0.0001)    
    logla.mer[t] <- log(lambda.mer[t])
    }
    
    mlam.mer <- exp((1/(nyears-1+K))*sum(logla.mer[1:(nyears-1+K)])) 
    # Geometric mean (all years)

    mlampast.mer <- exp((1/(nyears-1))*sum(logla.mer[1:(nyears-1)]))   	
    # Geometric mean for 1993-2015
    
    mlamfut.mer <- exp((1/(K-1))*sum(logla.mer[nyears:(nyears-1+K)]))   	  
    # Geometric mean for 2016-2026 

    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------

    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    for (t in 2:(nyears+K)){
    Ntot[t] <- NadSurv[t] + N1[t] + Nadimm[t]
    mean1[t] <- 0.5 * f[t-1] * phij[t-1] * Ntot[t-1]        
    N1[t] ~ dpois(mean1[t])
    NadSurv[t] ~ dbin(phia[t-1], Ntot[t-1])
    Nadimm[t] ~ dpois(omega[t-1])
    }
    
    # 4.1.2 Observation process
    for (t in 1:nyears){
    y[t] ~ dnorm(Ntot[t], tau.obs) 
    }
    
    # 4.2 Likelihood for capture-recapture data: CJS model (2 age classes)
    # Multinomial likelihood
    for (t in 1:(nyears-1)){
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t])
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-1)){
    q[t] <- 1-p[t]
    # Main diagonal
    pr.j[t,t] <- phij[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.j[t,j] <- phij[t]*prod(phia[(t+1):j])*prod(q[t:(j-1)])*p[j]
    } #j

    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults
    for (t in 1:(nyears-1)){
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    
    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t] * f[t]
    }

    #-------------------------------------------------------------------
    #  5. State-space model for merlin index (effect on adult survival)
    #-------------------------------------------------------------------

    # Priors and contraints
    logN.est[1] ~ dnorm(4.4, 0.01) # Prior for initial population size
    
    mean.r ~ dnorm(0, 0.01) # Prior for mean grown rate
    
    sigma.proc ~ dunif(0, 1) # Prior for SD of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    
    sigma.obs ~ dunif(0, 1) # Prior for SD of obs. process
    sigma2.obs <- pow(sigma.obs, 2)
    t.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){ # T is 34 years (24 + 10 prediction yrs)
    r[t] ~ dnorm(mean.r, tau.proc)
    logN.est[t+1] <- logN.est[t] + r[t]
    }
    
    # Observation process                     
    for (t in 1:T) {
    for (s in 1:S){
    x[t,s] ~ dnorm(logN.est[t], t.obs)
    }
    }
    
    # Population sizes on real scale
    for (t in 1:T) {                           
    N.est[t] <- exp(logN.est[t])
    N.cor[t] <- (N.est[t]-N.mean)/N.sd #standardize to use as covariate
    }
    }
    ",fill = TRUE)
sink()

#################################################################################
# Load packages---------------------------------------------------------------------
library(R2WinBUGS)
#library(rjags)
#library(jagsUI)
#library(lme4)
# Load data---------------------------------------------------------------------

M <- read.table("merlins.txt",header=TRUE)    #Hawk Mtn. and Whitefish Pt. counts

# First, alter data input
hawk <- c(M$HM)           #  Hawk Mountain counts
white <- c(M$WP)          #  Whitefish Point counts
mat <- matrix(c(hawk, white), nrow=length(hawk))

# CMR - Banded as juvenile
marray.j <- read.table("PIPL_asJUV_CMR.txt", header = F)
marray.j <- as.matrix(marray.j,
                      ncol = 24)

# CMR - Banded as adult
marray.a <- read.table("PIPL_asADULT_CMR.txt", header = F)
marray.a <- as.matrix(marray.a,
                      ncol = 24)
# **** Matrices need to be formatted for analyses - cf examples in Bayesian Population analyses  **** 

# Population count data
y <- read.table("PIPL_count.txt", h = T)
y <- as.vector(y[,2])

# Productivity data
prod <- read.table("PIPL_fecundity.txt", h = T)
J <- prod$FLEDGE # Number of offspring
R <- prod$SURV_BROOD # Number of surveyed broods
#-------------------------------------------------------------------------------
# Localization of WinBUGS.exe
bugs.dir <- "C:/Users/Etudiant/Documents/WinBUGS14/"

# Bundle data
K <- 10                           # Number of years with predictions
nyears <- ncol(marray.j) # Number of study years
N.mean = 114.9                    # Mean estimate of merlin population size
N.sd = 17.3                       # SD of merlin population size

# adjust merlin matrix for prediction years with NAs
v1 <- c(rep(NA, K))
v2 <- c(rep(NA, K))
mat.add <- matrix(c(v1, v2), nrow=length(v1))
mat.proj <-rbind(mat, mat.add)

jags.data <- list(nyears = nyears,
                  marray.j = marray.j,
                  marray.a = marray.a,
                  y = y,
                  J = J,
                  R = R,
                  r.j = rowSums(marray.j),
                  r.a = rowSums(marray.a),
                  K = K,
                  x = log(mat.proj),
                  T = nrow(mat.proj),
                  S = ncol(mat.proj),
                  N.mean = N.mean,
                  N.sd = N.sd)

# Initial values
inits <- function(){
    list(l.mphij = rnorm(1, 0.2, 0.5),
         l.mphia = rnorm(1, 0.2, 0.5),
         l.mfec = rnorm(1, 0.2, 0.5),
         l.p = rnorm(1, 0.2, 1),
         sig.phij = runif(1, 0.1, 10),
         sig.phia = runif(1, 0.1, 10),
         sig.fec = runif(1, 0.1, 10),
         n1 = round(runif(1, 1, 50), 0),
         nadSurv = round(runif(1, 5, 50), 0),
         beta.phia = runif(1, -1, 1),
         b0.omm = runif(1, 0, 10),
         sig.im = runif(1, 0.1, 10),
         nadimm = round(runif(1, 1, 50), 0),
         sigma.proc = runif(1, 0, 1),
         mean.r = rnorm(1),
         sigma.obs = runif(1, 0, 1),
         logN.est = c(rnorm(1, 4.4, 0.1),
                      rep(NA, (nrow(mat.proj) - 1))))
    }

# Parameters monitored
parameters <- c("phij",
                "phia",
                "f",
                "p",
                "lambda",
                "mphij",
                "mphia",
                "mfec",
                "mlam",
                "mlam.mer",
                "mlampast.mer",
                "mlamfut.mer",
                "beta.phia",
                "sig.phij",
                "sig.phia",
                "sig.fec",
                "sig.obs",
                "omega",
                "sig.im",
                "N1",
                "NadSurv",
                "Ntot",
                "Nadimm",
                "b0.omm",
                "r",
                "mean.r",
                "sigma2.obs",
                "sigma2.proc",
                "N.cor",
                "N.est")

# MCMC settings 
ni <- 400000    
nt <- 10
nb <- 200000
nc <- 3

# Call JAGS from R (jagsUI)
ipm.pva <- jagsUI::jags(jags.data,
                inits,
                parameters,
                "ipm.pva.jags",
                n.chains = nc,
                n.thin = nt,
                n.iter = ni,
                n.burnin = nb,
                parallel=FALSE,
                store.data=TRUE,
                bugs.format = TRUE)
