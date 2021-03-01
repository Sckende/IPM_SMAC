#*************************************************************
#IPM for Emperor penguins
# M5: {phij(sam), phia(.), f(t), p(t)}
# Pre breeding matrix population model with 6 age classes (1,2,3,4,5 5+(ad))
# SAM is calculated for the rearing period (JULY-DEC)
#**************************************************************
#**********************************************
# Set working directory
#**********************************************

setwd("C:/CNRS/emperor penguin/2016/output/May2016/M5-1-sam")

#******************************
# bugs directory
#******************************
bugs.dir <- c("C:/WinBUGS14")

#******************************
# Load libraries
#******************************
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("lattice","coda","R2WinBUGS","R2jags","rjags",'mcmcplots','runjags')
ipak(packages)

#*************************************************************************
# Function to create a m-array based on capture-histories (CH)
#*************************************************************************

marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data=0, ncol=n.occasions+1, nrow= n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,] != 0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z], pos[z+1]] <- m.array[pos[z], pos[z+1]] +1
    } # z
  } # i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions +1] <- m.array[t,1] -
      sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}


#************************
# data
#************************

#***********************************
# Covariates 
#**********************************

# Southern annular mode (1969-2011) - calculated for the rearing period (JULY-DEC)

sam <- read.csv("sam_july-dec.csv",header=T)

sam <- sam[3:29,2]   #1971-1997

zsam <- (sam-mean(sam))/sd(sam)     # Standardized 


#***************************
# count data 1971-1998
#***************************

CountData <- read.csv("Counts_fecundity.csv",sep=";",header=T)
BRAD <- CountData[-c(1:19,48:60),2]     
CHICK <- floor(CountData[-c(1:19,48:60),3]) 
BRAD <- as.vector(BRAD)
CHICK <- as.vector(CHICK)
Time <- length(BRAD)[1]

#***************************
# data of mark recapture  1971-1998
#***************************

hcr.adu <-read.table("ME_adults_v2_1971-1998.csv",sep=";",header=T)
hcr.adu <- as.matrix(hcr.adu)
dimnames(hcr.adu) <- NULL

hcr <- hcr.adu   
DI <- dim(hcr)

#find the first capture
get.first <- function(x) min(which(x!=0))
f <- apply(hcr, 1, get.first)

HC <- hcr
HC <- as.matrix(HC)

m.HC <- marray(HC)

#****************************************
# Covariate model
# M5: {phij(sam), phia(.), f(t), p(t)}
# Fixed effects model
#*****************************************

sink("ipm_emp_6agc_M5.jags")

cat("
    model{
    
    
    #***************************************
    # Define priors and constraints for the parameters
    #****************************************
    
    # Initial population sizes
    
    n1 ~ dnorm(500,0.01)T(0,)
    N1[1] <- round(n1)
    n2 ~ dnorm(500,0.01)T(0,)
    N2[1] <- round(n2)
    n3 ~ dnorm(500,0.01)T(0,)
    N3[1] <- round(n3)
    n4 ~ dnorm(500,0.01)T(0,)
    N4[1] <- round(n4)
    n5 ~ dnorm(500,0.01)T(0,)
    N5[1] <- round(n5)
    na ~ dnorm(5000,0.01)T(0,) 
    Na[1] <- round(na)
    
    #Survival, resighting and fecundity
    
     for(k in 1:2){
     beta[k] ~  dunif(-5,5)  
     } #k
    
    sig.phi1 ~ dunif(0,5) ; tau.phi1 <- pow(sig.phi1,-2)
    
    mu.phia ~ dunif(0,1)
    
    mu.fec ~ dunif(0,1) ; lmu.fec <- log(mu.fec/(1-mu.fec)) 
    sig.fec ~ dunif(0,5) ; tau.fec <- pow(sig.fec,-2)

    mu.p ~ dunif(0,1) ; lmu.p <- log(mu.p/(1-mu.p)) 
    sig.p ~ dunif(0,5) ; tau.p <- pow(sig.p,-2)
    
    
    # proportion of breeders fixed (see Jenouvrier et al.2005 / CB)
    
    pb5 <- 0.22      # age4
    pb6 <- 0.32      # age5
    
    
    for (t in 1:(n.occasions-1)){
    
    # juvenile survival     phi1  modeled as a function of covariates

    logit(phi1[t]) <- lphi1.lim[t]

    lphi1.lim[t] <- min(999, max(-999,lphi1[t]))

    lphi1[t] <- beta[1] + beta[2]*zsam[t] + eps.phi1[t]

    eps.phi1[t] ~ dnorm(0, tau.phi1)T(-5,5)

    }#t
    

    for (t in 1:(n.occasions-1)){ 
    
    # adult survival
    
     phia[t] <- mu.phia
    
    # productivity

    logit(fec[t]) <- lfec.lim[t]

    lfec.lim[t] <- min(999, max(-999,lfec[t]))
    
    lfec[t] <- lmu.fec + eps.fec[t]
    
    eps.fec[t] ~ dnorm(0, tau.fec)T(-5,5)

    #resighting
    
    logit(p[t]) <- lp.lim[t]
    
    lp.lim[t] <- min(999, max(-999,lp[t]))
    
    lp[t] <- lmu.p + eps.p[t]   

    eps.p[t] ~ dnorm(0, tau.p)T(-5,5)
  
    } #t
   

    #********************************
    # Derived parameters
    #********************************
    
    # mean juv. survival
    
    mu.phi1 <- exp(beta[1])/(1+exp(beta[1]))

    #***********************************
    # Likelihood for each data set
    #***********************************
    
    #***************************************
    # Likelihood for population population count data (state-space model)
    #****************************************
    
    #**********************
    # System process
    #***********************
    for (t in 2:Time){
    # Define the system process for the census data
    
    #age1
    mean1[t] <- 0.5*Na[t-1]*fec[t-1]*phi1[t-1]        
    N1[t] ~ dpois(mean1[t])   
    
    #age2
    N2[t] ~ dbin(phia[t-1],N1[t-1])
    
    #age3
    N3[t] ~ dbin(phia[t-1],N2[t-1])
    
    #age4
    N4[t] ~ dbin(phia[t-1],N3[t-1])
    
    #age5 
    mean5[t] <- round((1-pb5)*N4[t-1] + (1-pb6)*N5[t-1])
    N5[t] ~ dbin(phia[t-1],mean5[t])
    
    #breeding adults
    meana[t] <- round(pb5*N4[t-1] + pb6*N5[t-1] + Na[t-1]) 
    Na[t] ~ dbin(phia[t-1],meana[t])
    }
    
    for (t in 1:Time){
    Ntot[t] <- Na[t] + N5[t] + N4[t] + N3[t] + N2[t] + N1[t]
    } #t
    
    #******************************
    # Observation process
    #******************************
    
    for (t in 1:Time){
    Ka[t] ~ dpois(Na[t])
    }#t
    
    #*************************************
    # Likelihood for capture-recapture data:  m-array
    #****************************************
    
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
    marr[t,1:n.occasions] ~ dmulti(pr[t, ],r[t])
    }
    # Calculate the number of birds released each year
    #for (t in 1:(n.occasions-1)){
    #r[t] <- sum(marr[t, ])
    #}
    # Define the cell probabilities of the m-array
    # Main diagonal
    for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t] # Probability of non-recapture
    pr[t,t] <- phia[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } # j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr[t,j] <- 0
    } # j
    } # t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
    pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
    } # t
    
    #**************************************************
    # Likelihood for reproductive success data
    #**************************************************
    for (t in 1:(Time-1))
{
    Kp[t] ~ dbin(fec[t],Ka[t])
}#t
    
    
    } # End Model
    ",fill = TRUE)
sink()

#*********************
# Bundle data
#*********************

bugs.data <- list(n.occasions=dim(m.HC)[2],Ka=BRAD,Kp=CHICK,Time=length(BRAD),marr=m.HC,r= rowSums(m.HC),zsam=zsam) #,zsam=zsam


#**********************
# Initial values         ,
#*********************
#

inits <- function(){list(beta =runif(2,-2,2),mu.phia=runif(1,0.8,1),mu.fec=runif(1,0,1),mu.p=runif(1,0,1),
                         sig.phi1 =runif(1,0,5),sig.fec =runif(1,0,5),sig.p =runif(1,0,5),  
                         n1=rep(200,1),n2=rep(200,1),n3=rep(100,1),n4=rep(200,1),n5=rep(300,1),na=rep(2500,1))} #,w=rbinom(1,1,0.5)

#*********************************
# parameters to be monitored            
#**********************************

parameters <- c("beta","mu.phi1","mu.phia","mu.fec", "mu.p",
                "sig.phi1","sig.fec","sig.p")

#**************************
# MCMC settings
#***************************

ni <- 2000000
nb <- 1000000
nt <- 100
nc <- 3


#*******************************************
# Do the MCMC stuff calling WinBUGS from R
#*******************************************

start <- as.POSIXlt(Sys.time())
print(start)

ipm_emp_6agc_M5 <- jags(data = bugs.data, inits = inits, parameters.to.save = parameters, model.file = "ipm_emp_6agc_M5.jags",
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        working.directory=getwd())

end <- as.POSIXlt(Sys.time())

print(end)

duration <- print(difftime(end,start,units='mins'))


#**************************
# Summary
#*************************

print(ipm_emp_6agc_M5,digits=3)
