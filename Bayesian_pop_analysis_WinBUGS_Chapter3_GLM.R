# Bayesian population analysis using WinBUGS #
#### Chapter 3 - Intro to the GLM: The simplest model for count data ####

# how to represent a particular statistical distribution
# ?dnorm, ?dpois, ?dbinom, ?dmultinom, ?dunif
plot(density(rbeta(n = 10^6, shape1 = 2, shape2 = 4)))
hist(rbeta(10^6, 2, 4), nclass = 100, col = "gray")

plot(density(rnorm(n = 10000, mean = 10, sd = 3)))
hist(rnorm(n = 10000, mean = 10, sd = 3), nclass = 100, col = "grey")

plot(density(rpois(100, 5)), ylim = c(0,0.5))
lines(density(rpois(100, 2)), col = "red")
lines(density(rpois(100,1)), col = "green")


### Example with ANCOVA linear model ###
# Toy data set of nine data points
# Y = response
# A = factor with three levels
# X = continuous covariate

# Define and plot data
Y <- c(25, 14, 68, 79, 64, 139, 49, 119, 111) 
A <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))
X <- c(1, 14, 22, 2, 9, 20, 2, 13, 22)

plot(X, Y,
     col = c(rep("red", 3), rep("blue", 3), rep("green", 3)),
     xlim = c(-1, 25),
     ylim = c(0, 140))

# With R, we can fit an ANCOVA with parallel slopes
summary(fm <- lm(Y ~ A-1 + X))

# To obtain the design-matrix of the model with R
# Effects or treatment contrast parameterization
model.matrix(~ A + X) # the linear model in terms of a baseline response, which here is that for the first level of factor A, plus effects of the other levels of A relative to the first level and effects of each unit change in the covariate X
summary(fm <- lm(Y ~ A + X)) #the intercept is the mean response of units with level 1 of factor A at X = 0, and the parameters corresponding to the design matrix columns A2 and A3 quantify the difference in the mean response in levels 2 and 3, relative to that of level 1, for a given value of X.

# Or ...
# Means parameterization
model.matrix(~ A - 1 + X) # in the means parameterization, the first three parameters directly represent the mean response for each level of factor A at X = 0, while the meaning of the parameter represented by the column X is the same as before.
summary(fm <- lm(Y ~ A-1 + X))

# To obtain the values of the linear predictor or expected or "typical" response of this kind of model
model.matrix(~A-1+X) %*% fm$coef

### Generation and analysis of simulated data with Poisson GLM
# Definition of  a function that generates Poisson counts of peregrine falcons
# During n years
# Inspired by Monneret, 2006
# Count will follow a cubic polynomial function of time
# nu_i = alpha + beta1*Xi + beta2*Xi^2 + beta3*Xi^3
# 2 methods used here : frequentist (R) vs. Bayesian (WinBUGS)

# ----- Generation of simulated data

data.fn <- function(n, alpha = 3.5576, beta1 = -0.0912, beta2 = 0.0091, beta3 = -0.00014){
  # n: Number of years
  # alpha, beta1, beta2, beta3: coefficients of a cubic polynomial of count on year
  
  # Generate values of time covariate
  year <- 1:n
  
  # Signal: Build up systematic part of the GLM
  log.expected.count <- alpha + beta1*year + beta2*year^2 + beta3*year^3 # Equivalent to nu_i
  expected.count <- exp(log.expected.count) # Back-transformation for obtaining true values and what we wanted to find when modelling
  
  # Noise: Generate random part of GLM - Poisson noise around expected count
  C <- rpois(n = n, lambda = expected.count) # Data used for modelisation of peregrin counts
  
  # Plot simulated data
  plot(year,
       C,
       type = "b",
       lwd = 2,
       col = "black",
       main = "",
       las = 1,
       ylab = "Population size",
       xlab = "Year",
       cex.lab = 1.2,
       cex.axis = 1.2)
  lines(year,
        expected.count,
        type = "l",
        lwd = 3,
        col = "red")
  
  return(list(n = n,
              alpha = alpha,
              beta1 = beta1,
              beta2 = beta2,
              beta3 = beta3,
              year = year,
              expected.count = expected.count,
              C = C))
}

data <- data.fn(n = 40)
data

# ----- Analysing with frequentist framework

fm <- glm(formula = C ~ year + I(year^2) + I(year^3),
          family = poisson,
          data = data)
summary(fm)

# ----- Analysing with Bayesian framework & BUGS

# Specify model in BUGS language
sink("GLM_Poisson.txt")
cat("
model{

# Priors
alpha ~ dunif(-20, 20)
beta1 ~ dunif(-10, 10)
beta2 ~ dunif(-10, 10)
beta3 ~ dunif(-10, 10)

# Likelihood: Note key components of a GLM on one line each
for( i in 1:n){
  C[i] ~ dpois(lambda[i])   # 1. Distribution for random part
  log(lambda[i]) <- log.lambda[i]   # 2. Link function
  log.lambda[i] <- alpha + beta1*year[i] + beta2*pow(year[i], 2) + beta3*pow(year[i], 3)    # 3. Linear predictor
} # i
}
    ", fill = TRUE)
sink()

# Localization of WinBUGS.exe
bugs.dir <- "C:/Users/Etudiant/Documents/WinBUGS14/"

# Bundle into an R list the data needed for the analysis by WinBUGS
# win.data <- list(C = data$C,
#                  n = length(data$C),
#                  year = data$year) # *** DATA NEED TO BE STANDARDIZED cause value of parameter too far from 0 (year^3 = 40^3 = 64000)

# Standardization of data
mean.year <- mean(data$year)
sd.year <- sd(data$year)
win.data <- list(C = data$C,
                 n = length(data$C),
                 year = (data$year - mean.year)/sd.year) 

# Initials values - see comments in page 59
inits <- function(){
  list(alpha = runif(1, -2, 2),
       beta1 = runif(1, -3, 3))
}
# List of quantities we want to estimate - Parameters monitored
  params <- c("alpha", "beta1", "beta2", "beta3", "lambda")
  
# Before running the analysis, we set the MCMC characteristics: the number of draws per chain, thinning rate, burnin length, and the number of  chains. The burnin should be long enough to discard the initial part of  the Markov chains that have not yet converged to the stationary distribution.  We typically determine the required burnin in initial trials. Thinning  is useful to save computer disk space. We run multiple chains to check for  convergence.
ni <- 2000
nt <- 2
nb <- 1000
nc <- 3

# To observe a non convergence of chains
# ni <- 60
# nt <- 1
# nb <- 1
# nc <- 3

# Call WinBUGS from R
require(R2WinBUGS)
out <- R2WinBUGS::bugs(data = win.data,
                       inits = inits,
                       parameters.to.save = params,
                       model.file = "GLM_Poisson.txt",
                       n.chains = nc,
                       n.thin = nt,
                       n.iter = ni,
                       n.burnin = nb,
                       debug = T, # WinBUGS then remains open after the                       requested number of iterations has been produced, and we can visually inspect whether the chains have converged, or in the case of an error, directly read the log file
                       bugs.directory = bugs.dir,
                       working.directory = getwd())
# Summarize posteriors
print(out, digit = 3)
plot(out) # Graphical summary of Bayesian analysis + ***check the graphics in WinBUGS***

# What we learn from analysis of the model
# plot the Poisson means (the lambda parameters) for each year; these represent the expected peregrine counts in each year.
# And comparison of frequentist predictions vs. Bayesian pred.

R.predictions <- predict(fm, type = "response")
WinBUGS.predictions <- out$mean$lambda

plot(1:40,
     data$C,
     type = "b",
     lwd = 2,
     col = "black",
     main = "",
     las = 1,
     ylab = "Population size",
     xlab = "Year")
lines(1:40,
      R.predictions,
      type = "l",
      lwd = 3,
      col = "green")
lines(1:40,
      WinBUGS.predictions,
      type = "l",
      lwd = 3,
      col = "blue",
      lty = 2)
