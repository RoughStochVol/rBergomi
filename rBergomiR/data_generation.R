##################################################
## Data generation for ML
##################################################
## Provide data for the deep network for learning
## rough Bergomi.
## Code is provided here as many parameters enter
## the code directly.

source("rbergomi_console.R")

########################################################
## Use market data for (K,T)
########################################################

library(ks)
m.data <- read.csv(file = "christian.csv")
TK.matrix <- cbind(m.data$time.to.maturity..years., m.data$moneyness)
colnames(TK.matrix) <- c("T","K")
## Use open interest as rates, but capped at 1 and 400 
weights <- m.data$Open.Int
weights[weights > 400] <- 400
weights[weights < 1] <- 1
bandwidth <- Hpi(TK.matrix)
fhat <- kde(TK.matrix, bandwidth, w = weights)

## Generate sample data, conditioned on expiry > 0.
TK.sample <- function(M, fhat){
  library(ks)
  x <- rkde(M, fhat)
  while(sum(x[,1] < 0) > 0)
    x[x[,1] < 0,] <- rkde(sum(x[,1] < 0), fhat)
  return(x)
}

#################################################################
## Now use transformed, cut data with inverse spread as weights
## and higher bandwidth.
#################################################################
library(ks)
m.data <- read.csv(file = "processed_spx_data.csv")
TK.matrix <- cbind(m.data$time.to.maturity..years., m.data$moneyness)
colnames(TK.matrix) <- c("T","K")
## Use inverse spread, but capped at 8 
weights <- m.data$inv_spread
weights[weights > 8] <- 8
bandwidth <- Hpi(TK.matrix)
bandwidth <- 10*bandwidth ## artificially increase bandwidth to smooth further
fhat <- kde(TK.matrix, bandwidth, w = weights)
plot(fhat, cont = 10*1:10)
points(TK.matrix[,1], TK.matrix[,2])

## Generate sample data, conditioned on expiry > 0.
TK.sample <- function(M, fhat){
  library(ks)
  x <- rkde(M, fhat)
  while(sum(x[,1] < 0) > 0)
    x[x[,1] < 0,] <- rkde(sum(x[,1] < 0), fhat)
  return(x)
}


#########################################################
## Use "reasonable" distributions for other parameters
#########################################################
## Use truncated normal distributions for other parameters,
## possibly up to transformations
library(truncnorm)
rho.sample <- function(M) rtruncnorm(M, a = -1, b = -0.5, mean = -0.95, sd = 0.2)
eta.sample <- function(M) rtruncnorm(M, a = 1, b = 4, mean = 2.5, sd = 0.5)
H.sample <- function(M) rtruncnorm(M, a = 0.01, b = 0.5, mean = 0.07, sd = 0.05)
xi.sample <- function(M) rtruncnorm(M, a = 0.05, b = 1, mean = 0.3, sd = 0.1)^2

#########################################################
## Put all samplers together
#########################################################
par.sample <- function(M){
  TK <- TK.sample(M, fhat)
  rho <- rho.sample(M)
  eta <- eta.sample(M)
  H <- H.sample(M)
  xi <- xi.sample(M)
  par <- cbind(xi, H, eta, rho, TK)
  colnames(par) <- c("xi", "H", "eta", "rho", colnames(TK))
  return(data.frame(par))
}

#########################################################
## Provide wrapper function around rBergomi.pricer2
#########################################################
pricer <- function(par.data, num.steps, num.MC.samples, num.jobs){
  time.system <- system.time(res <<- rBergomi.pricer2(par.data$xi, par.data$H, 
                                                      par.data$eta, par.data$rho, 
                                                      par.data$T, par.data$K, num.steps, 
                                                      num.MC.samples, num.jobs))
  num.stat <- c(num.steps, num.MC.samples, num.jobs)
  names(num.stat) <- c("steps", "MC.samples", "num.jobs")
  num.stat <- c(num.stat, time.system)
  return(list(res = res, num.stat = num.stat))
}

############################################################
## Provide all numerical parameters
############################################################
# ## This is a trial run on 1 process to get an idea of the runtime
# ## required.
# 
# num.data <- 10000 ## number of prices to be provided
# num.steps <- 200 ## number of timesteps for Euler discretization
# num.MC.samples <- 1000 ## number of MC samples in the pricing routine
# num.jobs <- 30 ## number of parallel threads
# 
# ## Sample the parameters
# par.data <- par.sample(num.data)
# 
# ## Call the pricer (and measure timing)
# dat <- pricer(par.data, num.steps, num.MC.samples, num.jobs)
# 
# ## On erhard-06: steps = 200, MC.samples = 1000, num.jobs = 28, num.data = 10000 takes 45 seconds.

## This is a serious run.

num.data <- 1000000 ## number of prices to be provided
num.steps <- 200 ## number of timesteps for Euler discretization
num.MC.samples <- 50000 ## number of MC samples in the pricing routine
num.jobs <- 30 ## number of parallel threads

## Sample the parameters
par.data <- par.sample(num.data)

## save for Ben
head(par.data[,5:6])
write.csv(par.data[,5:6], file = "sampled.csv")

## Call the pricer (and measure timing)
dat <- pricer(par.data, num.steps, num.MC.samples, num.jobs)

save(dat, file = "rBergomi03.RData")

load("rBergomi03.RData")
write.csv(dat$res, file = "rBergomi03.csv")
