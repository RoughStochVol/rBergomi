#####################################################
## Test poor man's parallelization
#####################################################

source("rbergomi_console.R")

## Generate test data
n <- 1000 ## number of data points
xi <- runif(n, min = 0.01, max = 1)
H <- runif(n, min = 0.01, max = 0.2)
eta <- runif(n, min = 1, max = 3)
rho <- runif(n, min = -0.95, max = -0.5)
K <- runif(n, min = 0.8, max = 1.2)
expiry <- runif(n, min = 0.02, max = 1.5)

## numerical parameters
N <- 200
M <- 10000 ## expect high calibration error anyway
num.jobs <- 8

##  run the code
time.system <- system.time(res <<- rBergomi.pricer2(xi, H, eta, rho, 
                                                    expiry, K, N, M, num.jobs))
time.system ## Takes 74 seconds (M = 10000, n=1000)

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
hist(log(TK.matrix[,1])) ## Reject logarithmic transform, as it might lead to too few small samples
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
hist(TK.sample(1000, fhat)[,1]) ## looks good

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
