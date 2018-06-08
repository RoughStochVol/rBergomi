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

# library(ks)
# m.data <- read.csv(file = "christian.csv")
# TK.matrix <- cbind(m.data$time.to.maturity..years., m.data$moneyness)
# colnames(TK.matrix) <- c("T","K")
# ## Use open interest as rates, but capped at 1 and 400 
# weights <- m.data$Open.Int
# weights[weights > 400] <- 400
# weights[weights < 1] <- 1
# bandwidth <- Hpi(TK.matrix)
# fhat <- kde(TK.matrix, bandwidth, w = weights)
# 
# ## Generate sample data, conditioned on expiry > 0.
# TK.sample <- function(M, fhat){
#   library(ks)
#   x <- rkde(M, fhat)
#   while(sum(x[,1] < 0) > 0)
#     x[x[,1] < 0,] <- rkde(sum(x[,1] < 0), fhat)
#   return(x)
# }

#################################################################
## Now use transformed, cut data with inverse spread as weights
## and higher bandwidth.
#################################################################
library(ks)
m.data <- read.csv(file = "spx_liquid_calls.csv")
TK.matrix <- cbind(m.data$time.to.maturity..years., exp(m.data$log_moneyness))
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
num.steps <- 100 ## number of timesteps for Euler discretization
num.MC.samples <- 400000 ## number of MC samples in the pricing routine
num.jobs <- 20 ## # of jobs sequentially worked
num.threads <- 40 ## number of parallel threads

## function for providing file names
fname.fun <- function(i) paste("rBergomi04", i, sep = "_")

## Sample the parameters
par.data <- par.sample(num.data)

## save for Ben
#head(par.data[,5:6])
#write.csv(par.data[,5:6], file = "sampled_liquid.csv")

## Call the pricer on the i'th slice of data
count <- num.data / num.jobs
pricer.slice <- function(i){
  print(paste("Beginning of job", i))
  input <- par.data[((i-1)*count+1):(i*count),]
  dat <- pricer(input, num.steps, num.MC.samples, num.threads)
  save(dat, file = paste(fname.fun(i),".RData", sep = ""))
  write.csv(dat$res, file = paste(fname.fun(i),".csv", sep = ""))
  print(paste("End of job", i))
}

## Call the pricer (and measure timing)
dat.complete <- lapply(1:num.jobs, pricer.slice)
save(dat.complete, file = "rBergomi04.RData")

#####################################################################
## Plot a surface.
#####################################################################

K <- exp(seq(-0.1, 0.25, length.out = 50))
T <- seq(0.01, 0.5, length.out = 50)
K.vec <- expand.grid(K, T)[,1]
T.vec <- expand.grid(K, T)[,2]
xi <- rep(0.04, length(K.vec))
H <- rep(0.07, length(K.vec))
eta <- rep(2.2, length(K.vec))
rho <- rep(-0.9, length(K.vec))

num.steps <- 100 ## number of timesteps for Euler discretization
num.MC.samples <- 400000 ## number of MC samples in the pricing routine
num.threads <- 1 ## number of parallel threads

res <- rBergomi.pricer2(xi, H, eta, rho, T.vec, K.vec, num.steps, 
                 num.MC.samples, num.threads)
save(res, file="smile_res.RData")

## Plot the surface
persp(log(K), T, matrix(res$iv, ncol = length(K), byrow = FALSE), theta = -10,
      ticktype = "detailed", xlab = "log-moneyness", ylab = "maturity", zlab = "implied vol")

library(plot3D)
persp3D(log(K), T, matrix(res$iv, ncol = length(K), byrow = FALSE), 
      ticktype = "detailed", xlab = "log-moneyness", ylab = "maturity", zlab = "implied vol")
