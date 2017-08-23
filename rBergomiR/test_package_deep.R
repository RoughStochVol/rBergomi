####################################################################
## Deep test of RBergomi package
####################################################################

library(RBergomi)
source("test_package_aux.R")

## Test the RBergomi package against the old R-code.
xi <- 0.04
H <- 0.06
eta <- 2.5
rho <- -0.8

param <- list(xitu = function(t) rep(xi, length(t)), H = H, eta = eta, rho = rho)

## generate RBergomi samples using Cholesky factorization
n <- 11 ## number of attempts
tFinal <- 0.8
tsteps <- 2^(1:n)
paths <- 10*tsteps^2
K <- 1.0
#system.time(s <<- gen.rbergomi(param, tFinal, tsteps, paths)$s[tsteps,])
#res.chol <- get.price(s, K)

## generate RBergomi samples using the hybrid method in R
param <- c(xi, H, eta, rho)
Mc <- hybridScheme.vec(param)
#system.time(s <<- Mc(paths, tsteps, tFinal))
#res.hyb <- get.price(s, K)
## On erhard-06:
## Cholesky: price = 0.05784862, stat = 0.00011394, tsteps = 1000, paths = 100000
## Hybrid: price = 0.05761252, stat = 0.0001136558, tsteps = 1000, paths = 100000

## use the library, without RT trick
price <- rep(NA, n)
stat <- rep(NA, n)
time <- rep(NA, n)

for(i in 1:n){
   my.time <- system.time(res.C <<- computeIV(xi, H, eta, rho, tFinal, K, tsteps[i], paths[i], 
                                              c(12534+i, 45567+i, 53534+i, 25345+i), 4))
   stat[i] <- res.C$stat
   price[i] <- res.C$price
   time[i] <- my.time[3]
   cat(i,c(price[i], stat[i], tsteps[i], paths[i], time[i]),'\n')
}

res.C <- data.frame(price=price, stat=stat, time=time, tsteps=tsteps, paths=paths)
#save(res.C, file="res_C.RData")
#load(file="res_C.RData")

# ## use the library, with RT trick
# n <- 11
# price.RT <- rep(NA, n)
# stat.RT <- rep(NA, n)
# time.RT <- rep(NA, n)
# 
# for(i in 1:n){
#   my.time <- system.time(res.C <<- computeIV_RT(xi, H, eta, rho, tFinal, K, tsteps[i], paths[i], 
#                                              c(12534+i, 45567+i, 53534+i, 25345+i), 4))
#   stat.RT[i] <- res.C$stat
#   price.RT[i] <- res.C$price
#   time.RT[i] <- my.time[3]
#   cat(i,c(price.RT[i], stat.RT[i], tsteps[i], paths[i], time.RT[i]),'\n')
# }
# 
# res.RT.C <- data.frame(price=price.RT, stat=stat.RT, time=time.RT, 
#                        tsteps=tsteps[1:n], paths=paths[1:n])
# save(res.RT.C, file="resRT_C.RData")