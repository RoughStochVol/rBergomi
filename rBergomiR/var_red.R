################################################################
## Check variance reduction
################################################################

source("rbergomi_console.R")

xi <- 0.04
H <- 0.07
eta <- 2.5
expiry <- 1
rho <- -0.9
K <- 1
N <- 200
M <- 10000
num.jobs <- 1

var.red <- FALSE
system.time(res1 <<- rBergomi.pricer2(xi, H, eta, rho, expiry, K, N, M, num.jobs, var.red))
Sys.sleep(3)
var.red <- TRUE
system.time(res2 <<- rBergomi.pricer2(xi, H, eta, rho, expiry, K, N, M, num.jobs, var.red))

res1
res2

######################################################################
## Check again weak rate 1
######################################################################

N.vec <- 2^(1:10)
M.vec <- 16*N.vec^2
res <- sapply(1:length(N.vec), function(i) rBergomi.pricer2(xi, H, eta, rho, 
                                                            expiry, K, N.vec[i], M.vec[i], 
                                                            num.jobs, FALSE))

price <- unlist(res[7,])
stat <- unlist(res[9,])

err <- abs(price[1:(length(price)-1)] - price[length(price)])
mc.err <- stat[1:(length(stat)-1)]
N <- N.vec[1:(length(N.vec)-1)]

lin <- lm(log(err[5:9]) ~ log(N[5:9]))

plot(N, err, type="l", log="xy")
lines(N, err+mc.err, lty=3)
lines(N, err-mc.err, lty=3)
lines(N, exp(-2)/N, lty = 2)

## Still looks ok

######################################################################
## Check Richardson extrapolation
######################################################################

## Take the result from the previous computation as reference value
ref.price <- price[10]
N.vec.rich <- 2^(1:7)
M.vec.rich <- 500*N.vec.rich^2
res.rich <- sapply(1:length(N.vec.rich), function(i) rBergomi.richardson(xi, H, eta, rho, 
                                        expiry, K, N.vec.rich[i], M.vec.rich[i], 
                                        num.jobs, FALSE))

price.rich <- unlist(res.rich[7,])
stat.rich <- unlist(res.rich[9,])

err.rich <- abs(price.rich - ref.price)

## Check Richardson directly from the previous results

ind.rich <- 1:(length(price)-1)
N.vec.rich <- 2^(ind.rich)
price.rich <- 2*price[ind.rich] - price[-1]
stat.rich <- sqrt(4 * stat[ind.rich]^2 + stat[-1]^2)
err.rich <- abs(price.rich - ref.price)

lin <- lm(log(err.rich) ~ log(N.vec.rich))

plot(N.vec.rich, err.rich, type="l", log="xy")
lines(N.vec.rich, err.rich+stat.rich, lty=3)
lines(N.vec.rich, err.rich-stat.rich, lty=3)
lines(N.vec.rich, exp(-2)/N.vec.rich, lty = 2)

## Something is wrong: The error without Richardson is typically smaller than the error
## with Richardson
err[ind.rich]/err.rich

## Results are worse?? What is happening?