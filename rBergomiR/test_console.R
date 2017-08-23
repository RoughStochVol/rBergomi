#####################################################################
## Test console program
#####################################################################

source("rbergomi_console.R")

## Generate parameters
n.par <- 100
H <- seq(0.01, 0.2, length.out=n.par)
eta <- seq(0.1, 4.0, length.out=n.par)
rho <- seq(-0.4, -0.999, length.out=n.par)

n.op <- 100
T <- rep(1.0, n.op)
K <- seq(0.8, 1.2, length.out=n.op)

xi <- 0.04

## further parameters
N <- 500
M <- 500000
num.threads <- 8
output.name <- "out.txt"

## run code
res <- rBergomi.pricer(xi, H, eta, rho, T, K, N, M, num.threads, output.name)

save(res, file="test.RData")

## Computes 10000 prices with 500 timesteps and 500000 samples in 1.4 hours with 8 cores.
