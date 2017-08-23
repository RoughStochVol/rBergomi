## test calibration_objective
source("calibration_objective.R")

load("spxOptionMetricsIVols.rdata")
ivolData <- spxIvolList[["20100204"]]

H <- c(0.05, 0.1)
eta <- c(1.5, 2.5)
rho <- c(-0.8, -0.9)
xi <- 0.04
params <- expand.grid(H=H, eta=eta, rho=rho, xi=xi)

obj <- objectiveMC(ivolData, params, 100, 100000)

## parameters for Chebyshev interpolation
num <- c(10, 10, 10)
domain <- list(c(0.05, 0.15), c(0.5, 3), c(-0.99, -0.7))
xi <- 0.04
paths <- 500000
steps <- 500
obj <- objectiveMC.cheb(num, domain, xi, ivolData, steps, paths)

K <- 100
ret <- optim.objective(obj$cheb, domain, K, guess.method="sobol", method="L-BFGS-B")

