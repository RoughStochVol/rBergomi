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
