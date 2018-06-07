## Test for J. Nierste: T = 1, xi = 0,235^2, eta = 1,9, rho = -0,9 und H = 0,07
source("rbergomi_console.R")

rBergomi.pricer2(0.235^2, 0.07, 1.9, -0.9, 1, 1, 500, 200000, 1)

rBergomi.pricer2(0.235^2, 0.07, 1.9, 0.0, 1, 1, 500, 200000, 1)

rBergomi.pricer2(0.235^2, 0.07, 1.9, -0.9, 1, 1, 500, 200000, 1)

rBergomi.pricer2(0.235^2, 0.43, 1.9, -0.9, 1, 1, 500, 200000, 1)

rBergomi.pricer2(0.235^2, 0.43, 1.9, -0.9, 1, 1, 2, 200000, 1)

rBergomi.pricer2(0.235^2, 0.43, 1.9, -0.9, 1, 0.8, 500, 200000, 1)

rBergomi.pricer2(0.235^2, 0.07, 1.9, -0.9, 1, 1, 2000, 400000, 1, TRUE)
