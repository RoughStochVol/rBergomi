##############################################################
## Kernel desnity estimation and random number generation.
##############################################################

## Goal: Take data of maturities and strikes for SPX options,
## smooth them by a smoothing kernel and finally simulate from
## the corresponding density.

##############################################################
## Generate sample data
##############################################################

## Use mixture of normals

library(mvtnorm)
n <- 1000
mus <- rbind(c(1,2), c(0.4,-1), c(0,0))
sigmas <- rbind(0.1*diag(2), matrix(c(0.2,0.21,0.21,0.3), nrow = 2), 0.2*diag(2))
props <- c(0.2, 0.6, 0.2)
x <- rmvnorm.mixt(n = n, mus = mus, Sigmas = sigmas, props = props)
plot(x[,1], x[,2]) ## Scatterplot
x1.grid <- seq(-1, 1.5, length.out = 100)
x2.grid <- seq(-2.5, 3, length.out = 100)
f <- function(x)
  dmvnorm.mixt(x, mus, sigmas, props)
f.vec <- Vectorize(function(x,y) f(c(x,y)))
z <- outer(x1.grid, x2.grid,  f.vec)
contour(x1.grid, x2.grid, z)

## Candidate 1: ks
library(ks)
## Use two different "selectors" for bandwidth
H1 <- Hpi(x)
H2 <- Hpi.diag(x)
## Estimate the density
fhat1 <- kde(x, H1)
fhat2 <- kde(x, H2)
## Plots
plot(fhat1)
plot(fhat2)
## Sample from the estimated density
y <- rkde(n = 1000, fhat = fhat1)
plot(y[,1], y[,2]) ## Scatterplot

## Check the actual data from Ben. Date is 20180215
spx.data <- read.csv(file = 
    "/Users/bayerc/Documents/ben/deep-learning-finance/spx_20180215.csv", 
    header = TRUE)
## Calls and Puts contains both maturity and strike, which need to be extracted!