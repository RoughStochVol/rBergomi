#############################################################
## Check numerical performance of rBergomi implementation
#############################################################

## Parameters:
xi <- 0.04
H <- 0.06
eta <- 2.5
rho <- -0.8
tFinal <- 0.8
K <- 1.0

## First the standard scheme.

load(file="res_C.RData")

n <- length(res.C$price)


## take the last entry as reference value and do regression of the error to get
## the convergence rate
## Note that results suggest a weak rate of 1 instead of H+1/2...
## Note that the statistical error is too large!
##
## Note that we really need to take a much more accurate discretization as reference value!

t.res <- res.C$price[n]
error <- abs(res.C$price[-n] - t.res)
rel.error <- error/res.C$price[n]
lin <- lm(log(error[1:6]) ~ log(res.C$tsteps[1:6]))
## Even when taking only very few terms, the weak rate 1 is confirmed
## plot against the theoretically expected one
plot(res.C$tsteps[-n], error, log="xy", type="l", xlab="steps", ylab="error")
lines(res.C$tsteps[-n], 0.05*res.C$tsteps[-n]^(-1),#exp(lin$coefficients[1]) * res.C$steps[-n]^(-1),
      lty = 2)
lines(res.C$tsteps[-n], error + res.C$stat[-n], lty=3)
lines(res.C$tsteps[-n], error - res.C$stat[-n], lty=3)

## The same including the Romano-Touzi trick
load(file="resRT_C.RData")

n <- length(res.RT.C$price)
error <- abs(res.RT.C$price - t.res)

rel.error <- error/t.res
lin <- lm(log(error) ~ log(res.RT.C$tsteps))
## plot against the theoretically expected one
pdf(file="weak_RT.pdf", height=4, width=6)
par(mar=c(5, 4, 2, 2) + 0.1)
plot(res.RT.C$tsteps, error, log="xy", type="l", xlab="steps", ylab="error")
lines(res.RT.C$tsteps, 0.05*res.RT.C$tsteps^(-1),#exp(lin$coefficients[1]) * res.RT.C$steps[-n]^(-1),
      lty = 2)
lines(res.RT.C$tsteps, error + res.RT.C$sta, lty=3)
lines(res.RT.C$tsteps, error - res.RT.C$stat, lty=3)
dev.off()
## Again, we observe rate 1 of weak convergence and the statistical error is still too large

## Compare standard scheme with Romano-Touzi-scheme
n <- length(res.RT.C$stat)
mean(res.RT.C$stat / res.C$stat[1:n])
## So the statistical error is on average 30% smaller for the Romano-Touzi-trick
## Computational time:
n <- length(res.RT.C$stat)
mean(res.RT.C$time / res.C$time[1:n])
## The computational time required is about 10% less for the RT-version. Interestingly
## with incredibly high varition????
