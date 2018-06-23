################################################################
## Haar basis
################################################################
## Testing the C++ code.

h_1_0 <- function(t){
  x <- 0 * t
  x[(t >= 0) & (t <= 0.5)] <- 1
  x[(t <= 1) & (t >= 0.5)] <- -1
  return(x)
}

h <- function(t,n,k){
  if(n == 0)
    x <- rep(1, length(t))
  else
    x <- 2^((n-1)/2) * h_1_0((2^(n-1)) * t - k)
  return(x)
}

t.grid <- seq(0, 1, length.out = 100)

plot(t.grid, h(t.grid, 0, 0), type = "l", xlab = "t", ylab = "Haar", col = "black", 
     ylim = c(-2, 2))
lines(t.grid, h(t.grid, 1, 0), lty = 1, col = "red")
lines(t.grid, h(t.grid, 2, 0), lty = 1, col = "blue")
lines(t.grid, h(t.grid, 2, 1), lty = 1, col = "green4")
lines(t.grid, h(t.grid, 3, 0), lty = 2, col = "black")
lines(t.grid, h(t.grid, 3, 2), lty = 2, col = "black")

h.scaled <- function(t, n, k, tMax) h(t/tMax, n, k) / sqrt(tMax)

plot(t.grid, h.scaled(t.grid, 0, 0, 0.8), type = "l", xlab = "t", ylab = "Haar", 
     col = "black", ylim = c(-2, 2))
lines(t.grid, h.scaled(t.grid, 1, 0, 0.8), lty = 1, col = "red")
lines(t.grid, h.scaled(t.grid, 2, 0, 0.8), lty = 1, col = "blue")
lines(t.grid, h.scaled(t.grid, 2, 1, 0.8), lty = 1, col = "green4")
lines(t.grid, h.scaled(t.grid, 3, 0, 0.8), lty = 2, col = "black")
lines(t.grid, h.scaled(t.grid, 3, 2, 0.8), lty = 2, col = "black")

integrate(function(t) h.scaled(t, 2, 1, 0.8)^2, lower = 0, upper = 0.8)

## Now define the Schauder basis (this is done to check the formula)!
schauder <- function(t, n, k, tMax){
  f <- function(s) h.scaled(s, n, k, tMax)
  return(integrate(f, lower = 0, upper = t)$value)
}

schauder.vec <- Vectorize(schauder, vectorize.args = "t")

schauder(0.3, 3, 1, 0.8)
schauder.vec(c(0.3, 0.6), 3, 1, 0.8)

## Generate the numbers needed for testing the C++ code.
tGrid <- seq(0, 0.8, by = 0.1)
lVec <- c(0, 1, 2)
out.val <- c()
for(l in lVec){
  for(j in seq(0, ceiling(2^(l-1))-1, by = 1)){
    out.val <- c(out.val, schauder.vec(tGrid, l, j, 0.8))
  }
}

options(scipen = 999)
out.val
options(scipen = 0)

## The full integrated basis function
tMax <- 0.8
K <- 3
dt <- tMax / K
iota2alpha <- function(iota){
  alpha <- rep(NA, 4)
  alpha[3] <- iota %% 2
  iota1 <- (iota - alpha[3]) %/% 2
  alpha[4] <- K - (iota1 %% K) - 1
  iota2 <- (iota1 - (K - alpha[4] - 1)) %/% K
  if(iota2 == 0){
    alpha[1] <- 0
    alpha[2] <- 0
  }
  else{
    alpha[1] <- floor(log2(iota2)) + 1
    alpha[2] <- iota2 - 2^(alpha[1] - 1)
  }
  return(alpha)
}

Phi <- function(t, iota, i){
  alpha <- iota2alpha(iota)
  ret <- rep(NA, length(t))
  if(alpha[3] != i){
    ret <- rep(0, length(t))
  } else{
    k <- alpha[4]
    ret[t < k*dt] <- 0
    ret[t > (k+1) * dt] <- schauder.vec(dt, alpha[1], alpha[2], dt)
    ret[((k*dt) <= t) & (t <= (k+1)*dt)] <- 
      schauder.vec(t[((k*dt) <= t) & (t <= (k+1)*dt)] - k*dt, alpha[1], alpha[2], dt)
  }
  return(ret)
}

plot.Phi <- function(iota){
  t <- seq(0, tMax, length.out = 500)
  par(mfrow = c(2,1))
  plot(t, Phi(t, iota, 0), type="l", ylab = "Phi")
  plot(t, Phi(t, iota, 1), type="l", ylab = "Phi")
}

for(i in 0:23)
  plot.Phi(i)

## Evaluate each of these guys on a grid of length 12
tGrid <- seq(0, tMax, length.out = 13)
f <- function(iota) c(Phi(tGrid, iota, 0), Phi(tGrid, iota, 1))
vals <- as.vector(sapply(0:23, f))
write(vals, file = "Phivals.txt", ncolumns = length(vals), sep = ", ")

## Generate Brownian motion
dB.generate <- function(M, i){
  dB1 <- matrix(0, nrow = M, ncol = length(tGrid)-1)
  for(iota in 0:(2*length(tGrid)-1)){
    dPhi <- Phi(tGrid[-1], iota, i) - Phi(tGrid[1:(length(tGrid)-1)], iota, i)
    dB1 <- dB1 + outer(rnorm(M), dPhi)
  }
  return(dB1)
}

M <- 1000000
dB1 <- dB.generate(M, 0)
cov1 <- cov(dB1)
cov1

## The same problem as in the C++ code.

## Now compute the theoretical covariance.
dB.covariance <- function(i){
  cov.dB <- matrix(0, nrow = length(tGrid)-1, ncol = length(tGrid)-1)
  for(iota in 0:(2*length(tGrid)-1)){
    dPhi <- Phi(tGrid[-1], iota, i) - Phi(tGrid[1:(length(tGrid)-1)], iota, i)
    cov.dB <- cov.dB + outer(dPhi, dPhi)
  }
  return(cov.dB)
}

cov1.infty <- dB.covariance(0)
cov1.infty

max(abs(cov1 - cov1.infty))
## Theory and practise agree. But I do not understand the theory?
## One problem is the correlation between dB[1] and dB[5]
## I have no idea why I need unlist below -> should be a vector!
sapply(0:23, function(iota) unlist(Phi(tGrid[2], iota, 0)) - 
         unlist(Phi(tGrid[1], iota, 0)))
sapply(0:23, function(iota) unlist(Phi(tGrid[6], iota, 0)) - 
         unlist(Phi(tGrid[5], iota, 0)))

###########################################################################
## Haar representation of the fBm
###########################################################################

## The shifted and scaled Haar functions, i.e. the derivative of Phi
phi <- function(t, iota, i){
  alpha <- iota2alpha(iota)
  ret <- rep(NA, length(t))
  if(alpha[3] != i){
    ret <- rep(0, length(t))
  } else{
    k <- alpha[4]
    ret[t < k*dt] <- 0
    ret[t > (k+1) * dt] <- 0
    ret[((k*dt) <= t) & (t <= (k+1)*dt)] <- 
      h.scaled(t[((k*dt) <= t) & (t <= (k+1)*dt)] - k*dt, alpha[1], alpha[2], dt)
  }
  return(ret)
}

plot.phi <- function(iota){
  t <- seq(0, tMax, length.out = 500)
  par(mfrow = c(2,1))
  plot(t, phi(t, iota, 0), type="l", ylab = "phi")
  plot(t, phi(t, iota, 1), type="l", ylab = "phi")
}

for(i in 0:23)
  plot.phi(i)

## The fractional kernel
ker <- function(t, s, H) {
  ret <- rep(0, length(s))
  ret[s - t < -1e-12] <- sqrt(2*H) * (t - s)^(H - 0.5)
  return(ret);
}

## Integrate the kernel against the Haar function
Phi.tilde <- function(t, H, iota){
  f <- function(tt) {
    g <- function(s) ker(tt, s, H) * phi(s, iota, 0)
    return(integrate(g, lower = 0, upper = tt, subdivisions = 5000L,
                     rel.tol = .Machine$double.eps^0.5/1000,
                     abs.tol = .Machine$double.eps^0.5/1000,
                     stop.on.error = FALSE)$value)
  }
  ret <- sapply(t, f)
  return(ret)
}

## Try it out
Phi.tilde(0.5, 0.07, 2)
Phi.tilde(c(0.2,0.5), 0.07, 2)

## Generate values of Phi.tilde for comparison with C++-code
## Evaluate each of these guys on a grid of length 12
tGrid <- seq(0, tMax, length.out = 13)
H <- 0.07
f <- function(iota) Phi.tilde(tGrid, H, iota)
vals <- as.vector(sapply(0:23, f))
write(vals, file = "Phitildevals.txt", ncolumns = length(vals), sep = ", ")

## These integrals seem to be inaccurate compared to the onses provided
## by Mathematica. Hence, we will now switch to the Mathemtica implementation
## to provide reference values.

########################################################################
## rBergomi for H = 1/2
########################################################################
## Pedestrian, independent implementation of rBergomi
## for the special case H = 0.5

Bergomi <- function(xi, eta, rho, expiry, strike, N, M){
  dt <- expiry/N
  rhobar <- sqrt(1 - rho^2)
  ## dW1[1,] is the first trajectory of the increment process
  dW1 <- matrix(sqrt(dt) * rnorm(N*M), ncol=N)
  dW2 <- matrix(sqrt(dt) * rnorm(N*M), ncol=N)
  t.grid <- seq(0, expiry, length.out = N+1)
  W1 <- cbind(0,t(apply(dW1, 1, cumsum)))
  ## v; the strange t(...t(...)) form is needed for broadcasting
  v <- xi * exp(t(eta * t(W1) - 0.5 * eta^2 * t.grid))
  ## Euler scheme for X = log(S)
  X <- rep(0, M)
  for(n in 1:N)
    X <- X - 0.5 * v[,n] * dt + sqrt(v[,n]) * (rho * dW1[,n] + rhobar * dW2[,n])
  S <- exp(X)
  payoff <- (S - strike) * (S > strike)
  return(list(price = mean(payoff), stat = sd(payoff)/sqrt(M)))
}

xi <- 0.04
eta <- 1.9
rho <- -0.9
expiry <- 0.8
strike <- 0.95
N <- 1000
M <- 100000

ret <- Bergomi(xi, eta, rho, expiry, strike, N, M)