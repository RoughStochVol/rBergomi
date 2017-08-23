## Auxiliary stuff needed for testing the package
library(RBergomi)

#source("/home/bayerc/Documents/gatheral/fBM/r-code/EC2code/generateVolterraCB.R")

#############################################################
#
# Method 1: Generate the full covariance matrix.
#
#############################################################

# As before, we assume a uniform grid on [0,1]

# Generate Volterra covariance matrix
gExact <- function(h,x){
  require(hypergeo)
  a <- 4*h*x^(1/2 - h)/(1 + 2*h)
  hyp <- hypergeo(1, 1/2 - h, 3/2 + h, x)
  return(Re(a*hyp*(abs(hyp)>1e-8)))
}
vcv = function(n,h){
  uVect <- (1:n)/n
  g1 <- function(i,j){i^(2*h)*gExact(h,(i/j)^(j>=i))}
  m1 <- outer(uVect,uVect,g1)
  m1u <- upper.tri(m1,diag=F)*m1
  m1l <- t(m1u)
  res <- m1u+m1l
  diag(res) <- uVect^(2*h)
  return(res)
}

# Generate Brownian covariance matrix
# (Uses some regular expressions tricks.)
bm = function(n){
  tt = (1:n)/n
  vapply(seq_along(tt),
         function(i) c(tt[seq_len(i)], rep(tt[i], length(tt)-i)),
         numeric(length(tt)))
}

# generate cross-terms E[W^H * B]
cross = function(n, h, rho){
  a = rho*sqrt(2*h)/(h+1/2) # common factor
  # u = t, v=s
  f = function(u,v) u^(h+0.5) - (u - min(u,v))^(h+0.5)
  t.vec = (1:n)/n
  mat = outer(t.vec, t.vec, Vectorize(f))
  return(a*mat)
}
# generate cross-terms E[B * W^H]
cross2 = function(n, h, rho){
  a = rho*sqrt(2*h)/(h+1/2) # common factor
  # u = t, v=s
  f = function(v,u) u^(h+0.5) - (u - min(u,v))^(h+0.5)
  t.vec = (1:n)/n
  mat = outer(t.vec, t.vec, Vectorize(f))
  return(a*mat)
}

# the full matrix
VoltCov = function(n, h, rho){
  mat = matrix(NA, nrow=2*n, ncol=2*n)
  mat[1:n,1:n] = vcv(n,h)
  mat[(n+1):(2*n),(n+1):(2*n)] = bm(n)
  mat[1:n,(n+1):(2*n)] = cross(n, h, rho)
  mat[(n+1):(2*n),1:n] = t(mat[1:n,(n+1):(2*n)])
  return(mat)
}

# now generate the Volterra process and the Brownian motion
# m copies of the processes on an interval from 0 to tFinal
# are generated.
generate.VBM = function(n, h, rho, tFinal, m){
  BM = matrix(0, nrow=n+1, ncol=m)
  V = matrix(0, nrow=n+1, ncol=m)
  # generate the normalized covariance matrix
  S = VoltCov(n, h, rho)
  C = t(chol(S))
  # generate normal increments
  Z = matrix(rnorm(2*n*m), nrow=2*n, ncol=m)
  # apply cholesky decomposition
  Z = apply(Z, 2, function(x) C %*% x)
  V[2:(n+1),] = Z[1:n,]
  BM[2:(n+1),] = Z[(n+1):(2*n),]
  # scale according to time-interval
  V = tFinal^h * V
  BM = sqrt(tFinal) * BM
  return(list(Volt = V, BM = BM))
}

gen.rbergomi <- function(param, tFinal, tsteps, paths){
  ybm <- generate.VBM(n=tsteps, h=param$H, rho=param$rho, tFinal=1, m=paths)
  z <- ybm$BM
  wHatH <- ybm$Volt
  
  dtau <- tFinal/tsteps
  tau <- seq(0,tFinal,dtau)
  
  volterraScale <- tFinal^param$H
  tau2H <- tau^(2*param$H)
  vtnorm <- exp(param$eta*volterraScale*wHatH-param$eta^2*tau2H/2)
  
  xitu0 <- param$xitu(tau)
  
  vt <- vtnorm
  for (row in 1:(tsteps+1)){
    vt[row,] <- vtnorm[row,]*xitu0[row]  
  }
  sigt <- sqrt(vt)
  
  n1 <- tsteps+1
  totvars <- apply(sigt[-n1,]^2,2,cumsum)*dtau # Sum sig^2 for each column
  
  # Now we need to compute int sig dz
  dz <- apply(z,2,diff)*sqrt(tsteps) #dz have mean 0 and var=1
  sig.dz <- sigt[1:tsteps,]*dz*sqrt(dtau) # sigt needs to be adapted
  intSig.dz <- apply(sig.dz,2,cumsum)
  st <- exp(intSig.dz - 1/2*totvars)
  
  return(list(s=st, v=vt, totvar=totvars))
}

hybridScheme.vec <- function (params, use.fftw = FALSE){ # Assumes S0=1
  require(stats)
  is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
  use.fftw <- is.installed("fftwConvolve") && use.fftw
  my.convolve <- convolve
  if(use.fftw){
    require(fftwConvolve)
    my.convolve <- fftw.convolve
  }
  sigma <- sqrt(params[1])
  eta <- params[3]
  H <- params[2]
  alpha <- H - 1/2
  rho <- params[4]
  sqrt.2H <- sqrt(2 * H)
  rhoH <- sqrt.2H/(H + 1/2)
  covMatrix <- array(c(1, 1/(H + 1/2), 1/(H + 1/2), 1/(2 * H)), 
                     dim = c(2, 2))
  bstar <- function(k) {
    ((k^(alpha + 1) - (k - 1)^(alpha + 1))/(alpha + 1))^(1/alpha)
  }
  MC <- function(N, steps,expiries, seed=NULL) { # We assume kappa=1 and simulate with T=1 and effective parameters
    dt <- 1/steps
    dt1 <- dt^(H + 1/2)
    dt.2H <- dt^(2 * H)
    sqrt.dt <- sqrt(dt)
    dt.H <- dt^H
    x <- 1:steps
    Gamma <- bstar(x)^alpha
    Gamma[1] <- 0
    
    Simulation <- function(W1, W1hat, Z) {
      W1r <- W1[steps:1]
      Y2 <- my.convolve(Gamma, W1r, type = "open")[1:steps] # The convolution step is run once only
      S <- function(expiry){
        eta.t <- eta*expiry^H
        sigma.t <- sigma*sqrt(expiry)
        tseq <- eta.t * eta.t/2 * x^(2 * H) * dt.2H
        Y <- (W1hat + Y2) * dt.H
        v <- exp(eta.t * sqrt.2H * Y - tseq)
        v <- c(1, v[-steps])
        res <- exp(sigma.t*sum(v^0.5 * Z) * sqrt.dt - sigma.t^2/2 * sum(v) * dt) # Note that S0=1; Also cumsum replace by sum
        return(res)
      }
      return(S(expiries))
    }
    
    if(!is.null(seed))
      set.seed(seed)
    W1 <- rnorm(steps * N)
    W1perp <- rnorm(steps * N)
    Y1hat <- 1/sqrt.2H * (rhoH * W1 + sqrt(1 - rhoH^2) * W1perp)
    Wperp <- rnorm(steps * N)
    Z <- rho * W1 + sqrt(1 - rho * rho) * Wperp
    seq.path <- function(pathNum) {
      (1 + (pathNum - 1) * steps):(pathNum * steps)
    }
    f <- function(pathNum) {
      pick <- seq.path(pathNum)
      return(Simulation(W1[pick], Y1hat[pick], Z[pick]))
    }
    res <- sapply(1:N, f)
    return(res)
  }
  return(MC)
}

# Not vectorized!!
get.price <- function(s, K){
  payoff <- s - K
  payoff[payoff < 0] <- 0
  return(list(price=mean(payoff), stat=sd(payoff)/sqrt(length(s))))
}