####################################################################
## Debug C++ code for Premia
####################################################################

source("rbergomi_console.R")

## Check the Premia results:
xi <- 0.07
H <- 0.07
eta <- 2.2
rho <- -0.9
expiry <- 1.0
K <- 1.0
N <- 200
M <- 500000

## With these parameters, the Premia code gives the following results:
## price = 0.084 with hybrid scheme,
## price = 0.1264 with Cholesky.
## New version in Premia gives: 0.0860079 with Cholesky

###################################################################
## C++ implementation of hybrid scheme
###################################################################

output.name <- "out.txt"
time.cpp <- system.time(price.cpp 
                        <<- rBergomi.pricer(xi, H, eta, rho, expiry, K, 
                                            N, M, 1, output.name))
output.name <- "out.txt"
time.cpp.mult <- system.time(price.cpp.mult 
                        <<- rBergomi.pricer(xi, H, eta, rho, expiry, 
                                            K, N, M, 8, output.name))

## gives 0.08529019 as with Premia (hybrid), or 0.08518164 with M = 500000, time: 24.3
## Note: multi-threaded version gives 0.08532558, time 4.5. (But there is probably a bug.)

###################################################################
## R implementation of hybrid scheme
###################################################################

## Note that the Romano-Touzi trick is not used here!

hybridScheme.vec <- function (params, use.fftw = FALSE, max.cores = NULL){ 
  require(stats)
  require(parallel)
  require(foreach)
  require(doParallel)
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
  ## We assume kappa=1 and simulate with T=1 and effective parameters
  MC <- function(N, steps, expiries, W1, W1perp, Wperp) {
    ## W1, W1perp, Wperp are now assumed to be Nxsteps matrices
    stopifnot((dim(W1) == c(N,steps)) && 
                (dim(W1perp) == c(N,steps)) && 
                (dim(Wperp) == c(N,steps)))
    dt <- 1/steps
    dt1 <- dt^(H + 1/2)
    dt.2H <- dt^(2 * H)
    sqrt.dt <- sqrt(dt)
    dt.H <- dt^H
    x <- 1:steps
    Gamma <- bstar(x)^alpha
    Gamma[1] <- 0
    
    ## Simulation takes vectors of length steps.
    Simulation <- function(W1, W1hat, Z) {
      W1r <- W1[steps:1]
      # Y2 <- my.convolve(Gamma, W1r, type = "open")[1:steps] # The convolution step is run once only
      Y2 <- convolve(Gamma, W1r, type = "open")[1:steps] # The convolution step is run once only
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
    Y1hat <- 1/sqrt.2H * (rhoH * W1 + sqrt(1 - rhoH^2) * W1perp)
    Z <- rho * W1 + sqrt(1 - rho * rho) * Wperp
    if(is.null(max.cores))
      res <- sapply(1:N, function(n) Simulation(W1[n,], Y1hat[n,], Z[n,]))
    else{
      ## parSapply does not seem to work well
      distrib <- 100
      cl.num <- min(detectCores(), max.cores)
      cl <- makeCluster(cl.num)
      #       clusterExport(cl, varlist = c("W1", "Y1hat", "Z", "Simulation", "N"),
      #                     envir = environment())
      #       res <- parSapply(cl, 1:N, function(n) Simulation(W1[n,], Y1hat[n,], Z[n,]))
      registerDoParallel(cl)
      N.cl.num <- ceiling(N/(cl.num*distrib))
      pick <- function(i) (1 + (i-1) * N.cl.num):min(i*N.cl.num, N)
      res <- foreach(i = 1:cl.num, .combine = 'cbind', 
                     .export = ls(envir = environment(MC))) %dopar% {
                       sapply(pick(i), function(n) Simulation(W1[n,], Y1hat[n,], Z[n,]))
                     }
      stopCluster(cl)
    }
    return(res)
  }
  return(MC)
}

params <- c(xi, H, eta, rho)
pricer.R <- function(params, expiry, N, M){
  MC <- hybridScheme.vec(params)
  W1 <- matrix(rnorm(M*N), nrow = M)
  W1perp <- matrix(rnorm(M*N), nrow = M)
  Wperp <- matrix(rnorm(M*N), nrow = M)
  S <- MC(M, N, expiry, W1, W1perp, Wperp)
  mean((S - K) * (S > K))
}

time.R <- system.time(price.R <<- pricer.R(params, expiry, N, M))

## This also provides 0.08507181, or 0.08421449 for M = 500000, time: 97.1

###################################################################
## R implementation of Cholesky scheme
###################################################################

source("generateVolterraCB.R")

param <- list(H = H, rho = rho, eta = eta, xitu = function(t) rep(xi, length(t)))

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

time.chol <- system.time(S <<- gen.rbergomi(param = param, 
                                            tFinal = expiry, tsteps = N, 
                                            paths = M)$s[N,])
price.chol.R <- mean((S - K) * (S > K))

## Gives 0.08396558, close to the value of the hybrid scheme.
# For M = 500000: price = 0.08441938, time 209