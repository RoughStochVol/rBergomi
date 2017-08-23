
################################################################################
# Now a function to return one path of the Riemann-Liouville process
################################################################################

hybridScheme.Path <- function(params){
  
  S0 <- params$S0
  xi <- params$xi
  eta <- params$eta
  H <- params$H
  alpha <- H - .5
  rho <- params$rho
  
  sqrt.2H <- sqrt(2*H)
  rhoH <- sqrt.2H/(H+1/2) 
  
  covMatrix <- array(c(1,1/(H+1/2),1/(H+1/2),1/(2*H)),dim=c(2,2))
  
  bstar <- function(k){ # Function to find optimal integration point in discrete sum
    ((k^(alpha+1) - (k-1)^(alpha+1)) / (alpha+1)) ^ (1/alpha)
  }
  
  MC <- function(paths, steps, kappa, T=1, yFlag=F){ # If yFlag==T, returns the R-L process.
    
    N <- paths
    dt <- T/steps
    dt1 <- dt^(H+1/2)
    dt.2H <- dt^(2*H)
    sqrt.dt <- sqrt(dt)
    dt.H <- dt^H
    
    x <- 1:steps
    Gamma <- bstar(x)^alpha 
    Gamma[1] <- 0
    tseq <- eta*eta/2*x^(2*H)*dt.2H
    
    Simulation <- function(W1, W1hat, Z, yFlag=F){ # By default return the spots
      W1r <- W1[steps:1] # Reverse W1 by hand to save time
      Y2 <- convolve(Gamma, W1r, type="open")[1:steps]
      Y <- (W1hat+Y2)*dt.H ## The simulated R-L series
      v <- xi*exp(eta*sqrt.2H*Y - tseq)
      v <- c(xi, v[-steps])
      S <- S0 * exp(cumsum(v^0.5*Z)*sqrt.dt - 1/2*cumsum(v)*dt)
      if(yFlag) {return(Y)} else {return(S)}
    }
    
    W1 <- rnorm(steps*N)
    W1perp <- rnorm(steps*N)
    Y1hat <- 1/sqrt.2H*(rhoH * W1 + sqrt(1-rhoH^2) * W1perp)
    Wperp <- rnorm(steps*N)
    Z <- rho * W1 + sqrt(1-rho*rho)*Wperp
    
    seq.path <- function(pathNum){(1+(pathNum-1)*steps):(pathNum*steps)}
    
    f <- function(pathNum){ # This function generates one path
      pick <- seq.path(pathNum)
      return(Simulation(W1[pick],Y1hat[pick],Z[pick],yFlag=yFlag))
    }
    res <- sapply(1:N,f)  # Version with no parallelization
    
    return(t(res))
    
  }
  return(MC)
}

################################################################################
# Example
# 
# n <- 100000
# params <- list(S0=1, xi=0.235^2, eta=1.9, H=0.1, rho=-0.9)
# yPaths <- hybridScheme.Path(params)(paths=2, steps=n,kappa=1,T=1,yFlag=T)
# x <- (1:n)/n
# plot(x,yPaths[,1],type="l",col="red")
#
################################################################################

hybridScheme <- function (params){ # Assumes S0=1
  
  sigma <- sqrt(params$xi)
  eta <- params$eta
  H <- params$H
  alpha <- H - 1/2
  rho <- params$rho
  sqrt.2H <- sqrt(2 * H)
  rhoH <- sqrt.2H/(H + 1/2)
  covMatrix <- array(c(1, 1/(H + 1/2), 1/(H + 1/2), 1/(2 * 
                                                         H)), dim = c(2, 2))
  bstar <- function(k) {
    ((k^(alpha + 1) - (k - 1)^(alpha + 1))/(alpha + 1))^(1/alpha)
  }
  MC <- function(N, steps,expiries) { # We assume kappa=1 and simulate with T=1 and effective parameters
    
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


#####################################################################################################
## Hybrid scheme of Pakkanen et al. to simulate RBergomi
## Code amended by Christian Bayer

## Input:
## params ... vector (rather than a list) of parameters of the model (xi, H, eta, rho)
## use.fftw ... boolean controlling whether to use FFTW to compute the convolution
## Extra optional argument seed to use the same seed for each run
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

################################################################################
# Example
# 
# steps <- 200
# params <- list(xi=0.235^2, eta=1.9, H=0.1, rho=-0.9)
# res <- hybridScheme(params)(N=10000, steps=n,expiries = c(.25,.5,.75,1))
#
################################################################################


hybridSchemeH <- function (H){
  require(stats)
  alpha <- H - 1/2
  sqrt.2H <- sqrt(2 * H)
  rhoH <- sqrt.2H/(H + 1/2)
  bstar <- function(k) {
    ((k^(alpha + 1) - (k - 1)^(alpha + 1))/(alpha + 1))^(1/alpha)
  }
  function(N, steps) {
    dt <- 1/steps
    dt1 <- dt^(H + 1/2)
    dt.2H <- dt^(2 * H)
    sqrt.dt <- sqrt(dt)
    dt.H <- dt^H
    x <- 1:steps
    Gamma <- bstar(x)^alpha
    Gamma[1] <- 0
    Simulation <- function(W1, W1hat) {
      W1r <- W1[steps:1]
      Y2 <- convolve(Gamma, W1r, type = "open")[1:steps]
      Wtilde <- (W1hat + Y2) * dt.H
      return(Wtilde)
    }
    W1 <- rnorm(steps * N)
    W1perp <- rnorm(steps * N)
    W1hat <- 1/sqrt.2H * (rhoH * W1 + sqrt(1 - rhoH^2) * W1perp)
    seq.path <- function(pathNum) {
      (1 + (pathNum - 1) * steps):(pathNum * steps)
    }
    f <- function(pathNum) {
      pick <- seq.path(pathNum)
      return(Simulation(W1[pick], W1hat[pick]))
    }
    pathMatrix <- matrix(NA,nrow=steps*2,ncol=N)
    pathMatrix[1:steps,] <- matrix(W1,nrow=steps,ncol=N) # Paths of W
    pathMatrix[(1:steps)+steps,] <- sapply(1:N, f) # Paths of Wtilde
    return(pathMatrix) # Upper rows are W, lower rows are Wtilde
    
  }
}

################################################################################
# Recoded hybrid scheme
################################################################################

################################################################################
# Wtilde
################################################################################
Wtilde.sim <- function(W1, W1perp, H){
  
  library(stats)
  steps <- dim(W1)[1] #W1 and W1perp are steps x N matrices
  N <- dim(W1)[2]
  alpha <- H - 1/2
  stopifnot(dim(W1perp) == c(steps,N))
  
  bstar <- function(k) {
    ((k^(alpha + 1) - (k - 1)^(alpha + 1))/(alpha + 1))^(1/alpha)
  }
  
  dt <- 1/steps
  dt1 <- dt^(H + 1/2)
  dt.2H <- dt^(2 * H)
  sqrt.dt <- sqrt(dt)
  sqrt.2H <- sqrt(2 * H)
  dt.H <- dt^H
  rhoH <- sqrt.2H/(H + 1/2)
  x <- 1:steps
  
  Gamma <- bstar(x)^alpha
  Gamma[1] <- 0
  
  f <- function(j){ # j is the path number so 1 <= j <= N
    W1r <- W1[steps:1,j]
    Y2 <- convolve(Gamma, W1r, type = "open")[1:steps]
    W1hat <- 1/sqrt.2H * (rhoH * W1[,j] + sqrt(1 - rhoH^2) * W1perp[,j])
    return((W1hat + Y2) * dt.H)
  }
  
  Wtilde <- sqrt.2H * sapply(1:N, f)
  return(Wtilde)
}

################################################################################
# hybridScheme2
################################################################################
hybridScheme2 <- function(xi,params)function(N,steps,expiries) 
{
  eta <- params$eta
  H <- params$H
  rho <- params$rho
  
  W1 <- matrix(rnorm(N*steps),nrow = steps,ncol=N)
  W1perp <- matrix(rnorm(N*steps),nrow = steps,ncol=N)
  
  #     W1 <- W1 - apply(W1,1,mean)
  #     W1 <- W1/apply(W1,1,sd)
  
  #     W1perp <- W1perp - apply(W1perp,1,mean)
  #     W1perp <- W1perp/apply(W1perp,1,sd)
  
  Wperp <- matrix(rnorm(N*steps),nrow = steps,ncol=N)
  #     Wperp <- Wperp - apply(Wperp,1,mean)
  #     Wperp <- Wperp/apply(Wperp,1,sd)
  
  Z <- rho * W1 + sqrt(1 - rho * rho) * Wperp
  
  Wtilde <- Wtilde.sim(W1, W1perp, H)
  
  S <- function(expiry) {
    
    dt <- expiry/steps
    ti <- (1:steps)*dt # t-grid
    
    Wtilde.H <- expiry^H * Wtilde # rescale for the time-interval
    xi.t <- xi(ti) # xi evaluated on the grid
    v1 <- xi.t * exp(eta * Wtilde.H - 1/2 * eta^2 * ti^(2*H))
    v0 <- rep(xi(0),N)
    v <- rbind(v0,v1[-steps,])  
    logs <- apply(sqrt(v*dt) * Z - v/2*dt,2,sum)
    s <- exp(logs)
    #         s <- s/mean(s) # Normalize
    return(s)
  }
  
  st <- t(sapply(expiries,S))
  
  #return(list(W=W1,Wtilde=Wtilde,st=st))
  return(st)
}

################################################################################
# hybridScheme for cXM, cMM and zeta
################################################################################
zeta.MC <- function(xi,params)function(N,steps,expiries) 
{
  eta <- params$eta
  H <- params$H
  rho <- params$rho
  
  W1 <- matrix(rnorm(N*steps),nrow = steps,ncol=N)
  W1perp <- matrix(rnorm(N*steps),nrow = steps,ncol=N)
  
  Wtilde <- Wtilde.sim(W1, W1perp, H)
  
  compute.zeta <- function(expiry) {
    
    dt <- expiry/steps
    ti <- (1:steps)*dt # t-grid
    
    Wtilde.H <- expiry^H * Wtilde # rescale for the time-interval
    xi.t <- xi(ti) # xi evaluated on the grid
    v1 <- xi.t * exp(eta * Wtilde.H - 1/2 * eta^2 * ti^(2*H))
    v0 <- rep(xi(0),N)
    v <- rbind(v0,v1[-steps,])  
    int.sqrtv.dW <- apply(sqrt(v*dt) * W1,2,sum) # Vector of length N
    int.v.dt <- apply(v*dt,2,sum) # Vector of length N
    cXM <- rho * mean(int.sqrtv.dW * int.v.dt) 
    cMM <- var(int.v.dt)
    zeta <- -cXM + cMM/4
    return(list(cXM=cXM,cMM=cMM,zeta=zeta))
  }
  
  zeta.res <- t(sapply(expiries,compute.zeta))
  
  return(zeta.res)
}