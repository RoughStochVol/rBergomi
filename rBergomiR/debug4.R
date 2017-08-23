###########################################################################
## Debug quadrature version of C++ code.
###########################################################################
hybridScheme.vec.debug <- function (params, output = "S", use.fftw = FALSE, max.cores = NULL){ 
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
                (dim(Wperp) == c(N,steps)) && (length(expiries) == 1))
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
      qoi <- function(expiry){
        eta.t <- eta*expiry^H
        sigma.t <- sigma*sqrt(expiry)
        tseq <- 0.5 * eta^2 * expiry^(2*H) * x^(2 * H) * dt.2H
        Wtilde <- expiry^H * sqrt.2H * (W1hat + Y2) * dt.H
        v <- exp(eta * Wtilde - tseq)
        v <- c(1, v[-steps])
        if(output == "S")
          res <- exp(sigma.t*sum(v^0.5 * Z) * sqrt.dt - sigma.t^2/2 * sum(v) * dt)
        else if(output == "Wtilde")
          res <- Wtilde
        else if(output == "v")
          res <- params[1] * v[steps]
        return(res)
      }
      return(qoi(expiries))
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

params <- c(0.04, 0.1, 2.4, -0.85)
M <- 1
N <- 16
expiry <- 1

W1 <- c(-0.541023, -0.100824, -0.224605, -0.114091, 1.93134, -0.370674, -1.46428, 
        0.512231, 0.955507, -0.430132, -2.43461, -0.0773463, 0.030667, -1.89986, 
        -0.40995, 1.02311)
W1perp <- c(-0.297753, 0.391882, -0.516504, 1.20351, -0.76676, 0.478933, 0.804645, 
            -0.434525, 1.01745, -0.514849, 1.33975, 1.29001, -0.845606, 1.66612, 
            -1.68535, -0.939346) 
Wperp <- rep(NA, length(W1))
Wtilde.C <- c(-0.456047, -0.016567, -0.544785, 0.343098, 0.496035, 0.410568, -0.242436, 
              -0.204666, 0.991739, -0.284237, -0.654118, -0.0586418, -0.967448, 
              -0.705564, -2.06443, -0.858508) 
MC <- hybridScheme.vec.debug(params, output = "Wtilde")
Wtilde.R <- MC(M+1, N, expiry, rbind(W1,NA), rbind(W1perp,NA), rbind(Wperp,NA))[,1]
 