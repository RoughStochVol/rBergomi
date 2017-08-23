####################################################################
## Debug C++ code
####################################################################

####################################################################
## Check simulation of W and Wtilde
####################################################################

N <- 200 ## length of rows
H <- 0.07

Wtilde <- matrix(scan(file = "./Release/Wtilde.txt"), byrow=TRUE, ncol=N)
W1 <- matrix(scan(file = "./Release/W1.txt"), byrow=TRUE, ncol=N)
W1perp <- matrix(scan(file = "./Release/W1.txt"), byrow=TRUE, ncol=N)

## W1 is supposed to be an i.i.d. matrix of i.i.d. Gaussian vectors
## Check this!
shapiro.test(x = as.vector(W1)[1:5000]) ## p-value ~ 0.5
c.W1 <- cov(W1)
sum(sum(abs(c.W1 - diag(N))))/N^2

## Check Wtilde for normality
shapiro.test(as.vector(Wtilde)[1:5000]) ## p-value ~ 0.95
## let us try again for the 10th elements of each row -> should be normal then
shapiro.test(Wtilde[1:5000,10]) ## p-value ~ 0.5, the individual entries seems normal

## Try to test for multi-variate normality
mshapiro.test(t(Wtilde[1:5000,]))
mshapiro.test(t(W1[1:5000,]))
## This does not give normality, but not even for W1!

## Check the covariance matrix of Wtilde WRONG CALCULATIONS!!!
# c.Wtilde <- cov(Wtilde)
# source("generateVolterraCB.R")
# c.Volt <- vcv(N, H)
# sum(sum(abs(c.Wtilde - c.Volt)))/N^2 # Gives 0.0059, which seems reasonable. 

## Check the joint covariance matrix
c.joint <- cov(cbind(W1/sqrt(N), Wtilde))
c.VB <- VoltCov(N, H, 1)
sum(sum(abs(c.joint - c.VB)))/(2*N)^2 ## This is far off!
## For example:
c.joint[100,300]
c.VB[100,300]

## Recalculate Wtilde using R implementation of scheme
## Adjust to actually return Wtilde (= sqrt(2*H)*Y below)
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
        ##return(res)
        return(sqrt.2H*Y*expiry^H)
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

params <- c(0.04, H, 2, -0.9)
MC <- hybridScheme.vec(params)
steps <- N
expiries <- c(0.5)
M <- 100000
Wperp <- matrix(0, ncol=200, nrow=100000)

Wtilde.alt <- MC(M, steps, expiries, W1, W1perp, Wperp)

## Check the conditions on Wtilde.alt
c.Wtilde.alt <- cov(t(Wtilde.alt))
source("generateVolterraCB.R")
D <- diag(1:N/N*expiries)^H
c.Volt <- D%*%vcv(N, H)%*%D
sum(sum(abs(c.Wtilde.alt - c.Volt)))/N^2 # Gives 0.0059, which seems reasonable.
## Visual inspection; Both are far away!
c.Volt[1:5,1:5]
c.Wtilde[1:5,1:5]
c.Wtilde.alt[1:5,1:5]

c.Volt[195:200,195:200]
c.Wtilde[195:200,195:200]
c.Wtilde.alt[195:200,195:200]

## Check manually if c.Volt is correct
auto.cov <- function(s, t, H){
  f <- function(u) 2*H*abs(t-u)^(H-1/2)*abs(s-u)^(H-1/2)
  integrate(f, lower=0, upper=min(s,t))$value
}
s <- 75/N*expiries
t <- 175/N*expiries
auto.cov(s, t, H)
c.Volt[75,175]
c.Wtilde[75,175]
c.Wtilde.alt[75,175]
## they are not equal, but c.Wtilde seems closer
## Compute the full covariance matrix manually
## Note: Wtilde is for T=1!?
c.man <- outer(1:N/N*expiries, 1:N/N*expiries, Vectorize(function(s,t) auto.cov(s,t,H)))

## Comparison
c.Volt[1:5,1:5]
c.man[1:5,1:5]
c.Wtilde[1:5,1:5]
c.Wtilde.alt[1:5,1:5]
