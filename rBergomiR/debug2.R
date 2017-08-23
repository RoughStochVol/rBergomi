#######################################################
## The Riemann-Liouville fBm
#######################################################
## Systematic debugging.

## Fix the parameters 
N <- 200 ## length of rows = number of steps
H <- 0.07
expiry <- 0.5 # maturity in years = time point of last step

###########################################################
## Manual computations
###########################################################

## First compute the true covariance matrix manually and check against
## the covariance matrix as implemented 2 years ago by Jim and myself.
auto.cov <- function(s, t, H){
  f <- function(u) 2*H*abs(t-u)^(H-1/2)*abs(s-u)^(H-1/2)
  integrate(f, lower=0, upper=min(s,t), subdivisions=1000)$value
}
## Check against Mathematica: for s = 0.2625, t = 0.6125, 
## we get cov = 0.163774
s <- (75/N)*expiry
t <- (175/N)*expiry
auto.cov(s, t, H)
c.man <- outer((1:N/N)*expiry, (1:N/N)*expiry, Vectorize(function(s,t) auto.cov(s,t,H)))
## Note that the old manually computed covariance matrix was wrong!!

## Use the old code by Jim and myself
source("generateVolterraCB.R")
c.Volt <- vcv(N, H) ## Note that this is for expiry = 1, hence we need to rescale!
c.Volt <- expiry^(2*H) * c.Volt

## Compare both
sum(sum(abs(c.man-c.Volt))) ## 0.00003 ~ 0: they coincide.

#########################################################
## Simulations in the C++ code
#########################################################
## Note that Wtilde should not need any re-scalings
W1 <- matrix(scan(file = "./Release/W1.txt"), byrow=TRUE, ncol=N)
W1perp <- matrix(scan(file = "./Release/W1perp.txt"), byrow=TRUE, ncol=N)
Wperp <- matrix(scan(file = "./Release/Wperp.txt"), byrow=TRUE, ncol=N)
Wtilde <- expiry^H * matrix(scan(file = "./Release/Wtilde.txt"), byrow=TRUE, ncol=N)
#Wtilde.noscaled <- matrix(scan(file = "./Release/Wtilde.txt"), byrow=TRUE, ncol=N)
#Wtilde.scaled <- matrix(scan(file = "./Release/WtildeScaled.txt"), byrow=TRUE, ncol=N)
v <- matrix(scan(file = "./Release/v.txt"), byrow=TRUE, ncol=N)
S <- scan(file = "./Release/S.txt")
c.Wtilde <- cov(Wtilde)

## compare
(c.Wtilde[1:5,1:5] - c.man[1:5,1:5])/c.man[1:5,1:5]
(c.Wtilde[195:200,195:200] - c.man[195:200,195:200])/c.man[195:200,195:200]
## OLD: c.Wtilde is very far off on the diagonal, otherwise not that bad.
## NEW: The values are usually off by about 10^{-3}, which could
## possibly be explained by the Monte Carlo error.

## Motivated by that observation, we look at the correlation matrices.
rho.man <- cov2cor(c.man)
rho.Wtilde <- cov2cor(c.Wtilde)

rho.man[1:5,1:5]
rho.Wtilde[1:5,1:5]
## Correlation matrix is off by a factor 2 approximately (off-diagonally)

#########################################################
## R implementation of BLP-algorithm
########################################################

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
        return(res)
        ##return(sqrt.2H*Y*expiry^H)
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

params <- c(0.04, H, 2, -0.9) ## arbitrary values for eta, xi_0 and rho
MC <- hybridScheme.vec(params)
steps <- N
M <- 100000
Wtilde.alt <- MC(M, N, expiry, W1, W1perp, Wperp)
c.Wtilde.alt <- cov(t(Wtilde.alt))
c.man[1:5,1:5]
c.Wtilde.alt[1:5,1:5]
c.Wtilde[1:5,1:5]

## Conclusion: all three are different!!

############################################################
## Alternative R implementation of BLP algorithm
############################################################
hybrid.scheme <- function(W1, W1perp, H){
  library(stats)
  N <- dim(W1)[1]
  steps <- dim(W1)[2]
  stopifnot(dim(W1perp) == c(N, steps))
  alpha <- H - 1/2
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
  f <- function(n){
    Y2 <- convolve(Gamma, W1[n,steps:1], type = "open")[1:steps]
    W1hat <- 1/sqrt.2H * (rhoH * W1[n,] + sqrt(1 - rhoH^2) * W1perp[n,])
    ## W1hat = linearComb(rhoH/s2H, W1, sqrt(1.0 - rhoH*rhoH)/s2H, W1perp); 
    ## Seems equal
    return((W1hat + Y2) * dt.H)
  }
  Wtilde <- sqrt.2H * t(sapply(1:N, f)) ## NOTE: This is for T = 1! (Just like the C++ code)
  ## Wtilde = linearComb(sqrt(2.0*H)*pow(1.0/N, H), Y2, sqrt(2.0*H)*pow(1.0/N, H), W1hat)
  return(Wtilde)
}
Wtilde.3 <- expiry^H * hybrid.scheme(W1, W1perp, H)
## Note that t(Wtilde.3) is exactly equal to Wtilde.alt

## The first goal is to have both the R- and the C++ version produce
## The same values...
test.convolve <- function(H, steps){
  library(stats)
  alpha <- H - 1/2
  bstar <- function(k) {
    ((k^(alpha + 1) - (k - 1)^(alpha + 1))/(alpha + 1))^(1/alpha)
  }
  W1 <- 0:(steps-1)
  x <- 1:steps
  Gamma <- bstar(x)^alpha
  Gamma[1] <- 0
  W1r <- W1[steps:1]
  Y2 <- convolve(Gamma, W1r, type = "open")[1:steps]
  return(list(Gamma=Gamma, W1=W1, Y2=Y2))
}
test <- test.convolve(H, N)
test$Gamma ## the same as in the corresponding C++ test case
test$W1 ## the same as in the corresponding C++ test case
test$Y2 ## the same as in the C++ test case (up to very small differences)
## This resulted in a bug found in the C++ code.

###########################################################
## R and C++ implementation
###########################################################
## Both R implementations give the same paths when the input
## is equal. This is not true for the C++ implementation.
## The Convolution seems to give the same results pathwise,
## though.
Wtilde[1:5,1:5]
Wtilde.3[1:5,1:5]

## Run the hybridscheme for W1 = 0:(N-1) and W1perp = W1^2
x <- matrix(0:(N-1), nrow=2, ncol=N, byrow=TRUE)
hybrid.scheme(x,x^2, H)[1,]
## This is exactly the same as the output of the C++ code for the same input!
x <- W1[1:2,]
y <- W1perp[1:2,]
expiry^H * hybrid.scheme(x,y,H)[1,1:5]
Wtilde.3[1,1:5]
## Output (first 5 entries) copied from C++:
expiry^H * c(0.0670043, 1.06309, -1.54809, -0.752725, 0.634369)
## This is essentially the same!!! Why is the W1 or W1perp that we get not the one that is used internally??
## Running the C++ code again, it seems that different W1perp is used... (W1 and Wtilde seem equal to before)
## UPDATE: Now it works, the problem was that instead of W1perp.txt W1.txt was read twice!

###################################################################
## Double check the option prices using the R implementation of BLP
###################################################################
## Old implementation used reference price 0.2149403
K <- 0.8
M <- 100000
params <- c(0.07, H, 2.2, -0.9) 
MC <- hybridScheme.vec(params)
S.mat <- MC(M, N, expiry, W1, W1perp, Wperp)
S <- S.mat
price.R <- sum(S[S > K] - K)/M # = 0.2153002 => A bit different, but probably in statistical error range.
## This is very different from the "true price" computed earlier!!!
## Double check by generating new random numbers
MC <- hybridScheme.vec(params)
S.mat <- MC(M, N, expiry, matrix(rnorm(M*N), nrow=M), 
            matrix(rnorm(M*N), nrow=M), matrix(rnorm(M*N), nrow=M))
S <- S.mat#[,N]
price.R <- sum(S[S > K] - K)/M # = 0.2149029 => The same as before!
## Gives the same result! (Which is very different from other R implementations...)

## Note the refernece price was computed using the Cholesky pricer...

#########################################################################
## Debug information for hybrid scheme.
#########################################################################
## output either Wtilde, v, or S; Note that the returned values cannot be seen as
## trajectories of the processes Wtilde etc.: only the marginal distributions should
## fit, but not the joint ones.
## IS THIS THE CURRENT PROBLEM????
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
          res <- Wtilde[steps]
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

## Given the trajectories of W1, W1perp, compute the values of the scaled Wtilde at time = expiry
MC <- hybridScheme.vec.debug(params, output = "Wtilde")
Wtilde.alt <- MC(M, N, expiry, W1, W1perp, Wperp)
head(Wtilde[,N])/head(Wtilde.alt)
## They coincide almost exactly -> The above R code produces the same W-tilde as the C++-code

#############################################################
## Compare scaled Wtilde
#############################################################
## This time use the version of Wtilde already scaled in C++
Wtilde[1:5,1:5]/Wtilde.scaled[1:5,1:5]
Wtilde.scaled[1,1:5]/Wtilde.noscaled[1,1:5]
Wtilde.scaled[2,1:5]/Wtilde.noscaled[2,1:5]
## WtildeScaled seems to work fine.

#############################################################
## Compare v.
#############################################################
MC <- hybridScheme.vec.debug(params, output = "v")
v.alt <- MC(M, N, expiry, W1, W1perp, Wperp)
head(v[,N])/head(v.alt)
## These values are almost equal! (Relative differences of the order of 1e-3. Why?)

##############################################################
## Compare S
###############################################################
MC <- hybridScheme.vec.debug(params, output = "S")
S.alt <- MC(M, N, expiry, W1, W1perp, Wperp)
head(S)/head(S.alt)
## Seems fine: Relative differences of the order 3e-3

#############################################################
## Check the prices
#############################################################
price.R
mean((S.alt-K) * (S.alt > K))
mean((S-K) * (S > K))
## prices seem ok!