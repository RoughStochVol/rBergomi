#############################################################################
# This requires the following source files:
# hybridScheme.R
# BlackScholes.R


#############################################################################
# The following function outputs computes the weighted squared distance between
# empirical implied volatilities (from ivolData) and the implied volatilities 
# corresponding to the Monte Carlo paths in mcMatrix
objectiveMC <- function (ivolData, mcMatrix, slices = NULL, plot = F) 
{
  require(stinepack)
  gam <- 1/2
  bidVols <- as.numeric(ivolData$Bid)
  askVols <- as.numeric(ivolData$Ask)
  expDates <- unique(ivolData$Texp)
  nSlices <- length(expDates)
  if (!is.null(slices)) {
    nSlices <- length(slices)
  }
  else {
    slices <- 1:nSlices
  }
  if (plot) {
    layout <- objectiveMC.plotLayout(nSlices)
    par(mfrow = c(layout$rows, layout$columns), mex = 0.5)
  }
  atmVols <- numeric(nSlices)
  atmSkews <- numeric(nSlices)
  sq.dist <- 0
  l.dist <- 0
  for (slice in slices) {
    t <- expDates[slice]
    texp <- ivolData$Texp
    bidVol <- bidVols[texp == t]
    askVol <- askVols[texp == t]
    midVol <- (bidVol + askVol)/2
    f <- (ivolData$Fwd[texp == t])[1]
    k <- log(ivolData$Strike[texp == t]/f)
    include <- !is.na(midVol)
    kIn <- k[include]
    volInterp <- function(xout) {
      stinterp(x = kIn, y = midVol[include], xout)$y
    }
    atmVol <- volInterp(0)
    midVoln <- midVol[include]/atmVol
    spreadn <- askVol[include]/atmVol - bidVol[include]/atmVol
    z <- k[include]/(atmVol * t^gam)
    volInterpn <- function(xout) {
      stinterp(x = z, y = midVoln, xout)$y
    }
    atmVol <- volInterp(0)
    atmSkewn <- (volInterpn(0.01) - volInterpn(-0.01))/0.02
    atmVols[slice] <- atmVol
    atmSkews[slice] <- atmSkewn
    if (plot) 
      objectiveMC.plotBasic()
    spots <- mcMatrix[slice, ]
    s0 <- mean(spots)
    volMC <- function(k) {
      bsOut(spots, t, s0 * exp(k))$BSV
    }
    atmVolMC <- volMC(0)
    kMC <- kIn * atmVolMC/atmVol
    volsMCn <- volMC(kMC)/atmVolMC
    admit <- (volsMCn > 0.1) & (volsMCn < 10)
    if (plot) 
      points(z, volsMCn, col = "orange", pch = 20, cex = 0.5)
    dist <- (volsMCn - midVoln)/spreadn * admit
    sq.dist <- sq.dist + sum(dist^2)
    l.dist <- l.dist + sum(admit)
  }
  obj <- sq.dist/l.dist
  return(list(expiries = expDates, atmVols = atmVols, atmSkews = atmSkews, 
              obj = obj))
}


#############################################################################
# The following function outputs a matrix of MC paths
# Steps in rows; paths in columns

gen.paths <- function(params, expiries, paths, steps, distrib = 100, 
                      max.cores = 4, seed = NULL){
  require(foreach)
  require(doParallel)
  iters<- max(1,ceiling(paths/distrib))
  cl.num <- min(detectCores(), max.cores) # number of cores
  cl<-makeCluster(cl.num)
  registerDoParallel(cl)
  ## run simulation in parallel, output everything in mcMatrix
  f <- hybridScheme.vec(params)
  mcMatrix <- foreach(icount(iters),.combine='cbind') %dopar% {
    f(N=distrib, steps=steps, expiries=expiries, seed = seed)
  }
  stopCluster(cl)
  return(mcMatrix)
}

#############################################################################
# Given H, the following function outputs a matrix of MC paths of W and Wtilde
# Steps in rows; paths in columns

gen.WtildeH <- function(H, paths, steps, distrib = 100, max.cores = 4){
  require(foreach)
  require(doParallel)
  iters<- max(1,ceiling(paths/distrib))
  cl.num <- min(detectCores(), max.cores) # number of cores
  cl<-makeCluster(cl.num)
  registerDoParallel(cl)
  ## run simulation in parallel, output everything in mcMatrix
  f <- hybridSchemeH(H)
  mcMatrix <- foreach(icount(iters),.combine='cbind') %dopar% {
    f(N=distrib, steps=steps)
  }
  stopCluster(cl)
  return(mcMatrix)
}


#############################################################################
## Compute the objective function given model parameters
#
## Input parameters
## params ... vector containing xi, H, eta, rho
## ivolData ... data frame containing market data
## paths ... number of paths for MC
## steps ... time discretization
## slices ... vector of indices to expiries, or NULL 
## distrib ... number of iterations for foreach 
## (each iteration computing paths/contrib trajectories)
## max.cores ... maximal number of cores
gen.objectiveMC <- function(params, ivolData, paths, steps, slices = NULL, 
                            distrib = 100, max.cores = 4, seed = NULL){
  expiries <- unique(sort(ivolData$Texp))
  if(!is.null(slices))
    expiries <- expiries[slices]
  mcMatrix <- gen.paths(params, expiries, paths, steps, distrib, 
                        max.cores, seed)
  return(objectiveMC(ivolData, mcMatrix, slices)$obj)
}

#############################################################################
## Compute the objective function given model parameters
#
## Input parameters
## params ... vector containing xi, H, eta, rho
## ivolData ... data frame containing market data
## paths ... number of paths for MC
## steps ... time discretization
## slices ... vector of indices to expiries, or NULL 
## distrib ... number of iterations for foreach 
## (each iteration computing paths/contrib trajectories)
## max.cores ... maximal number of cores
gen.objectiveMCtilde <- function(paramsTilde, ivolData, paths, steps, slices = NULL, 
                            distrib = 100, max.cores = 4, seed = NULL){
  expiries <- unique(sort(ivolData$Texp))
  if(!is.null(slices))
    expiries <- expiries[slices]
  params <- paramsTilde
  params[3] <- paramsTilde[3]/sqrt(2*paramsTilde[2]) # Convert \tilde \eta to \eta
  mcMatrix <- gen.paths(params, expiries, paths, steps, distrib, 
                        max.cores, seed)
  return(objectiveMC(ivolData, mcMatrix, slices)$obj)
}

# Test the code:
# setwd("/Users/JGatheral/Documents/rFiles")
# source("hybridScheme.R")
# source("BlackScholes.R")
# load("/Users/JGatheral/Google Drive/SPX Data/spxOptionMetricsIVols.rdata")
# ivolData <- spxIvolList[["20100204"]]
# params.test <- c(0.04,0.05,1.82,-.94)
# system.time(obj.test <- gen.objectiveMC(params.test, ivolData, paths=10000, steps=200, slices = NULL, 
#                                         distrib = 100, max.cores = 8, seed=123))
# obj.test