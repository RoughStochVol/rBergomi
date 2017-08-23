################################################################
## Compute the objective function for the calibration problem
################################################################

source("BlackScholes.R")
source("rbergomi_console.R")

## Extract all strikes and expiries in terms of moneyness from ivolData,
## a dataframe of option price data
## Strikes and expiries are only considered when there are both bid and ask
## prices available.
## Also keep the prices in this dataframe.
## ATM vols are computed by interpolation
iv.surface <- function(ivolData){
  library(stinepack)
  texp <- ivolData$Texp
  expDates <- unique(texp)
  f <- ivolData$Fwd
  K <- ivolData$Strike/f
  ## now check if prices are available
  bid.available <- !is.na(ivolData$Bid)
  ask.available <- !is.na(ivolData$Ask)
  both.available <- bid.available & ask.available
  midVols <- 0.5 * (ivolData$Bid + ivolData$Ask)
  atmVols <- numeric(length(texp))
  ## interpolate midVols to get ATMVol
  for(t in expDates){
    include <- (!is.na(midVols)) & (texp == t)
    volInterp <- function(xout) {
      stinterp(x = log(K[include]), y = midVols[include], xout)$y
    }
    atmVols[texp == t] <- volInterp(0)
  }
  ## gather expiries and strikes in dataframe
  df <- data.frame(Texp = texp[both.available], K = K[both.available], 
                   Bid = ivolData$Bid[both.available], Ask = ivolData$Ask[both.available],
                   Mid = midVols[both.available], ATM = atmVols[both.available])
  ## order according to expiries
  df <- df[with(df, order(Texp)), ]
  return(df)
}

## soft.equal is an equality function taking into account rounding errors
## due to C++ IO.
soft.equal <- function(x, y, eps = 1e-07){
  return(abs(x - y) < eps)
}

## objectiveMC.new compoutes the objective function at some parameter
## specification. It uses the auxiliary functions objectiveMC.plotLayout
## and objectiveMC.plotBasic.
## Input:
## ivolData ... dataframe of option price date
## plot ... boolean controlling whether the model IVs should be plotted.
## params ... dataframe of parameter values. If the length is > 1 and plot == TRUE,
##            then the first one is printed. 
## Note: slices are not implemented, as the user is expected to slice the input ("ivolData")
## instead.


objectiveMC <- function (ivolData, params, tsteps, paths) 
{
  ## Isolate H, eta, rho from params
  params.new <- data.frame(H = params$H, eta = params$eta, rho = params$rho)
  library(parallel) ## used to detect the number of cores
  ## first clean the ivolData
  ivolData.clean <- iv.surface(ivolData)
  num.threads <- detectCores()
  output.name <- ".out.txt"
  ## compute all the option prices in the model
  modelData <- rBergomi.pricer(params$xi, params$H, params$eta, params$rho, 
                               ivolData.clean$Texp, ivolData.clean$K, tsteps, paths, 
                               num.threads, output.name)

  ## Normalize both market (mid-) vols and model vols by the observed ATM vols.
  market.vols.norm <- ivolData.clean$Mid / ivolData.clean$ATM
  ## Which row of ivolData.clean corresponds to time t and moneyness k
  translate <- function(t, k){
    temp <- (t == ivolData.clean$Texp) & (k == ivolData.clean$K)
    return(which.max(temp))
  }
  f <- function(i) modelData$iv[i] / ivolData.clean$ATM[translate(modelData$T[i], 
                                                                  modelData$K[i])]
  model.vols.norm <- sapply(1:length(modelData$iv), f)
  ## compute the normalized spreads
  spread.norm <- (ivolData.clean$Ask - ivolData.clean$Bid) / ivolData.clean$Mid
  
  ## compute the loss for each row param in the matrix of parameter values params.new
  ## Here the order of columns in params.new is important, as we cannot address them by
  ## name, only by index!
  loss <- function(param){
    ## isolate all entries in model.vols.norm corresponding to the chosen parameters
    index <- soft.equal(modelData$H, param[1]) & soft.equal(modelData$eta, param[2]) &
      soft.equal(modelData$rho, param[3])
    ## sanity check: sum(index) = length(ivolData.clean$Texp)
    stopifnot(sum(index) == length(ivolData.clean$Texp))
    ## check that the order was not destroyed in the C++-code
    stopifnot(sum(abs(modelData$T[index] == ivolData.clean$Texp)) < 0.001)
    stopifnot(sum(abs(modelData$K[index] == ivolData.clean$K)) < 0.001)
    ## check that the model-vols are not too far away from the market ATM vols
    admit <- (0.1 < model.vols.norm[index]) & (model.vols.norm[index] < 10)
    ## Add sum of squares of errors, weighted by number of admitted guys
    dist <- (model.vols.norm[index] - market.vols.norm) / spread.norm * admit
    sq.dist <- sum(dist^2)
    weight <- sum(admit)
    return(sq.dist / weight)
  }
  obj <- apply(params.new, 1, loss)
  
  return(data.frame(params, obj = obj))
}


## In case that plotting is also required...

## heuristic function for the layout of nSlices plots arranged in one
## frame.
objectiveMC.plotLayout <- function(nSlices){
  colnum <- sqrt(nSlices * 2)
  rows <- round(colnum/2, 0)
  columns <- round(colnum, 0)
  while (rows * columns < nSlices) {
    rows <- rows + 1
  }
  return(list(columns=columns, rows=rows))
}

## This is a lazy hack -> should rather get data as proper arguments!!!
objectiveMC.plotBasic <- function(env = parent.frame()){
  bidVoln <- env$bidVol[env$include]/env$atmVol ## only used for plotting
  askVoln <- env$askVol[env$include]/env$atmVol ## only used for plotting
  ybottom <- 0.8 * min(bidVoln) ## moved inside if for better readability
  ytop <- 1.2 * max(askVoln)
  xrange <- c(min(env$z), max(env$z))
  yrange <- c(ybottom, ytop)
  plot(env$z, bidVoln, col = "red", pch = 24, cex = 0.5, xlim = xrange, 
       ylim = yrange, 
       main = paste("T =", format(t, digits = 2, nsmall = 2)), 
       xlab = "Log-Strike", ylab = "Implied Vol.")
  points(env$z, askVoln, col = "blue", pch = 25, cex = 0.5)
  abline(h = 1, lty = 2)
  abline(v = 0, lty = 2)
}


## Compute Chebyshev interpolation of the objective function
## num ... vector containing the number of Chebyshev polynomials in each dimension
## domain ... domain of the parameters in each dimension
## xi ... xi parameter of rBergomi, not calibrated
## ivolData ... a dataframe contianing the option price data
## paths ... number of paths for MC
## steps ... time discretization
## Return value: a function taking three arguments (H, eta, rho) and returning the value
## of the objective function.
objectiveMC.cheb <- function(num, domain, xi, ivolData, steps, paths){
  library(chebpol)
  ## First generate the grid-points
  knots <- chebknots(c(H = num[1], eta = num[2], rho = num[3]), intervals=domain)
  knots.total <- expand.grid(knots)
  ## re-order so that H varies slowest and rho fastest
  ord <- with(knots.total, order(H, eta))
  knots.total.ordered <- knots.total[ord, ]
  ## add xi
  params <- cbind(knots.total.ordered, xi)
  ## Compute the objective function for all the parameter values
  ## in knots.total
  obj <- objectiveMC(ivolData, params, steps, paths)
  ## Compute Chebyshev approximation
  re.ord <- order(ord)
  obj.re.order <- obj$obj[re.ord]
  dim(obj.re.order) <- num
  cheb <- chebappx(obj.re.order, intervals = domain)
  return(list(knots=knots, cheb=cheb))
}

## Interactive plotting of the objective function
## The objective is plotted as function in either H and rho, 
## with eta fixed, but adjustable on a slider, or in H and eta.
## The only inputs are the Chebyshev polynomial and
## the plot range for the values.
## CAUTION: ONLY WORKS IN RSTUDIO.
## Input:
## obj.cheb ... Chebyshev interpolation of objective function,
##              as returned by objectiveMC.cheb
## domain ... domain of the Chebyshev interpolation, a list of intervals
## val.lim ... plotting range in the z (= value) axis
## params ... string encoding the choice of parameters to be plotted. 
##            Implemented are "H-rho" and "H-eta"
plot.objective <- function(obj.cheb, domain, val.lim = c(0,10), 
                           params = "H-rho"){
  require(manipulate)
  H.int <- domain[[1]]
  eta.int <- domain[[2]]
  rho.int <- domain[[3]]
  if(params == "H-rho"){
    g <- function(eta){
      H.grid <- seq(H.int[1], H.int[2], length.out=100)
      rho.grid <- seq(rho.int[1], rho.int[2], length.out=100)
      f <- function(H, rho) obj.cheb(c(H, eta, rho))
      z <- outer(H.grid, rho.grid, Vectorize(f))
      persp(H.grid, rho.grid, z, xlab="H", ylab="rho", 
            phi=20, theta=-25, ticktype="detailed",
            zlim=val.lim)
    }
    manipulate(g(eta), eta=slider(eta.int[1], eta.int[2]))  
  }
  else if(params == "H-eta"){
    g <- function(rho){
      H.grid <- seq(H.int[1], H.int[2], length.out=100)
      eta.grid <- seq(eta.int[1], eta.int[2], length.out=100)
      f <- function(H, eta) obj.cheb(c(H, eta, rho))
      z <- outer(H.grid, eta.grid, Vectorize(f))
      persp(H.grid, eta.grid, z, xlab="H", ylab="eta", 
            phi=20, theta=-25, ticktype="detailed",
            zlim=val.lim)
    }
    manipulate(g(rho), rho=slider(rho.int[1], rho.int[2]))  
  }
}

## Optimize the parameters
## Input:
## obj.cheb .. Chebyshev interpolation of objective function
## domain ... domain of interpolation
## guess.method ... "sobol" for choosing initial guesses according to
##                  Sobol points, "chebyshev" for using the Chebyshev knots.
## K ... number of Sobol points used as intial guesses. Otherwise,
##       the number of polynomials in each coordinate (a vector).
## method ... method for optim
## Output: the full output of optim at the minimizing initail guess, 
## together with the minimizing initial guess and the var of the
## minimal values for all the initial guesses.
## Remark: obj.cheb is a polynomial, so we should be able to take the
## derivatives explicitly to use in the optimization. This does not seem
## possible in chebpol.
## A dedicated multivariate polynomial optimizer does not seem to exist in R.

## First an auxiliary function to construct the initial guesses
optim.objective.param0 <- function(guess.method, domain, K){
  if(guess.method == "sobol"){ ## use Sobol points
    stopifnot(length(K)==1) # stop if K is not a number
    param0 <- sobol(K, dim=3) ## these are in [0,1]
    to.domain <- function(x, domain){
      y <- rep(NA, length(domain))
      for(i in 1:length(domain))
        y[i] <- domain[[i]][1] + x[i]*(domain[[i]][2] - domain[[i]][1])
      return(y)
    }
    param0 <- t(apply(param0, 1, function(x) to.domain(x, domain)))
  }
  else{ ## use Chebyshev grid points
    stopifnot(length(K) == 3)
    knots <- chebknots(K, intervals=domain)
    param0 <- expand.grid(knots)
  }
  colnames(param0) <- c("H", "eta", "rho")
  return(param0)
}

## Now the main routine
optim.objective <- function(cheb.fun, domain, K, 
                            guess.method="sobol", 
                            method="L-BFGS-B"){
  library(randtoolbox)
  library(chebpol)
  lower <- sapply(domain, function(x) x[1]) # lower bound in form assumed by optim
  upper <- sapply(domain, function(x) x[2])
  ## Choose initial guesses
  param0 <- optim.objective.param0(guess.method, domain, K)
  n <- dim(param0)[1] ## number of initial guesses
  ## The optimization routine
  f <- function(param0) optim(param0, cheb.fun, method=method, 
                              lower=lower, upper=upper)
  ## apply for each choice of intial parameter
  res <- apply(param0, 1, f)
  ## extract the values, find the arg.min and return the corresponding values
  res.value <- sapply(1:n, function(i) res[[i]]$value)
  arg.min.index <- which.min(res.value)
  ret <- list(par=res[[arg.min.index]]$par, 
              value=res[[arg.min.index]]$value,
              counts=res[[arg.min.index]]$counts,
              convergence=res[[arg.min.index]]$convergence,
              message=res[[arg.min.index]]$message,
              param0=param0[arg.min.index,],
              var.values=var(res.value))
  return(ret)
}