#setwd("/Users/JGatheral/Google Drive/BayerFrizPaper/EC2code")

# Generate the Volterra process together with a correlated Bm
# In the sense of the manuscript.

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

# ybm <- generate.VBM(n=1000, h=0.1, rho=-.85, tFinal=1, m=100000) # This generates w1 and wHat1 in the notation of rSabrSmile

# More elegantly: creat a R6 class for generating samples of the
# volterra process
library(R6)
GenVolt = R6Class("GenVolt",
  public = list(
    initialize = function(n, h, rho, tFinal){
      private$n = n
      private$tFinal = tFinal
      private$h = h
      private$rho = rho
      private$C = t(chol(VoltCov(n, h, rho)))
    },
    generate = function(m = 1){
      BM = matrix(0, nrow=private$n+1, ncol=m)
      V = matrix(0, nrow=private$n+1, ncol=m)
      Z = matrix(rnorm(2*private$n*m), nrow=2*private$n, ncol=m)
      if(m == 1)
        Z = private$C %*% Z
      else
        Z = apply(Z, 2, function(x) private$C %*% x)
      V[2:(private$n+1),] = Z[1:private$n,]
      BM[2:(private$n+1),] = Z[(private$n+1):(2*private$n),]
      V = private$tFinal^private$h * V
      BM = sqrt(private$tFinal) * BM
      if(m == 1){ # return vectors then
        V = as.vector(V) # may be slow!!!!
        BM = as.vector(BM)
      }
      return(list(Volt = V, BM = BM))
    }
    ),
  private = list(
    n = NA, # length of trajectory
    tFinal = NA,
    h = NA,
    rho = NA,
    C = NA # Cholesky matrix
    )
)

#############################################################
#
# Method 2: Use Euler type discretization
#
#############################################################

# More precisely, we use Decreusefond and Ustunel
generate.VBM.DU = function(n, h, rho, tFinal, m){
  t.vec = tFinal*(0:n)/n
  dt.vec = diff(t.vec)
  sdt.vec = sqrt(dt.vec)
  tt.vec = (tFinal - t.vec)^(h+1/2)
  a = sqrt(2*h)/((h+1/2)*1/n)
  # function mapping BM increments dw (of length n) to Volt-path
  f = function(dw) cumsum((tt.vec[1:n]-tt.vec[2:(n+1)])*dw)
  # generate increments of W
  dw = matrix(rnorm(n*m), nrow = n, ncol = m)
  dw = apply(dw, 2, function(x) sdt.vec*x)
  # generate Volterra process
  V = matrix(0, nrow=n+1, ncol=m)
  V[2:(n+1),] = apply(dw, 2, f)
  # generate the Brownian motion
  BM = matrix(0, nrow=n+1, ncol=m)
  BM[2:(n+1),] = rho*apply(dw, 2, cumsum)
  dw = matrix(rnorm(n*m), nrow = n, ncol = m)
  dw = apply(dw, 2, function(x) sdt.vec*x)
  BM[2:(n+1),] = BM[2:(n+1),] + sqrt(1-rho^2)*apply(dw, 2, cumsum)
  return(list(Volt=V, BM=BM))
}
# The above version is incorrect! Try the one described in Jim's write-up...
# (See RSVVsimulation.pdf)
generate.VBM.DU.Jim = function(n, h, rho, tFinal, m){
  t.vec = tFinal*(0:n)/n
  dt = tFinal/n
  sdt = sqrt(dt)
  # first generate dW
  dw = sdt*matrix(rnorm(n*m), nrow = n, ncol = m)
  # for a given column dwc compute the vector W^H
  f = function(i, t.vec, h, dt, dwc){
    tCur = t.vec[i+1]
    integrand = sqrt((tCur - t.vec[1:i])^(2*h) - (tCur - t.vec[2:(i+1)])^(2*h))
    sum(integrand*dwc[1:i])/dt
  }
  f2 = function(dwc, t.vec, h, dt)
    sapply(1:n, function(i) f(i, t.vec, h, dt, dwc))
  # now apply over columns of dw
  V = matrix(0, nrow=n+1, ncol=m)
  V[2:(n+1),] = sapply(1:m, function(j) f2(dw[,j], t.vec, h, dt))
    # generate the Brownian motion
  BM = matrix(0, nrow=n+1, ncol=m)
  BM[2:(n+1),] = rho*apply(dw, 2, cumsum)
  dw = sdt*matrix(rnorm(n*m), nrow = n, ncol = m)
  BM[2:(n+1),] = BM[2:(n+1),] + sqrt(1-rho^2)*apply(dw, 2, cumsum)
  return(list(Volt=V, BM=BM))
}

