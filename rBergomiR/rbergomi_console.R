###########################################################################
## Call the console version of the RBergomi program from within R
###########################################################################

## xi, H, eta, rho ... model parameters. Except for xi, these are assumed to be given
##                     as vectors.
## T, K ... option parameters.
## Prices are computed for each choice H[i], eta[i], rho[i] and each choice T[j], K[j].
## I.e., no new combinations of H[i] and eta[j], for instance, is considered.
## M ... number of samples to be used
## N ... number of time steps for the Euler scheme.
## num.threads ... number of threads to be used.
## output.name ... file name
## Prices and IVs are then computed by:
## 1) Save parameters in text files in the folder path contining the rBergomi executable.
## 2) The executable is called, computes prices and IVs and saves them in a text file
##    with name output.name in paths.
## 3) The results are read from the output file.
rBergomi.pricer <- function(xi, H, eta, rho, T, K, N, M, num.threads, output.name,
                            debugP = FALSE){
  ## force non-scientific notation of numbers
  #scipen.old <- options()$scipen
  #options(scipen = 999)
  ## Convert number to string containing first 10 digits
  f10 <- function(x) formatC(x, digits=10, format="f")
  
  ## Save in files expected by console
  path <- "/Users/bayerc/Documents/workspace/roughBergomi/rBergomi/Release/"
  if(debugP)
    path <- "/Users/bayerc/Documents/workspace/roughBergomi/rBergomi/Debug/"
  
  write.table(f10(H), file=paste(path,".H.txt", sep=""), col.names=FALSE, row.names=FALSE, 
              sep="\n", quote=FALSE)
  write.table(f10(eta), file=paste(path,".eta.txt", sep=""), col.names=FALSE, row.names=FALSE, 
              sep="\n", quote=FALSE)
  write.table(f10(rho), file=paste(path,".rho.txt", sep=""), col.names=FALSE, row.names=FALSE, 
              sep="\n", quote=FALSE)
  write.table(f10(T), file=paste(path,".T.txt", sep=""), col.names=FALSE, row.names=FALSE, 
              sep="\n", quote=FALSE)
  write.table(f10(K), file=paste(path,".K.txt", sep=""), col.names=FALSE, row.names=FALSE, 
              sep="\n", quote=FALSE)
  write.table(f10(xi), file=paste(path,".xi.txt", sep=""), col.names=FALSE, row.names=FALSE, 
              sep="\n", quote=FALSE)
  
  finf <- function(x) formatC(x, format="d")
  ## call console version of rBergomi pricer
  sys.command <- paste(paste(path, "rBergomi", sep=""), finf(N), finf(M), 
                       finf(num.threads), path, output.name)
  system(sys.command)
  
  ## restore old options
  #options(scipen=scipen.old)
  
  ## read the result in again
  res <- read.table(paste(path, output.name, sep=""), header=TRUE)
  
  return(res)
}

## rBergomi.recover re-reads the output-file from step 3 in rBergomi.pricer.
## Note that no new prices are computed!!
rBergomi.recover <- function(output.name){
  path <- "/Users/bayerc/Documents/workspace/roughBergomi/rBergomi/Release/"
  res <- read.table(paste(path, output.name, sep=""), header=TRUE)
  return(res)
}

## Split vector x in list of n vectors of size length(x)/n. x[i] is put in 
## z[[i %% n]]
split.vec <- function(x, n){
  split(x, seq_along(x) %% n)
}

## Single threaded version. Here, multiple processes are run in parallel
## using sime " & " syntax.
## xi, H, eta, rho ... model parameters. All assumed to be vectors.
## T, K ... option parameters.
## Prices are computed for each choice xi[i], H[i], eta[i], rho[i], T[i], K[i].
## M ... number of samples to be used
## N ... number of time steps for the Euler scheme.
## num.jobs ... number of jobs to be used.
## var.red ... boolean: should variance reduction be used?
## Prices and IVs are then computed by:
## 1) Save parameters in text files in the folder path contining the rBergomi executable.
## 2) The executable is called, computes prices and IVs and saves them in a text file
##    with name output.name in paths.
## 3) The results are read from the output file.
##
## BUG!
## There seems to be a bug when the number of jobs is higher than previously, i.e.,
## when one of the output files is generated for the first time.
## The bug happens in line 145.
rBergomi.pricer2 <- function(xi, H, eta, rho, T, K, N, M, num.jobs, var.red = FALSE){
  ## force non-scientific notation of numbers
  #scipen.old <- options()$scipen
  #options(scipen = 999)
  ## Convert number to string containing first 10 digits
  f10 <- function(x) formatC(x, digits=10, format="f")
  
  ## Save in files expected by console
  if(Sys.info()[1] == "Darwin")
    path <- "~/Documents/workspace/roughBergomi/ML_console_ST/"
  else
    path <- "~/workspace/roughBergomi/ML_console_ST/"
  
  ## Seperate all the input vectors in num.jobs different vectors
  xi.list <- split.vec(xi, num.jobs)
  H.list <- split.vec(H, num.jobs)
  eta.list <- split.vec(eta, num.jobs)
  rho.list <- split.vec(rho, num.jobs)
  T.list <- split.vec(T, num.jobs)
  K.list <- split.vec(K, num.jobs)
  
  ## Save H.list[['0']] in file ".in0.H.txt" and so on
  for(ind in names(xi.list)){
    write.table(f10(xi.list[[ind]]), file=paste(path, ".in", ind, "xi.txt", sep=""), 
                col.names=FALSE, row.names=FALSE, sep="\n", quote=FALSE)
    write.table(f10(H.list[[ind]]), file=paste(path, ".in", ind, "H.txt", sep=""), 
                col.names=FALSE, row.names=FALSE, sep="\n", quote=FALSE)
    write.table(f10(eta.list[[ind]]), file=paste(path, ".in", ind, "eta.txt", sep=""), 
                col.names=FALSE, row.names=FALSE, sep="\n", quote=FALSE)
    write.table(f10(rho.list[[ind]]), file=paste(path, ".in", ind, "rho.txt", sep=""), 
                col.names=FALSE, row.names=FALSE, sep="\n", quote=FALSE)
    write.table(f10(T.list[[ind]]), file=paste(path, ".in", ind, "T.txt", sep=""), 
                col.names=FALSE, row.names=FALSE, sep="\n", quote=FALSE)
    write.table(f10(K.list[[ind]]), file=paste(path, ".in", ind, "K.txt", sep=""), 
                col.names=FALSE, row.names=FALSE, sep="\n", quote=FALSE)
  }
  
  ## method of variance reduction
  if(var.red)
    method <- "t"
  else
    method <- "s"
  
  ## vector of output file names and input filename modifiers
  in.names <- as.list(outer("in", names(xi.list), function(x,y) paste(x, y, sep="")))
  names(in.names) <- names(xi.list)
  out.names <- as.list(outer(".out", names(xi.list), function(x,y) paste(x, y, sep="")))
  names(out.names) <- names(xi.list)
  
  finf <- function(x) formatC(x, format="d")
  
  ## Construct each individual system command and concatenate them by paste
  sys.command.list <- lapply(names(xi.list), function(ind) paste(paste(path, 
                                            "Release/ML_console_ST", sep=""), 
                                             finf(N), finf(M), path, out.names[[ind]],
                                             in.names[[ind]], method))
  ## bugfix
  ##print(sys.command.list)
  
  sys.command <- paste(sys.command.list, collapse = " & ")
  system(sys.command)
  
  ## restore old options
  #options(scipen=scipen.old)
  #browser()
  ## read the result in again
  res.list <- lapply(names(xi.list), function(ind) read.table(
    paste(path, out.names[[ind]], sep=""), header=TRUE))
  ## Concatenate the data.frames
  res <- do.call("rbind", res.list)
  
  return(res)
}

## Richardson extrapolation based on rBergomi pricer
## N now denotes the smaller N, i.e., the extrapolation is based on N and 2*N
## Extrapolation is done on both prices and implied vols, i.e., the
## estimator for the implied vol does not correspond to the estimator for the price.
## As both estimators are independent, the variance of the estimator is assumed to be the
## sum of the variances.
## The extrapolation is based on the assumption that the true weak rate is 1.
rBergomi.richardson <- function(xi, H, eta, rho, T, K, N, M, num.jobs, var.red = FALSE){
  res1 <- rBergomi.pricer2(xi, H, eta, rho, T, K, N, M, num.jobs, var.red)
  res2 <- rBergomi.pricer2(xi, H, eta, rho, T, K, 2*N, M, num.jobs, var.red)
  res <- res1
  res$price <- 2*res1$price - res2$price
  res$iv <- 2*res1$iv - res2$iv
  res$stat <- sqrt(4*res1$stat^2 + res2$stat^2)
  return(res)
}