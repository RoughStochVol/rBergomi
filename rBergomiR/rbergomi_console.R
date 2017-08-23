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
rBergomi.pricer <- function(xi, H, eta, rho, T, K, N, M, num.threads, output.name){
  ## force non-scientific notation of numbers
  #scipen.old <- options()$scipen
  #options(scipen = 999)
  ## Convert number to string contining first 10 digits
  f10 <- function(x) formatC(x, digits=10, format="f")
  
  ## Save in files expected by console
  path <- "/Users/bayerc/Documents/workspace/roughBergomi/rBergomi/Release/"
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