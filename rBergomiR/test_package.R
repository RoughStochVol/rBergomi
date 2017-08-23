###############################################################
## Installaation instructions:
## 1) Compile the dynamic library librBergLib.dynlib with a compiler supporting
##    openmp. We use a plain C interface, which should work across compilers.
## 2) Copy the library into ~/local/lib
## 3) Add ~/local/lib to the library serach path. For R, this id one by adding 
##    the following to .Rprofile:
##    .libPaths(c(.libPaths(), "~/local/lib"))
##    You can check success by calling .libPaths(). Note that RStudio by default
##    does not source .Rprofile, unless it is explicitly run from the terminal.
##    I also found it necessary to add ~/local/lib to the library serach path of the
##    shell. Under OSX, I add the following to .bash_profile
##    export DYLD_FALLBACK_LIBRARY_PATH=${HOME}/local/lib/:${DYLD_FALLBACK_LIBRARY_PATH}
## THIS STILL DOES N OT WORK! UNLIKE WITH LINUX...
## CLAIMS NBOT TO FIND THE LIBRARY...
## Question: How to find out the library search path during install.packages??


# ## Test the hello-word package
library(Rcpp)
# 
# ## construct automatically generated functions
compileAttributes("./RBergomi")
# 
# ## Remove it again
detach("package:RBergomi", unload=TRUE)
remove.packages("RBergomi")
install.packages("./RBergomi", repos=NULL, type="source", clean=TRUE)

#dyn.load("/home/bayerc/workspace/roughBergomi/rBergLib/Release/librBergLib.so")
dyn.load("~/local/lib/librBergLib.dylib")

## Test the library
library(RBergomi)

## Error message:
## Error in .Call("rBergomi_stupid_two", PACKAGE = "rBergomi") : 
## "rBergomi_stupid_two" not available for .Call() for package "rBergomi"
## But: nm -gU applied to rBergomi.so shows that rcpp_hello_word() is exposed
## by the library exactly the same way as stupid_two(), for instance.
## Note that you may have to restart R or somehow detach the underlying
## shared library for changes to get into effect.

xi <- 0.04
H <- c(0.05, 0.1, 0.15)
eta <- c(1.5, 2.5, 3.0)
rho <- c(-0.8, -0.9, -0.9)
T <- c(0.5, 1.0)
K <- c(1.0, 1.0)
N <- 100
M <- 100
num.threads <- 1
seed <- c(1234, 5678)
res <- computeIV(xi, H, eta, rho, T, K, N, M, seed, num.threads)
res
# eta and rho are switched??

xi <- 0.07
H <- c(0.07, 0.1)
eta <- c(2.2, 2.4)
rho <- c(-0.9, -0.85)
T <- c(0.5, 1.0 )
K <- c(0.8, 1.0)
truePrice <- c(0.2149403, 0.08461882, 0.2157838, 0.07758564)

N <- 1000
M <- 10000
num.threads <- 1
seed <- c(1234, 5678)
res <- computeIV(xi, H, eta, rho, T, K, N, M, seed, num.threads)
res

N <- 1000
M <- 10000
num.threads <- 1
seed <- c(1234, 5678)
res <- computeIV_RT(xi, H, eta, rho, T, K, N, M, seed, num.threads)
res

N <- 1000
M <- 10000
num.threads <- 4
seed <- c(1234, 5678)
res <- computeIV(xi, H, eta, rho, T, K, N, M, seed, num.threads)
res
## Results look reasonable
## For some reason, computeIV_RT even gets the eta and rho right in the output file!