/*
 * BlackScholes.h
 *
 *  Created on: Apr 13, 2017
 *      Author: bayerc
 */

#pragma once

#include<cmath>
#include <boost/math/tools/roots.hpp>

// [[Rcpp::export]]
double Phi(double value);

double BS_call_price(double S0, double K, double tau, double sigma, double r = 0.0);

double BS_put_price(double S0, double K, double tau, double sigma, double r = 0.0);

double IV_call(double price, double S0, double K, double tau, double r = 0.0);

double IV_put(double price, double S0, double K, double tau, double r = 0.0);

/*
 * Error message:
 * Error in dyn.load(file, DLLpath = DLLpath, ...) : 
 unable to load shared object '/Library/Frameworks/R.framework/Versions/3.3/Resources/library/rBergomi/libs/rBergomi.so':
dlopen(/Library/Frameworks/R.framework/Versions/3.3/Resources/library/rBergomi/libs/rBergomi.so, 6): Library not loaded: librBergLib.dylib
Referenced from: /Library/Frameworks/R.framework/Versions/3.3/Resources/library/rBergomi/libs/rBergomi.so
Reason: image not found
Error: loading failed
Execution halted
ERROR: loading failed
* removing '/Library/Frameworks/R.framework/Versions/3.3/Resources/library/rBergomi'
Warning in install.packages :
installation of package '/Users/bayerc/Documents/workspace/roughBergomi/rBergomiR/rBergomi' had non-zero exit status
 */