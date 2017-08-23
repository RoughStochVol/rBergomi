/*
 * Here we provide the R-interface to the rBergomi library. Essentially, we wrap
 * input and outputs for the functions declared in interface.h.
 * 
 * 
 */
#include <vector>
#include <Rcpp.h>
//#include "interface.h"
#include "cinterface.h"

using namespace Rcpp;

// try to simply use the adress of the first element of a NumericVector
// as pointer
// [[Rcpp::export]]
DataFrame computeIV(double xi, NumericVector H, NumericVector eta, 
                    NumericVector rho, NumericVector T, NumericVector K, 
                    int N, long M, IntegerVector seed, int num_threads){
  // construct vectors for the solutions
  int par_size = H.size();
  int op_size = T.size();
  int out_size = par_size * op_size;
  NumericVector price(out_size);
  NumericVector iv(out_size);
  NumericVector stat(out_size);
  NumericVector H_tot(out_size);
  NumericVector eta_tot(out_size);
  NumericVector rho_tot(out_size);
  NumericVector T_tot(out_size);
  NumericVector K_tot(out_size);
  // Automatic conversion does not seem to work between R-integer and uint64_t
  std::vector<uint64_t> seedv(seed.begin(), seed.end());
  
  // call the C-style interface, assuming contiguous memory in NumericVector
  computeIV_tensor(xi, &(H[0]), &(eta[0]), &(rho[0]), par_size, &(T[0]), &(K[0]), 
                   op_size, N, M, &(seedv[0]), seed.size(), num_threads, 
                   &(price[0]), &(iv[0]), &(stat[0]), &(H_tot[0]), &(eta_tot[0]),
                   &(rho_tot[0]), &(T_tot[0]), &(K_tot[0]), out_size);
                   
  // construct dataframe
  return DataFrame::create(_["price"]= price, _["iv"]= iv, _["H"]= H_tot, 
                           _["eta"]= eta_tot, _["rho"]= rho_tot, _["T"]= T_tot,
                           _["K"]= K_tot, _["xi"]= xi, _["stat"]= stat, _["N"]= N, 
                           _["M"]= M, _["num.threads"]= num_threads);
}

// [[Rcpp::export]]
DataFrame computeIV_RT(double xi, NumericVector H, NumericVector eta, 
                    NumericVector rho, NumericVector T, NumericVector K, 
                    int N, long M, IntegerVector seed, int num_threads){
  // construct vectors for the solutions
  int par_size = H.size();
  int op_size = T.size();
  int out_size = par_size * op_size;
  NumericVector price(out_size);
  NumericVector iv(out_size);
  NumericVector stat(out_size);
  NumericVector H_tot(out_size);
  NumericVector eta_tot(out_size);
  NumericVector rho_tot(out_size);
  NumericVector T_tot(out_size);
  NumericVector K_tot(out_size);
  // Automatic conversion does not seem to work between R-integer and uint64_t
  std::vector<uint64_t> seedv(seed.begin(), seed.end());
  
  // call the C-style interface, assuming contiguous memory in NumericVector
  computeIVRT_tensor(xi, &(H[0]), &(eta[0]), &(rho[0]), par_size, &(T[0]), &(K[0]), 
                   op_size, N, M, &(seedv[0]), seed.size(), num_threads, 
                   &(price[0]), &(iv[0]), &(stat[0]), &(H_tot[0]), &(eta_tot[0]),
                   &(rho_tot[0]), &(T_tot[0]), &(K_tot[0]), out_size);
                   
                   // construct dataframe
                   return DataFrame::create(_["price"]= price, _["iv"]= iv, _["H"]= H_tot, 
                                            _["eta"]= eta_tot, _["rho"]= rho_tot, _["T"]= T_tot,
                                              _["K"]= K_tot, _["xi"]= xi, _["stat"]= stat, _["N"]= N, 
                                              _["M"]= M, _["num.threads"]= num_threads);
}