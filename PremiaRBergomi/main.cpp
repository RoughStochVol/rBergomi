/*
 * main.cpp
 *
 *  Created on: 26 Jan 2018
 *      Author: bayerc
 */

#include <iostream>
#include "rBergomi.h"

template <typename T>
inline std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  if ( !v.empty() ) {
    out << '(';
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b)";
  }
  return out;
}

int main() {

	// sample parameters
	const double xi = 0.07;
	const double H = 0.07;
	const double eta = 2.2;
	const double rho = -0.9;
	const double T = 1.0;
	const double K = 0.9;
	const double S0 = 1.0;
	const double ref = 0.2149403;

	double price, stat;
	double price_cholesky, stat_cholesky;

	int N = 200;
	int M = 500000;

	/*
	RNorm rnorm(12345678);
	RfBm rfbm(5, H, &rnorm);

	std::vector<Vector> L = rfbm.GetL();
	std::vector<Vector> A = rfbm.GetA();
	std::cout << "Covariance matrix:\n";
	for(size_t i = 0; i<L.size(); ++i)
		std::cout << A[i] << "\n"; // Does not seem positive definite!!!!

	std::cout << "\nCholesky decomposition of the covariance matrix:\n";
	for(size_t i = 0; i<L.size(); ++i)
		std::cout << L[i] << "\n";
		*/


	mc_bayer_roughbergomi(S0, eta, H, rho, xi, K, T, M, N, &price, &stat);
	std::cout << "The computed price is " << price << " with standard deviation " << stat
				<< ".\nThe reference price is " << ref << ".\n";

	mc_bayer_roughbergomi_cholesky(S0, eta, H, rho, xi, K, T, M, N, &price_cholesky, &stat_cholesky);
	std::cout << "Using Cholesky, the computed price is " << price_cholesky << " with standard deviation " << stat_cholesky
				<< ".\nThe reference price is " << ref << ".\n";


	return 0;
}
