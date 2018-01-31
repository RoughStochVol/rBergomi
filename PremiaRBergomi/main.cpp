/*
 * main.cpp
 *
 *  Created on: 26 Jan 2018
 *      Author: bayerc
 */

#include <iostream>
#include "rBergomi.h"

int main() {

	// sample parameters
	const double xi = 0.07;
	const double H = 0.07;
	const double eta = 2.2;
	const double rho = -0.9;
	const double T = 0.5;
	const double K = 0.8;
	const double S0 = 1.0;
	const double ref = 0.2149403;

	double price, stat;
	double price_cholesky, stat_cholesky;

	int N = 500;
	int M = 100000;

	mc_bayer_roughbergomi(S0, eta, H, rho, xi, K, T, M, N, &price, &stat);
	std::cout << "The computed price is " << price << " with standard deviation " << stat
				<< ".\nThe reference price is " << ref << ".\n";

	mc_bayer_roughbergomi_cholesky(S0, eta, H, rho, xi, K, T, M, N, &price_cholesky, &stat_cholesky);
	std::cout << "Using Cholesky, the computed price is " << price_cholesky << " with standard deviation " << stat_cholesky
				<< ".\nThe reference price is " << ref << ".\n";

	return 0;
}
