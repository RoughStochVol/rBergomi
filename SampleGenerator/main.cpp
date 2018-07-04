/*
 * main.cpp
 *
 *  Created on: 4 Jul 2018
 *      Author: bayerc
 */

#include<iostream>
#include"RBergomi.h"

int main(){
	const double H = 0.1;
	const double rho = -0.9;
	const double eta = 1.9;
	const double S0 = 1.0;
	const double r = 0.0;
	std::function<double(double)> xi = [](double t){return 0.04;};

	// option parameters
	const double T = 1.0;
	const double strike = 0.95;

	// numerical parameters
	const uint64_t seed = 12345678;
	const int N = 2000;
	double ds = T/N;
	const long M = 100000;

	const double refPrice = 0.09744944;

	RBergomi rberg(S0, r, H, eta, rho, xi, seed, N, ds);

	// Get samples from S
	auto STraj = std::get<0>(rberg(M));
	// extract the last elements of S
	Vector S(M);
	std::transform(STraj.begin(), STraj.end(), S.begin(), [](Vector x){return x.back();});
	// Compute the payoff
	auto payoff = [strike](double s){return s > strike ? s - strike : 0.0;};
	Vector payoffVec(M);
	std::transform(S.begin(), S.end(), payoffVec.begin(), payoff);

	// Compute mean and variance
	double price = std::accumulate(payoffVec.begin(), payoffVec.end(), 0.0) / M;
	double var = std::inner_product(payoffVec.begin(), payoffVec.end(), payoffVec.begin(), 0.0) / M - price * price;

	std::cout << "Price = " << price << ", error = " << fabs(price - refPrice)
			<< ", stat = " << sqrt(var) / sqrt(M) << std::endl;
}
