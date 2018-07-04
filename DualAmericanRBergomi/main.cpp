/*
 * main.cpp
 *
 * Check the convergence of the Haar approximation to the actual rBergomi model prices for European options.
 *
 *  Created on: 13 Jun 2018
 *      Author: bayerc
 */

#include <cmath>
#include "RBergomi.h"
#include "HaarRep2D.h"
#include "aux.h"

int main(){
	// model parameters
	const double H = 0.4;
	const double rho = -0.9;
	const double eta = 1.9;
	const double S0 = 1.0;
	const double r = 0.0;
	std::function<double(double)> xi = [](double t){return 0.04;};

	// option parameters
	const double T = 0.8;
	const double strike = 0.95;
	const int K = 1; // only one exercise possibility

	// numerical parameters
	const int numRuns = 5; // number of runs
	std::vector<int> Nvec = range<int>(2, 2 + numRuns);
	std::transform(Nvec.begin(), Nvec.end(), Nvec.begin(),
			[](int n){return static_cast<int>(pow(2, n));});
	std::vector<long> Mvec(numRuns);
	std::transform(Nvec.begin(), Nvec.end(), Mvec.begin(),
			[](int n){return n*n;});
	const std::vector<uint64_t> seed{1234677, 8525684, 26422367};
	Vector price(numRuns);
	Vector stat(numRuns);

	// reference price computed using the hybrid code.
	// Actually, use a different implementation for H = 0.5
	// Note: the hybrid scheme does not seem to work well for H = 1/2
	//const double refPrice = 0.09104141; // H = 1/2
	const double refPrice = 0.09434211; // H < 1/2 does not work at the moment

	// Generate objects
	HaarRep2D haar(H, T, K, Nvec[0], seed);
	RBergomi rBergomi(&haar, S0, r, eta, rho, xi);

	// Run the simulations
	for(int n = 0; n<numRuns; ++n){
		haar = HaarRep2D(H, T, K, Nvec[n], seed);
		Vector payoffs(Mvec[n]);
		for(long m=0; m<Mvec[n]; ++m)
			payoffs[m] = rBergomi.generateEuler().back(); // the asset price at the terminal time
		// now compute the payoffs
		std::transform(payoffs.begin(), payoffs.end(), payoffs.begin(),
				[strike](double s){return s > strike ? s - strike : 0.0;});
		// compute mean and variance of payoffs
		price[n] = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / Mvec[n];
		stat[n] = sqrt(std::inner_product(payoffs.begin(), payoffs.end(), payoffs.begin(), 0.0) / Mvec[n] -
				price[n] * price[n]) / sqrt(Mvec[n]);
	}

	Vector error(numRuns);
	std::transform(price.begin(), price.end(), error.begin(),
			[refPrice](double x) {return fabs(x - refPrice);});

	// report results
	std::cout << "N = " << Nvec << "\nPrice = " << price << "\nError = " << error << "\nStat = " << stat << std::endl;

	return 0;
}
