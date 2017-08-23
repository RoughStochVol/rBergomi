// Test the shared library version of the rBergomi code

#include<iostream>
#include<fstream>
#include<chrono>

#include"aux.h"
#include"interface.h"
#include "cinterface.h"

using Clock = std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

int main(){

	const double xi = 0.07;
	Vector H { 0.07, 0.1 };
	Vector eta { 2.2, 2.4 };
	Vector rho { -0.9, -0.85 };
	Vector T { 0.5, 1.0 };
	Vector K { 0.8, 1.0 };
	Vector truePrice { 0.2149403, 0.08461882, 0.2157838, 0.07758564 };
	std::vector<int> steps { 4, 8, 16, 32, 64};//, 128, 256, 512, 1024};
	std::vector<long> M(steps.size());
	// asymptotically, we expect to use M ~ steps^(1+2H)
	for(size_t i=0; i<M.size(); ++i)
		M[i] = 80*static_cast<long>(pow(static_cast<double>(steps[i]), 2.0)); // this is in order to make the stat. error much smaller
	const int numThreads = 4;

	std::string fname = "temp.txt";
	std::fstream file;
	file.open(fname.c_str(), std::fstream::out | std::fstream::trunc);
	file
			<< "steps M price1 price2 price3 price4 iv1 iv2 iv3 iv4 stat1 stat2 stat3 stat4 numthreads time(ms)\n";
	file.close();

	// for C-interface, we also need the return values.
	Vector price(H.size() * T.size());
	Vector iv(H.size() * T.size());
	Vector stat(H.size() * T.size());
	Vector H_tot(H.size() * T.size());
	Vector eta_tot(H.size() * T.size());
	Vector rho_tot(H.size() * T.size());
	Vector T_tot(H.size() * T.size());
	Vector K_tot(H.size() * T.size());

	std::vector<uint64_t> seed { 1, 2, 3, 4, 5, 6, 7, 8 };

	std::cout << "RBergomi-pricer with " << numThreads << " threads." << std::endl;

	for (size_t i = 0; i < steps.size(); ++i) {
		auto start = Clock::now();
		//ResultSimple resPar1 = ComputeIV_multiThread(xi, H, eta, rho, T, K, steps[i], M[i], numThreads, seed);
		computeIV_tensor(xi, &H[0], &eta[0], &rho[0], H.size(), &T[0], &K[0], T.size(),
				steps[i], M[i], &seed[0], seed.size(), numThreads, &price[0], &iv[0],
				&stat[0], &H_tot[0], &eta_tot[0], &rho_tot[0], &T_tot[0], &K_tot[0],
				H.size() * T.size());

		auto end = Clock::now();
		auto diffPar = duration_cast<milliseconds> (end - start);
		// put results in ResultSimple in order to re-use following code
		ResultSimple resPar2{price, iv, H_tot, eta_tot, rho_tot, T_tot, K_tot,
			xi, stat, steps[i], M[i], numThreads};

		// Print the results
		std::cout << '\n';
/*
		std::cout << "C++ Interface: \n";
		std::cout << "Price = " << resPar1.price << "\nIV = " << resPar1.iv
				<< "\nStatistical error = "
				<< resPar1.stat << "\nExact price = " << truePrice
				<< "\nH = " << resPar1.H
				<< "\neta = " << resPar1.eta
				<< "\nrho = " << resPar1.rho
				<< "\nTime elapsed = " << diffPar.count() << "ms" << std::endl;*/
		std::cout << "C Interface: \n";
		std::cout << "Price = " << resPar2.price << "\nIV = " << resPar2.iv
						<< "\nStatistical error = "
						<< resPar2.stat << "\nExact price = " << truePrice
						<< "\nH = " << resPar2.H
						<< "\neta = " << resPar2.eta
						<< "\nrho = " << resPar2.rho
						<< "\nTime elapsed = " << diffPar.count() << "ms" << std::endl;

		file.open(fname.c_str(), std::fstream::out | std::fstream::app);
		file << steps[i] << " " << M[i] << " " << resPar2.price[0] << " "
				<< resPar2.price[1] << " " << resPar2.price[2] << " "
				<< resPar2.price[3] << " "<< resPar2.iv[0] << " "
				<< resPar2.iv[1] << " " << resPar2.iv[2] << " "
				<< resPar2.iv[3] << " " << resPar2.stat[0] << " "
				<< resPar2.stat[1] << " " << resPar2.stat[2] << " "
				<< resPar2.stat[3] << " " << resPar2.numThreads << " "
				<< diffPar.count() << std::endl;
		file.close();
	}

/*
	double H[2] = {1.0,2.9};
	ParSimple par(H, 2);
	ResultSimple x = stupid(par);
	std::cout << x.H[0] << " " << x.H[1] << std::endl;
	//double x = stupidS(H);
	//std::cout << x << std::endl;
	return 0;
	*/
}
