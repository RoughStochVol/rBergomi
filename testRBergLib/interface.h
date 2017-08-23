/*
 * interface.h
 *
 * Provides simpler functions to be used as interface to R (and, potentially python).
 * Specifically, the return values are given by a simple struct which contains vectors of
 * the parameter values used (of the same size as the vector of resulting prices/vols).
 *
 * We just add declarations for the functions used, thereby avoiding the need for including the full headers -- in particular, omp.h.
 *
 *  Created on: May 13, 2017
 *      Author: bayerc
 */
#pragma once

#include<vector>
#include"BlackScholes.h"
#include"ParamTot.h"

typedef std::vector<double> Vector;

struct ResultSimple {
	Vector price; // prices for each combination of (H, eta, rho) and (T,K)
	Vector iv; // implied vols; note that we always fill both of them
	Vector H; // The same size as price. H[i] is used for price[i], and similarly for the other parameters.
	Vector rho;
	Vector eta;
	Vector T;
	Vector K;
	double xi;
	Vector stat; // statistical errors
	int N;
	long M;
	int numThreads;
};

struct Result {
	Vector price; // prices for each combination of (H, eta, rho) and (T,K)
	Vector iv; // only one of price and iv is computed! The other is set to NA or something like it.
	ParamTot par; // vectors of all the parameter combinations corresponding to prices and iv's
	Vector stat; // statistical errors
	int N;
	long M;
	int numThreads;
	double time; // run time? Maybe more sophisticated?
};

ResultSimple ComputeIV_singleThread(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K, int N, long M,
			std::vector<uint64_t> seed);

ResultSimple ComputeIV_multiThread(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K, int N, long M, int numThreads,
			std::vector<uint64_t> seed);

// Both functions with the Romano-Touzi-trick
ResultSimple ComputeIVRT_singleThread(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K, int N, long M,
			std::vector<uint64_t> seed);

ResultSimple ComputeIVRT_multiThread(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K, int N, long M, int numThreads,
			std::vector<uint64_t> seed);

// auxiliary function doing the work
ResultSimple Result2ResultSimple(const Result& r);


struct ParSimple{
	// a C-style array for H
	double *H;
	size_t n_H;
	// constructor
	ParSimple();
	ParSimple(double* h, size_t n_h);
	~ParSimple();
	// output the corresponding vectors
	Vector getH() const;
};
// test this
ResultSimple stupid(ParSimple x);
