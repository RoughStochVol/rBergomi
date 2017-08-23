/*
 * interface.h
 *
 * Provides simpler functions to be used as interface to R (and, potentially python).
 * Specifically, the return values are given by a simple struct which contains vectors of
 * the parameter values used (of the same size as the vector of resulting prices/vols).
 *
 * Also use a simple C-style interface for easier compatability.
 *
 *  Created on: May 13, 2017
 *      Author: bayerc
 */
#pragma once

#include"RBergomi.h"
#include"rBergomiMT.h"

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
