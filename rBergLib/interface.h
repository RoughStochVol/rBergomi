/*
 * interface.h
 *
 * Provides simpler functions to be used as interface to R (and, potentially python).
 * Specifically, the return values are given by a simple struct which contains vectors of
 * the parameter values used (of the same size as the vector of resulting prices/vols).
 *
 * Note that there seems to be a problem with linking this library as built with g++  with
 * another code built with clang -> something related to std::vector??? Maybe we should already
 * dumb down everything here and replace the vectors by arrays? (This could also help the integration with
 * R and python).
 *
 *  Created on: May 13, 2017
 *      Author: bayerc
 */
#pragma once

#include<algorithm>

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

// appareently, we need simpler objects for input parameters
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

// test case
ResultSimple stupid(ParSimple x);
