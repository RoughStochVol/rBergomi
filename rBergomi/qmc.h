/*
 * qmc.h
 *
 * Header file for implementing Sobol' numbers.
 *
 *  Created on: Nov 3, 2017
 *      Author: bayerc
 */

#pragma once

#include <vector>
#include <iostream>
#include <cstdlib>
#include <boost/math/distributions/normal.hpp>
#include "aux.h"
#include "sobol.hpp"

typedef std::vector<double> Vector;

const int maxDim = 1111;

// Generate an array of M normal random vectors of size N each.
// (I.e., M i.i.d. random samples of an N(0, I_N) random vector.)
// The samples are based on Sobol' numbers, with the assumption that
// N <= 1111. Note that two calls to this function with the same parameters
// will give back exactly the same results.
std::vector<Vector> normalQMC(long M, int N);

// Generate one sample of an N(0, I_N) random vector based on the i'th
// Sobol' vector in dimension N.
// Auxiliary function for normalQMC.
void normalQMC_sample(Vector& Z, int N, long long i);

// Standard normal quantile function.
inline double Phi_inv(double x){
	boost::math::normal dist(0.0, 1.0);
	return boost::math::quantile(dist, x);
}

// The same, but vectorized. Note that input vector is overwritten!
inline void Phi_inv(Vector& x){
	boost::math::normal dist(0.0, 1.0);
	for(size_t i = 0; i < x.size(); ++i)
		x[i] = boost::math::quantile(dist, x[i]);
}
