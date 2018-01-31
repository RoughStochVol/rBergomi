/*
 * fftOMP.h
 *
 * Get to the source of my FFT-bug!
 *
 *  Created on: 8 Dec 2017
 *      Author: bayerc
 */

#pragma once

#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <numeric>
#include <algorithm>
#include <fftw3.h>
#include <omp.h>
#include "aux.h"
#include "rBergomiMT.h"

typedef std::vector<double> Vector;
typedef std::vector<Vector> Array;

// The Single threaded convolution, duly vectorized Z[i] = X[i] * Y[i]
Array convolveST(const Array& X, const Array& Y);
// multi-threaded version
Array convolveMT(const Array& X, const Array& Y, int numThreads);

// Functions needed for FFT
void copyToComplexNew(const Vector& x, fftw_complex* xc, size_t nDFT);
// copy from complex array to Vector (only fill up to length of x
void copyToRealNew(Vector& x, const fftw_complex* xc);
// pointwise complex multiplication z = x * y
void fftw_c_multNew(const fftw_complex x, const fftw_complex y, fftw_complex z);
// element-wise complex multiplication z = x .* y
void complexMultNew(const fftw_complex* x, const fftw_complex* y, fftw_complex* z, size_t nDFT);

// Provide sample arrays for the convolution story
Array sampleArrayCos(size_t M, size_t N);
Array sampleArraySin(size_t M, size_t N);

// Compare Arrays row-wise
std::vector<bool> compareArrays(const Array& X, const Array& Y);
// compare floating point
bool pred(double x, double y); //used in above
