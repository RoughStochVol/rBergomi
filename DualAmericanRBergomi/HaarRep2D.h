/*
 * HaarRep2D.h
 *
 * Two-dimensional Haar basis for the rBergomi model. We follow the notation from the "plan".
 * As N is fixed, we can precompute once and for all the terms appearing in the formula for Wtilde
 * as well as in the formulas for the increments of the Brownian motions, so that we eventually
 * only need to compute cumulative sums. These precomputations are done in the constructor.
 *
 * Note that the basis functions are indexed by a multi-index \alpha = (l, j, i, k), where
 * - (l, j) denote the standard index of the underlying Haar basis, i.e., n denotes the level
 *   of the basis function (scale) and l the shift. Hence, we have l \in {0, ..., n}, j \in {0, ..., max(2^(l-1)-1, 0)}.
 * - i denotes the dimension, i.e., i \in {0, 1}.
 * - k denotes the sub-interval along the grid of exercise dates, k \in {0, ..., K-1}.
 * Following the order in the write-up, we need to map \iota \in {0, ..., I-1} to \alpha and conversely.
 * Note that I = K * 2 * (2^n + 1), and N = I / 2. Hence, K and I determine n.
 * We have: \iota = i + (K - k - 1)*d + K * d * int(2^(l-1) + j)
 * i = \iota % d; \iota_1 := (\iota - i)/d
 * k = K - (\iota_1 % K) - 1; \iota_2 := (\iota_1 - (K - k - 1)) / K = int(2^(l-1) + j)
 * if \iota_2 == 0 then l = j = 0;
 * else l = int(log_2(\iota_2)) + 1, j = \iota_2 - 2^(l-1)
 *
 *  Created on: 8 May 2018
 *      Author: bayerc
 */

#ifndef HAARREP2D_H_
#define HAARREP2D_H_

#include <array>
#include <iostream>
#include <numeric>
#include "HierarchicalRep.h"
#include "aux.h"

// as described above, multi-indices always have 4 elements here.
typedef std::array<int, 4> MultiIndex;

class HaarRep2D: public HierarchicalRep {
private:
	int n; // maximal level of the Haar basis used
	// auxiliary variables for all deterministic terms needed for computing Wtilde or dB.
	// these are vectors of size I, whose inner product with X determines the output.
	std::vector<Vector> WtildeKer;
	std::vector<Vector> dB1Ker;
	std::vector<Vector> dB2Ker;
	double sdt; // = sqrt(dt)
	// auxiliary functions for constructing those kernels
	void getWtildeKer();
	void getdBKer();
public:
	HaarRep2D();
	HaarRep2D(double h, double t, int k, int n, std::vector<uint64_t> seed);
	// The following functions should be private, but are public to allow testing.
	// transform multi-indices to indices and vice versa.
	int alpha2iota(MultiIndex alpha) const;
	MultiIndex iota2alpha(int iota) const;
	// The integrated Haar basis function based on [0, 1]
	double Psi(double t, int l, int j) const;
	// The integrated basis function for 2D based on the exercise grid. Return the i'th component.
	double Phi(double t, int iota, int i) const;
	// The Haar basis integrated against the fractional kernel;
	//double PsiTilde(double t, int l, int j) const;
	double PhiTilde(double t, int iota) const;
	// return a sample of \tilde{W} on a grid of size N
	// note that this depends on the choice of hierarchical representation!
	// It is assumed that N is a multiple of K.
	virtual Vector Wtilde() const;
	// return a sample of the increments of Brownian motion j along a grid of size N (assumed to be a multiple of K).
	// note that this depends on the choice of hierarchical representation!
	virtual Vector dB(int j) const;
};

#endif /* HAARREP2D_H_ */
