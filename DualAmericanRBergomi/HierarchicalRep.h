/*
 * HierarchicalRep.h
 *
 * Abstract base class encoding hierarchical representations of Brownian motion,
 * based on d independent Brownian motions and K subintervals of
 * [0,T].
 *
 * This class is responsible for simulating the standard normal random variables X_i, i=1, ..., I,
 * and the integrals of \phi_i^1 w.r.t. the fractional kernel K.
 * We may also want to provide the values of \tilde{W} and the increments of all Brownian motions
 * W^j along grids s_n which are refinements of the grid t_k. Both grids are assumed to be uniform.
 *
 * Note that only generate() will actually change the underlying random numbers!
 *
 *  Created on: 8 May 2018
 *      Author: bayerc
 */

#ifndef HIERARCHICALREP_H_
#define HIERARCHICALREP_H_

#include<vector>
#include<random>
#include<cmath>

typedef std::vector<double> Vector;
typedef std::mt19937_64 MTGenerator;
typedef std::normal_distribution<double> normDist;

class HierarchicalRep {
protected:
	double H; // Hurst index
	double T; // terminal time
	double sT; // 1/sqrt(T)
	int K; // number of exercise dates
	double dt; // mesh size of exercise grid
	int N; // size of the discretization grid. Assumed a multiple of K;
	double ds; // mesh size of the discretization grid.
	int d; // dimension of the Brownian motion
	long I; // size of the truncated basis
	Vector X; // vector for the standard norm numbers of size I
	MTGenerator gen;
	normDist dist;
	void setGen(std::vector<uint64_t> seed); // set up generators
public:
	HierarchicalRep();
	HierarchicalRep(double h, double t, int k, int n, int D, long i, std::vector<uint64_t> seed);
	// generate a new random sample X, do nothing else
	virtual void generate();
	// return X
	Vector getX() const;
	// return a sample of \tilde{W} on a grid of size N
	// note that this depends on the choice of hierarchical representation!
	// It is assumed that N is a multiple of K.
	// The sample starts from the value at 0, but does not include the value at t_N = T.
	virtual Vector Wtilde() const = 0;
	// return a sample of the increments of Brownian motion j along a grid of size N (assumed to be a multiple of K).
	// note that this depends on the choice of hierarchical representation!
	virtual Vector dB(int j) const = 0;
};

#endif /* HIERARCHICALREP_H_ */
