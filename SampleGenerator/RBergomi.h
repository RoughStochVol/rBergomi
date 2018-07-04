/*
 * RBergomi.h
 *
 * Class for computing option prices in the R-Bergomi model in a single-threaded way.
 *
 * Input:
 * -) The parameters of the rBergomi model. Note that xi_0 can be a function of time.
 * -) Number of time steps N, and the time increment ds = T/N (but not T).
 * -) Seed of the random number generators.
 * -) Number of samples M to be generated.
 *
 * Output:
 * The function generate returns a tuple (S, W1, Wperp), where
 * -) S is of type std::vector<std::vector<double>> of size M x (N+1). S[m] corresponds to one trajectory of the
 *    underlying price process S.
 * -) W1, Wperp are of the same type as S, but of size M x N. W1[m] and Wperp[m] contain the increments of two
 *    independent Brownian motions driving the stock price. Note that the variance process only depends on
 *    W1, not on Wperp.
 *
 *  Created on: Nov 16, 2016
 *      Author: bayerc
 */

#ifndef RBERGOMI_H_
#define RBERGOMI_H_

#include<vector>
#include<random>
#include<cmath>
#include<iostream>
#include<fstream>
#include<iterator>
#include<string>
#include<numeric>
#include<tuple>
#include<fftw3.h>

#include"aux.h"
#include"RNorm.h"

typedef std::vector<double> Vector;
typedef std::vector<Vector> Matrix;
typedef std::tuple<Matrix, Matrix, Matrix> SampleCollection;

class RBergomi {
private:
	int N; // # of timesteps
	double ds; // timestep
	double S0; // spot price
	double r; // interest rate respectively negative dividend rate
	double H; // Hurst index
	double eta; // vol of vol
	double rho; // correlation
	double rhobar; // = sqrt(1 - rho^2)
	std::function<double(double)> xi0; // forward variance curve
	int nDFT; // length of arrays for which DFT needs to be computed; = 2*N - 1 ???
	RNorm ran;
	fftw_complex *xC; // one complex array (of size 2*N-1???)
	fftw_complex *xHat;
	fftw_complex *yC;
	fftw_complex *yHat;
	fftw_complex *zC;
	fftw_complex *zHat;
	fftw_plan fPlanX;
	fftw_plan fPlanY;
	fftw_plan fPlanZ;
	Vector Gamma;
	// some private methods
	// compute Wtilde as in terms of W1, W1perp, H by Pakkanenen et al.
	void updateWtilde(Vector& Wtilde, const Vector& W1, const Vector& W1perp);
	// compute Wtilde scaled to the interval [0,T]
	void scaleWtilde(Vector& WtildeScaled, const Vector& Wtilde) const;
	// compute Z scaled to the interval [0,T]
	void scaleZ(Vector& ZScaled, const Vector& Z) const;
	// update v
	void updateV(Vector& v, const Vector& WtildeScaled) const;
	// update Z
	void updateZ(Vector& Z, const Vector& W1, const Vector& Wperp) const;
	// compute S
	void updateS(Vector& S, const Vector& v, const Vector& Z) const;
	// Compute bstar
	Vector getGamma();
	// copy real vector to complex array; all remaining entries up to nDFT are set to 0
	void copyToComplex(const Vector& x, fftw_complex* xc);
	// copy from complex array to Vector (only fill up to length of x
	void copyToReal(Vector& x, const fftw_complex* xc) const;
	// pointwise complex multiplication z = x * y
	void fftw_c_mult(const fftw_complex x, const fftw_complex y, fftw_complex z);
	// element-wise complex multiplication z = x .* y
	void complexMult(const fftw_complex* x, const fftw_complex* y, fftw_complex* z);
	// Generate one sample
	void generateSample(Vector& S, Vector& W1, Vector& Wperp);
public:
	//*structors
	RBergomi();
	// seed is an optional parameter
	RBergomi(double S0I, double rI, double HI, double etaI, double rhoI,
			std::function<double(double)> xi0I, uint64_t seed, int NI, double dsI);
	~RBergomi();

	// methods
	SampleCollection generate(long M);
	SampleCollection operator()(long M);
};

#endif /* RBERGOMI_H_ */
