/*
 * RBergomiST.h
 *
 * Class for computing option prices in the R-Bergomi model in a single-threaded way.
 *
 * Input:
 * -) S_0 and xi_0 are known constants (xi_0 is assumed to be constant!), r is assumed to be 0 (forward prices).
 *    Actually, S_0 = 1. (K in terms of moneyness).
 * -) Values (T,K) of maturities and strikes for which prices are to be computed.
 *    These are assumed to be ordered according to T.
 * -) Parameter values (H, eta, rho), for which prices are to be computed.
 *    These are assumed to be ordered according to H, then eta, then rho.
 * -) Number of time steps N and number of samples M. (Note that the same samples of the underlying Gaussian
 *    processes are used for all different (T,K) and all different (H, eta, rho).)
 * -) Seeds of the random number generators, to be determined.
 *
 * Output:
 * A vector of call prices for all combinations of (H, eta, rho) and (T,K).
 * We implement both the standard Monte Carlo algorithm and the one based on the Romano-Touzi argument.
 *
 * In a class we cannot really parallelize using openmp! Do it now in the class for a single-thread case, and then
 * define a single standing function for the general case.
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
#include<fftw3.h>
#include<omp.h>

#include"ParamTot.h"
#include"aux.h"
#include"BlackScholes.h"

typedef std::vector<double> Vector;
typedef std::mt19937_64 MTGenerator;
typedef std::normal_distribution<double> normDist;

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

class RBergomiST {
private:
	double xi;
	Vector H;
	Vector eta;
	Vector rho;
	Vector T;
	Vector K;
	int N;
	int nDFT; // length of arrays for which DFT needs to be computed; = 2*N - 1 ???
	long M;
	int numThreads;
	MTGenerator gen;
	normDist dist;
	fftw_complex *xC; // one complex array (of size 2*N-1???)
	fftw_complex *xHat;
	fftw_complex *yC;
	fftw_complex *yHat;
	fftw_complex *zC;
	fftw_complex *zHat;
	fftw_plan fPlanX;
	fftw_plan fPlanY;
	fftw_plan fPlanZ;
	// some private methods
	void setGen(std::vector<uint64_t> seed); // set up generators
	// compute Wtilde as in terms of W1, W1perp, H by Pakkanenen et al.
	void updateWtilde(Vector& Wtilde, const Vector& W1, const Vector& W1perp, double H);
	// compute Wtilde scaled to the interval [0,T]
	void scaleWtilde(Vector& WtildeScaled, const Vector& Wtilde, double T, double H) const;
	// compute Z scaled to the interval [0,T]
	void scaleZ(Vector& ZScaled, const Vector& Z, double sdt) const;
	// update v
	void updateV(Vector& v, const Vector& WtildeScaled, double h, double e, double dt) const;
	// update Z
	void updateZ(Vector& Z, const Vector& W1, const Vector& Wperp, double r) const;
	// compute S
	double updateS(const Vector& v, const Vector& Z, double dt) const;
	// Compute bstar
	void getGamma(Vector& Gamma, double H) const;
	// copy real vector to complex array; all remaining entries up to nDFT are set to 0
	void copyToComplex(const Vector& x, fftw_complex* xc);
	// copy from complex array to Vector (only fill up to length of x
	void copyToReal(Vector& x, const fftw_complex* xc) const;
	// pointwise complex multiplication z = x * y
	void fftw_c_mult(const fftw_complex x, const fftw_complex y, fftw_complex z);
	// element-wise complex multiplication z = x .* y
	void complexMult(const fftw_complex* x, const fftw_complex* y, fftw_complex* z);
	// integrate v w.r.t. time
	double intVdt(const Vector& v, double dt) const;
	// integrate sqaure root of v w.r.t. dW, with \Delta W = sdt * W1
	double intRootVdW(const Vector& v, const Vector& W1, double sdt) const;
public:
	//*structors
	RBergomiST();
	// seed is an optional parameter
	RBergomiST(double x, Vector HIn, Vector e, Vector r, Vector t, Vector k, int NIn, long MIn,
			std::vector<uint64_t> seed);
	~RBergomiST();

	// methods
	//Result ComputePrice();
	Result ComputePrice(); // Single threaded version
	//Result ComputeIV();
	Result ComputeIV(); // Single threaded version
	//Result ComputePrice();
	Result ComputePriceRT(); // Single threaded version, using Romano-Touzi
	//Result ComputeIV();
	Result ComputeIVRT(); // Single threaded version, using Romano-Touzi
	// auxiliary methods
	void genGaussian(Vector& X); // generate Gaussian vector X.
	//void genGaussianMT(Vector& X); // generate Gaussian vector X. Use generator corresponding to right thread.
	// getters and setters
	long getM() const;
	void setM(long m);
	int getN() const;
	void setN(int n);
	int getNumThreads() const;
	double getXi() const;
	void setXi(double xi);
	// test the generation of Wtilde; simulate Wtilde and W and put in file
	void testScaleWtilde();
	void testWtilde();
	void testConvolve();
	//void test_genGaussian();
};

#endif /* RBERGOMI_H_ */
