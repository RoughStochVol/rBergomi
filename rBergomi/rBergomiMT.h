/*
 * rBergomiMT.h
 *
 * For the multi-threaded version using OMP, we would need to copy
 * essentially all the data members in the rBergomi class. This does not make sense,
 * hence we do not use a class for this, but functions. We copy the functionality
 * already implemented in the single-threaded version.
 *  Created on: Feb 14, 2017
 *      Author: bayerc
 */

#ifndef RBERGOMIMT_H_
#define RBERGOMIMT_H_

#include<map>
#include<cstdlib>
#include"RBergomi.h"

// Define a struct for holding the complex arrays used by fftw3.
// We actually use pointers to pointers since we implement one array
// for each thread.
// ATTENTION: We have seen before that managing the copied resources
// manually like this is not very efficient! I do not think that this would
// be easy to do because of the fftw_plan-problem....
struct fftData{
	int numThreads;
	fftw_complex **xC; // one complex array (of size nDFT = 2*N-1) for each thread
	fftw_complex **xHat;
	fftw_complex **yC;
	fftw_complex **yHat;
	fftw_complex **zC;
	fftw_complex **zHat;
	std::vector<fftw_plan> fPlanX; // one for each thread
	std::vector<fftw_plan> fPlanY; // one for each thread
	std::vector<fftw_plan> fPlanZ; // one for each thread
	// constructor and destructor
	fftData(int nDFT, int numThreads);
	~fftData();
};

// Define a struct for the random number generators and the distribution
struct RNG{
	std::vector<MTGenerator> gen; // one for each thread
	normDist dist; // one is enough
	RNG(int numThreads, std::vector<uint64_t> seed);
};

// compute Wtilde as in terms of W1, W1perp, H by Pakkanenen et al.
void updateWtilde(Vector& Wtilde, const Vector& W1, const Vector& W1perp, double H,
		fftData& fft, int threadID, int nDFT);

// compute Wtilde scaled to the interval [0,T]
void scaleWtilde(Vector& WtildeScaled, const Vector& Wtilde, double T, double H);

// compute Z scaled to the interval [0,T]
void scaleZ(Vector& ZScaled, const Vector& Z, double sdt);

// update v
void updateV(Vector& v, const Vector& WtildeScaled, double xi, double H, double eta, double dt);

// update Z
void updateZ(Vector& Z, const Vector& W1, const Vector& Wperp, double rho);

// compute S
double updateS(const Vector& v, const Vector& Z, double dt);

// Compute bstar
void getGamma(Vector& Gamma, double H);

// copy real vector to complex array; all remaining entries up to nDFT are set to 0
void copyToComplex(int nDFT, const Vector& x, fftw_complex* xc);

// copy from complex array to Vector (only fill up to length of x
void copyToReal(Vector& x, const fftw_complex* xc);

// pointwise complex multiplication z = x * y
void fftw_c_mult(const fftw_complex x, const fftw_complex y, fftw_complex z);

// element-wise complex multiplication z = x .* y
void complexMult(int nDFT, const fftw_complex* x, const fftw_complex* y, fftw_complex* z);

// generate Gaussians
void genGaussianMT(Vector& X, RNG& rng, int threadID);

// Compute the price
Result ComputePriceMT(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K, int N, long M,
			int numThreads, std::vector<uint64_t> seed);

// Compute implied vols
Result ComputeIVMT(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K, int N, long M,
			int numThreads, std::vector<uint64_t> seed);

// Compute the price (with Romano-Touzi trick)
Result ComputePriceRTMT(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K, int N, long M,
			int numThreads, std::vector<uint64_t> seed);

// Compute implied vols (with Romano-Touzi trick)
Result ComputeIVRTMT(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K, int N, long M,
			int numThreads, std::vector<uint64_t> seed);

// Compute implied vols (with Romano-Touzi trick) based on quadrature points provided from user
// W1 and W1perp contain the quadrature points, the weights are actually not needed here.
std::vector<Vector> ComputePayoffRTsamples(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K,
		int numThreads, const std::vector<Vector> & W1Arr, const std::vector<Vector> & W1perpArr);

// For bug-fixing a single threaded variant...
std::vector<Vector> ComputePayoffRTsamples_ST(double xi, Vector H, Vector eta,
		Vector rho, Vector T, Vector K, const std::vector<Vector> & W1Arr,
		const std::vector<Vector> & W1perpArr);

// auxiliary functions needed by Romano-Touzi approach
// integrate v w.r.t. time
double intVdt(const Vector& v, double dt);

// integrate sqaure root of v w.r.t. dW, with \Delta W = sdt * W1
double intRootVdW(const Vector& v, const Vector& W1, double sdt);


#endif /* RBERGOMIMT_H_ */
