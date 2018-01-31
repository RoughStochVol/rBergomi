/*
 * rBergomi.h
 *
 *  Created on: 26 Jan 2018
 *      Author: bayerc
 */

#ifndef RBERGOMI_H_
#define RBERGOMI_H_

#include <vector>
#include <numeric>
#include <fftw3.h>
#include "RNorm.h"
#include "RfBm.h"

typedef std::vector<double> Vector;

/*
 Compute the call option price in the rough Bergomi model, with r = 0.
 */
void mc_bayer_roughbergomi(double S0, double eta, double H, double rho,
		double xi, double K, double T, int M, int N, double* price,
		double* stat);
void mc_bayer_roughbergomi_cholesky(double S0, double eta, double H, double rho,
		double xi, double K, double T, int M, int N, double* price,
		double* stat);
/*
 * Re-interpreting K as moneyness, we may have S0 = 1, without loss of generality.
 */
void mc_bayer_roughbergomi_moneyness(double eta, double H, double rho,
		double xi, double K, double T, int M, int N, double* price,
		double* stat);

/*
 * As above, but using Cholesky.
 */
void mc_bayer_roughbergomi_moneyness_cholesky(double eta, double H, double rho,
		double xi, double K, double T, int M, int N, double* price,
		double* stat);

// Sample (Wtilde, W1). Here, Wtilde is a fBm, and W corresponds to a correlated Brownian motion.
// More precisely, we have Wtilde[i] = \int_0^{t_{i+1}} K(t_{i+1},s) dW_s for t_i = i/N, i=0,..,N-1.
// On the other hand, W1[i] = (W_{t_{i+1}} - W_{t_i}) * \sqrt{N}, i.e., scaled to be
// standard normal.
// This function uses the hybrid scheme.
void sample_WtildeW_hybrid(Vector& Wtilde, Vector& W1, double H, Vector& Gamma,
		int N, RNorm& rnorm, int nDFT, fftw_complex* xC, fftw_complex* xHat, fftw_complex* yC,
		fftw_complex* yHat, fftw_complex* zC, fftw_complex* zHat,
		fftw_plan& fPlanX, fftw_plan& fPlanY, fftw_plan& fPlanZ);

// Use the hybrid scheme to generate Wtilde from standard Gaussian samples W1, W1perp
// Note that Wtilde here corresponds to samples of the (fBm) process Wtilde on a
// uniform grid of size N on [0,1], whereas W1, W1perp are i.i.d. N(0,1) gaussian vectors.
// Hence, both require rescaling.
void compute_Wtilde(Vector& Wtilde, const Vector& W1, const Vector& W1perp,
		double H, Vector& Gamma, int N, int nDFT, fftw_complex* xC,
		fftw_complex* xHat, fftw_complex* yC, fftw_complex* yHat,
		fftw_complex* zC, fftw_complex* zHat, fftw_plan& fPlanX,
		fftw_plan& fPlanY, fftw_plan& fPlanZ);

// Compute the Vector Gamma needed for the hybrid scheme
Vector getGamma(int N, double H);

// Compute value on one single trajectory
double updatePayoff(Vector& Wtilde, const Vector& W1, const Vector& W1perp,
		Vector& v, Vector& Gamma, fftw_complex* xC, fftw_complex* xHat,
		fftw_complex* yC, fftw_complex* yHat, fftw_complex* zC,
		fftw_complex* zHat, fftw_plan& fPlanX, fftw_plan& fPlanY,
		fftw_plan& fPlanZ, double eta, double H, double rho, double xi,
		double T, double K, int N);

// Compute value on one single trajectory
double updatePayoff_cholesky(Vector& Wtilde, const Vector& W1,
		Vector& v, double eta, double H, double rho, double xi,
		double T, double K, int N);

// compute v
void compute_V(Vector& v, const Vector& Wtilde, double H, double eta, double xi,
		double dt);

// compute \int_0^T v_t dt
double intVdt(const Vector & v, double dt);

// compute \int_0^T \sqrt{v_t} dW_t, with W = W1
double intRootVdW(const Vector & v, const Vector & W1, double sdt);

// BS Call price formula
double BS_call_price(double S0, double K, double tau, double sigma, double r =
		0.0);

// Auxiliary functions

// Normal CDF
double pnorm(double value);

// z = a*x+b*y
template<typename T>
inline std::vector<T> linearComb(T a, const std::vector<T>& x, T b,
		const std::vector<T>& y) {
	std::vector<T> z(x.size());
	for (size_t i = 0; i < z.size(); ++i)
		z[i] = a * x[i] + b * y[i];
	return z;
}

// x = s*x
template<typename T>
inline void scaleVector(std::vector<T>& x, T s) {
	for (auto& v : x)
		v = v * s;
}

// Copy real Vector to ftw_complex array. xc is of size nDFT, x is of size <= nDFT
void copyToComplex(const Vector& x, fftw_complex* xc, int nDFT);

// Copy real part of fftw_complex array to Vector
void copyToReal(Vector& x, const fftw_complex* xc);

// element-wise multiplication of arrays of fftw_complex numbers
void complexMult(const fftw_complex* x, const fftw_complex* y, fftw_complex* z,
		int nDFT);

// Multiply two fftw_complex numbers
void fftw_c_mult(const fftw_complex a, const fftw_complex b, fftw_complex c);

#endif /* RBERGOMI_H_ */

