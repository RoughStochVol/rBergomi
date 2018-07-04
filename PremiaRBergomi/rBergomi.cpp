/*
 * rBergomi.cpp
 *
 *  Created on: 26 Jan 2018
 *      Author: bayerc
 */

#include "rBergomi.h"

void mc_bayer_roughbergomi(double S0, double eta, double H, double rho,
		double xi, double K, double T, int M, int N, double* price,
		double* stat) {
	mc_bayer_roughbergomi_moneyness(eta, H, rho, xi, K / S0, T, M, N, price,
			stat);
	*price = *price * S0;
	*stat = *stat * S0;
}

void mc_bayer_roughbergomi_cholesky(double S0, double eta, double H, double rho,
		double xi, double K, double T, int M, int N, double* price,
		double* stat) {
	mc_bayer_roughbergomi_moneyness_cholesky(eta, H, rho, xi, K / S0, T, M, N, price,
			stat);
	*price = *price * S0;
	*stat = *stat * S0;
}

void mc_bayer_roughbergomi_moneyness(double eta, double H, double rho,
		double xi, double K, double T, int M, int N, double* price,
		double* stat) {
	// Prepare random number generator
	const uint64_t seed = 123456789;
	RNorm rnorm(seed);

	// Prepare FFT-Plans
	int nDFT = 2 * N - 1; // size of arrays for DFT.
	// generate arrays of complex numbers for DFT
	fftw_complex* xC = new fftw_complex[nDFT];
	fftw_complex* xHat = new fftw_complex[nDFT];
	fftw_complex* yC = new fftw_complex[nDFT];
	fftw_complex* yHat = new fftw_complex[nDFT];
	fftw_complex* zC = new fftw_complex[nDFT];
	fftw_complex* zHat = new fftw_complex[nDFT];
	// generate FFTW plans
	fftw_plan fPlanX = fftw_plan_dft_1d(nDFT, xC, xHat, FFTW_FORWARD,
	FFTW_ESTIMATE);
	fftw_plan fPlanY = fftw_plan_dft_1d(nDFT, yC, yHat, FFTW_FORWARD,
	FFTW_ESTIMATE);
	fftw_plan fPlanZ = fftw_plan_dft_1d(nDFT, zHat, zC, FFTW_BACKWARD,
	FFTW_ESTIMATE);

	// Compute Gamma
	Vector Gamma = getGamma(N, H);

	// Allocate memory for Gaussian random numbers and the random vector v (instantaneous variance)
	// Note that W1, W1perp correspond to UNNORMALIZED increments of Brownian motions,
	// i.e., are i.i.d. standard normal.
	Vector W1(N);
	Vector W1perp(N);
	Vector Wtilde(N);
	Vector v(N);

	double mean = 0.0; // will eventually be the mean
	double mu2 = 0.0; // will become second moment
	double var; // will eventually become variance (i.e., MC error).

	// The big loop which needs to be parallelized in future
	for (int m = 0; m < M; ++m) {
		// generate the fundamental Gaussians
		rnorm(W1);
		rnorm(W1perp); // Wperp is  not needed!

		double payoff = updatePayoff(Wtilde, W1, W1perp, v, Gamma, xC, xHat, yC,
				yHat, zC, zHat, fPlanX, fPlanY, fPlanZ, eta, H, rho, xi, T, K,
				N);
		mean += payoff;
		mu2 += payoff * payoff;
	}

// compute mean and variance
	mean = mean / M;
	mu2 = mu2 / M;
	var = mu2 - mean * mean;

// price = mean, stat = sqrt(var) / sqrt(M)
	*price = mean;
	*stat = sqrt(var / M);

// free arrays of complex numbers again
	delete[] xC;
	delete[] xHat;
	delete[] yC;
	delete[] yHat;
	delete[] zC;
	delete[] zHat;
}

void mc_bayer_roughbergomi_moneyness_cholesky(double eta, double H, double rho,
		double xi, double K, double T, int M, int N, double* price,
		double* stat) {
	// Prepare random number generator
	const uint64_t seed = 123456789;
	RNorm rnorm(seed);
	RfBm rfbm(N, H, &rnorm);

	// Allocate memory for Gaussian random numbers and the random vector v (instantaneous variance)
	// Note that W1, W1perp correspond to UNNORMALIZED increments of Brownian motions,
	// i.e., are i.i.d. standard normal.
	Vector W1(N);
	Vector Wtilde(N);
	Vector v(N);

	double mean = 0.0; // will eventually be the mean
	double mu2 = 0.0; // will become second moment
	double var; // will eventually become variance (i.e., MC error).

	// The big loop which needs to be parallelized in future
	for (int m = 0; m < M; ++m) {
		// generate W and Wtilde
		rfbm(W1, Wtilde);

		double payoff = updatePayoff_cholesky(Wtilde, W1, v, eta, H, rho, xi, T, K,
				N);
		mean += payoff;
		mu2 += payoff * payoff;
	}

// compute mean and variance
	mean = mean / M;
	mu2 = mu2 / M;
	var = mu2 - mean * mean;

// price = mean, stat = sqrt(var) / sqrt(M)
	*price = mean;
	*stat = sqrt(var / M);

}

void compute_Wtilde(Vector& Wtilde, const Vector& W1, const Vector& W1perp,
		double H, Vector& Gamma, int N, int nDFT, fftw_complex* xC,
		fftw_complex* xHat, fftw_complex* yC, fftw_complex* yHat,
		fftw_complex* zC, fftw_complex* zHat, fftw_plan& fPlanX,
		fftw_plan& fPlanY, fftw_plan& fPlanZ) {
	double s2H = sqrt(2.0 * H);
	double rhoH = s2H / (H + 0.5);
	Vector W1hat = linearComb(rhoH / s2H, W1, sqrt(1.0 - rhoH * rhoH) / s2H,
			W1perp);
	Vector Y2(N); // see R code
// Convolve W1 and Gamma
// Copy W1 and Gamma to complex arrays
	copyToComplex(W1, xC, nDFT);
	copyToComplex(Gamma, yC, nDFT);
// DFT both
	fftw_execute(fPlanX); // DFT saved in xHat[0]
	fftw_execute(fPlanY); // DFT saved in yHat[0]
// multiply xHat and yHat and save in zHat
	complexMult(xHat, yHat, zHat, nDFT);
// inverse DFT zHat
	fftw_execute(fPlanZ);
// read out the real part, re-scale by 1/nDFT
	copyToReal(Y2, zC);
	scaleVector(Y2, 1.0 / nDFT);
// Wtilde = (Y2 + W1hat) * sqrt(2*H) * dt^H ??
	Wtilde = linearComb(sqrt(2.0 * H) * pow(1.0 / N, H), Y2,
			sqrt(2.0 * H) * pow(1.0 / N, H), W1hat);
}

double updatePayoff(Vector& Wtilde, const Vector& W1, const Vector& W1perp,
		Vector& v, Vector& Gamma, fftw_complex* xC, fftw_complex* xHat,
		fftw_complex* yC, fftw_complex* yHat, fftw_complex* zC,
		fftw_complex* zHat, fftw_plan& fPlanX, fftw_plan& fPlanY,
		fftw_plan& fPlanZ, double eta, double H, double rho, double xi,
		double T, double K, int N) {
	double dt = T / N;
	double sdt = sqrt(dt);
	int nDFT = 2 * N - 1;
	compute_Wtilde(Wtilde, W1, W1perp, H, Gamma, N, nDFT, xC, xHat, yC, yHat,
			zC, zHat, fPlanX, fPlanY, fPlanZ);
	scaleVector(Wtilde, pow(T, H)); // scale Wtilde for time T
	compute_V(v, Wtilde, H, eta, xi, dt); // compute instantaneous variance v
// now compute \int v_s ds, \int \sqrt{v_s} dW_s with W = W1
	double Ivdt = intVdt(v, dt);
	double IsvdW = intRootVdW(v, W1, sdt);
// now compute the payoff by inserting properly into the BS formula
	double BS_vol = sqrt((1.0 - rho * rho) * Ivdt);
	double BS_spot = exp(-0.5 * rho * rho * Ivdt + rho * IsvdW);
	return BS_call_price(BS_spot, K, 1.0, BS_vol);
}

double updatePayoff_cholesky(Vector& Wtilde, const Vector& W1,
		Vector& v, double eta, double H, double rho, double xi,
		double T, double K, int N){
	double dt = T / N;
	double sdt = sqrt(dt);
	scaleVector(Wtilde, pow(T, H)); // scale Wtilde for time T
	compute_V(v, Wtilde, H, eta, xi, dt); // compute instantaneous variance v
// now compute \int v_s ds, \int \sqrt{v_s} dW_s with W = W1
	double Ivdt = intVdt(v, dt);
	double IsvdW = intRootVdW(v, W1, sdt);
// now compute the payoff by inserting properly into the BS formula
	double BS_vol = sqrt((1.0 - rho * rho) * Ivdt);
	double BS_spot = exp(-0.5 * rho * rho * Ivdt + rho * IsvdW);
	return BS_call_price(BS_spot, K, 1.0, BS_vol);
}

// Note that Wtilde plays the role of the old WtildeScaled!
void compute_V(Vector& v, const Vector& Wtilde, double H, double eta, double xi,
		double dt) {
	v[0] = xi;
	for (int i = 1; i < v.size(); ++i)
		v[i] = xi
				* exp(
						eta * Wtilde[i - 1]
								- 0.5 * eta * eta * pow(i * dt, 2 * H));
}

double intVdt(const Vector & v, double dt) {
	return dt * std::accumulate(v.begin(), v.end(), 0.0);
}

double intRootVdW(const Vector & v, const Vector & W1, double sdt) {
	double IsvdW = 0.0;
	for (size_t i = 0; i < v.size(); ++i)
		IsvdW += sqrt(v[i]) * sdt * W1[i];
	return IsvdW;
}

Vector getGamma(int N, double H) {
	Vector Gamma(N);
	double alpha = H - 0.5;
	Gamma[0] = 0.0;
	for (int i = 1; i < N; ++i)
		Gamma[i] = (pow(i + 1.0, alpha + 1.0) - pow(i, alpha + 1.0))
				/ (alpha + 1.0);
	return Gamma;
}

double pnorm(double value) {
	return 0.5 * erfc(-value * M_SQRT1_2);
}

double BS_call_price(double S0, double K, double tau, double sigma, double r) {
	double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * tau)
			/ (sigma * sqrt(tau));
	double d2 = d1 - sigma * sqrt(tau);
	return pnorm(d1) * S0 - pnorm(d2) * K * exp(-r * tau);
}

void copyToComplex(const Vector& x, fftw_complex* xc, int nDFT) {
	for (size_t i = 0; i < x.size(); ++i) {
		xc[i][0] = x[i]; // real part
		xc[i][1] = 0.0; // imaginary part
	}
// fill up with 0s
	for (size_t i = x.size(); i < nDFT; ++i) {
		xc[i][0] = 0.0; // real part
		xc[i][1] = 0.0; // imaginary part
	}
}

void copyToReal(Vector& x, const fftw_complex* xc) {
	for (size_t i = 0; i < x.size(); ++i)
		x[i] = xc[i][0]; // real part
}

void complexMult(const fftw_complex* x, const fftw_complex* y, fftw_complex* z,
		int nDFT) {
	for (size_t i = 0; i < nDFT; ++i)
		fftw_c_mult(x[i], y[i], z[i]);
}

void fftw_c_mult(const fftw_complex a, const fftw_complex b, fftw_complex c) {
	c[0] = a[0] * b[0] - a[1] * b[1];
	c[1] = a[0] * b[1] + a[1] * b[0];
}
