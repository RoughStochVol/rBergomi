/*
 * RBergomiST.cpp
 *
 *  Created on: Nov 16, 2016
 *      Author: bayerc
 */

#include "RBergomi.h"

RBergomi::RBergomi() {
	S0 = 0.0;
	r = 0.0;
	H = 0.0;
	rho = 0.0;
	rhobar = 0.0;
	eta = 0.0;
	xi0 = [](double t) {return 0.0;};
	ran = RNorm(0);
	N = 0;
	ds = 0.0;
	nDFT = 0;
	xC = new fftw_complex[0];
	xHat = new fftw_complex[0];
	yC = new fftw_complex[0];
	yHat = new fftw_complex[0];
	zC = new fftw_complex[0];
	zHat = new fftw_complex[0];
	Gamma = Vector(0);
}

// seed is optional
RBergomi::RBergomi(double S0I, double rI, double HI, double etaI, double rhoI,
		std::function<double(double)> xi0I, uint64_t seed, int NI, double dsI) {
	// Some safety tests/asserts?
	S0 = S0I;
	r = rI;
	H = HI;
	eta = etaI;
	rho = rhoI;
	rhobar = sqrt(1.0 - rho * rho);
	xi0 = xi0I;
	ran = RNorm(seed);
	ds = dsI;
	N = NI;
	nDFT = 2 * N - 1;
	xC = new fftw_complex[nDFT];
	xHat = new fftw_complex[nDFT];
	yC = new fftw_complex[nDFT];
	yHat = new fftw_complex[nDFT];
	zC = new fftw_complex[nDFT];
	zHat = new fftw_complex[nDFT];
	fPlanX = fftw_plan_dft_1d(nDFT, xC, xHat, FFTW_FORWARD,
	FFTW_ESTIMATE);
	fPlanY = fftw_plan_dft_1d(nDFT, yC, yHat, FFTW_FORWARD,
	FFTW_ESTIMATE);
	fPlanZ = fftw_plan_dft_1d(nDFT, zHat, zC, FFTW_BACKWARD,
	FFTW_ESTIMATE);
	Gamma = getGamma();
}

// delete allocated arrays
RBergomi::~RBergomi() {
	delete[] xC;
	delete[] xHat;
	delete[] yC;
	delete[] yHat;
	delete[] zC;
	delete[] zHat;
}

// Note that Wtilde[0] = 0!
void RBergomi::updateWtilde(Vector& Wtilde, const Vector& W1,
		const Vector& W1perp) {
	double s2H = sqrt(2.0 * H);
	double rhoH = s2H / (H + 0.5);
	Vector W1hat = linearComb(rhoH / s2H, W1, sqrt(1.0 - rhoH * rhoH) / s2H,
			W1perp);
	Vector Y2(N); // see R code
	// Convolve W1 and Gamma
	// Copy W1 and Gamma to complex arrays
	copyToComplex(W1, xC);
	copyToComplex(Gamma, yC);
	// DFT both
	fftw_execute(fPlanX); // DFT saved in xHat[0]
	fftw_execute(fPlanY); // DFT saved in yHat[0]
	// multiply xHat and yHat and save in zHat
	complexMult(xHat, yHat, zHat);
	// inverse DFT zHat
	fftw_execute(fPlanZ);
	// read out the real part, re-scale by 1/nDFT
	copyToReal(Y2, zC);
	scaleVector(Y2, 1.0 / nDFT);
	// Wtilde = (Y2 + W1hat) * sqrt(2*H) * dt^H ??
	Wtilde = linearComb(sqrt(2.0 * H) * pow(1.0 / N, H), Y2,
			sqrt(2.0 * H) * pow(1.0 / N, H), W1hat);
}

void RBergomi::scaleWtilde(Vector& WtildeScaled, const Vector& Wtilde) const {
	double T = N * ds;
	for (int i = 0; i < N; ++i)
		WtildeScaled[i] = pow(T, H) * Wtilde[i];
}

void RBergomi::scaleZ(Vector& ZScaled, const Vector& Z) const {
	for (int i = 0; i < N; ++i)
		ZScaled[i] = sqrt(ds) * Z[i];
}

void RBergomi::updateV(Vector& v, const Vector& WtildeScaled) const {
	v[0] = xi0(0.0);
	for (int i = 1; i < N; ++i)
		v[i] = xi0(i * ds)
				* exp(
						eta * WtildeScaled[i - 1]
								- 0.5 * eta * eta * pow(i * ds, 2 * H));
}

void RBergomi::updateZ(Vector& Z, const Vector& W1, const Vector& Wperp) const {
	Z = linearComb(rho, W1, rhobar, Wperp);
}

void RBergomi::updateS(Vector& S, const Vector& v,
		const Vector& ZScaled) const {
	Vector X(N + 1, log(S0));
	for (size_t i = 0; i < v.size(); ++i)
		X[i + 1] = X[i] + sqrt(v[i]) * ZScaled[i] - 0.5 * v[i] * ds;
	std::transform(X.begin(), X.end(), S.begin(),
			[](double x) -> double {return exp(x);});
}

Vector RBergomi::getGamma() {
	double alpha = H - 0.5;
	Vector Gamma(N);
	Gamma[0] = 0.0;
	for (int i = 1; i < N; ++i)
		Gamma[i] = (pow(i + 1.0, alpha + 1.0) - pow(i, alpha + 1.0))
				/ (alpha + 1.0);
	return Gamma;
}

void RBergomi::copyToComplex(const Vector& x, fftw_complex* xc) {
	for (size_t i = 0; i < x.size(); ++i) {
		xc[i][0] = x[i]; // real part
		xc[i][1] = 0.0; // imaginary part
	}
	// fill up with 0es
	for (size_t i = x.size(); i < nDFT; ++i) {
		xc[i][0] = 0.0; // real part
		xc[i][1] = 0.0; // imaginary part
	}
}

void RBergomi::copyToReal(Vector& x, const fftw_complex* xc) const {
	for (size_t i = 0; i < x.size(); ++i)
		x[i] = xc[i][0]; // real part
}

void RBergomi::complexMult(const fftw_complex* x, const fftw_complex* y,
		fftw_complex* z) {
	for (size_t i = 0; i < nDFT; ++i)
		fftw_c_mult(x[i], y[i], z[i]);
}

void RBergomi::fftw_c_mult(const fftw_complex a, const fftw_complex b,
		fftw_complex c) {
	c[0] = a[0] * b[0] - a[1] * b[1];
	c[1] = a[0] * b[1] + a[1] * b[0];
}

void RBergomi::generateSample(Vector& S, Vector& W1, Vector& Wperp) {
	Vector W1perp(N);
	Vector Wtilde(N);
	Vector WtildeScaled(N); // Wtilde scaled according to time
	Vector ZScaled(N);
	Vector v(N);
	Vector Z(N);
	double sds = sqrt(ds); // square root of time increment

	// generate the fundamental Gaussians
	ran(W1);
	ran(W1perp);
	ran(Wperp);

	updateWtilde(Wtilde, W1, W1perp);
	scaleWtilde(WtildeScaled, Wtilde);
	updateV(v, WtildeScaled);

	updateZ(Z, W1, Wperp);
	scaleZ(ZScaled, Z);
	updateS(S, v, ZScaled);

	scaleVector(W1, sds);
	scaleVector(Wperp, sds);
}

SampleCollection RBergomi::generate(long M) {
	Matrix S(M, Vector(N + 1));
	Matrix W1(M, Vector(N));
	Matrix Wperp(M, Vector(N));
	for (long m = 0; m < M; ++m)
		generateSample(S[m], W1[m], Wperp[m]);
	return SampleCollection { S, W1, Wperp };
}

SampleCollection RBergomi::operator ()(long M) {
	return generate(M);
}
