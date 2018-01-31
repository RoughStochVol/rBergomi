/*
 * RBergomiST.cpp
 *
 *  Created on: Nov 16, 2016
 *      Author: bayerc
 */

#include "RBergomi.h"

RBergomiST::RBergomiST() {
	N = 0;
	M = 0;
	nDFT = 0;
	par = ParamTot();
	xC = new fftw_complex[0];
	xHat = new fftw_complex[0];
	yC = new fftw_complex[0];
	yHat = new fftw_complex[0];
	zC = new fftw_complex[0];
	zHat = new fftw_complex[0];
}

// seed is optional
RBergomiST::RBergomiST(double x, Vector HIn, Vector e, Vector r, Vector t,
		Vector k, int NIn, long MIn,
		std::vector<uint64_t> seed = std::vector<uint64_t>(0)) {
	// Some safety tests/asserts?
	N = NIn;
	nDFT = 2 * N - 1;
	M = MIn;
	//numThreads = numThreadsIn;
	par = ParamTot(HIn, e, r, t, k, x);
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
}

// delete allocated arrays
RBergomiST::~RBergomiST() {
	delete[] xC;
	delete[] xHat;
	delete[] yC;
	delete[] yHat;
	delete[] zC;
	delete[] zHat;
}

// Note that Wtilde[0] = 0!
void RBergomiST::updateWtilde(Vector& Wtilde, const Vector& W1,
		const Vector& W1perp, double H) {
	Vector Gamma(N);
	getGamma(Gamma, H);
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

void RBergomiST::scaleWtilde(Vector& WtildeScaled, const Vector& Wtilde,
		double T, double H) const {
	for (int i = 0; i < N; ++i)
		WtildeScaled[i] = pow(T, H) * Wtilde[i];
}

void RBergomiST::scaleZ(Vector& ZScaled, const Vector& Z, double sdt) const {
	for (int i = 0; i < N; ++i)
		ZScaled[i] = sdt * Z[i];
}

void RBergomiST::updateV(Vector& v, const Vector& WtildeScaled, double h,
		double e, double dt) const {
	v[0] = par.xi();
	for (int i = 1; i < N; ++i)
		v[i] = par.xi() * exp(
				e * WtildeScaled[i - 1] - 0.5 * e * e
						* pow((i - 1) * dt, 2 * h));
}

void RBergomiST::getGamma(Vector& Gamma, double H) const {
	double alpha = H - 0.5;
	Gamma[0] = 0.0;
	for (int i = 1; i < N; ++i)
		Gamma[i] = (pow(i + 1.0, alpha + 1.0) - pow(i, alpha + 1.0)) / (alpha
				+ 1.0);
}

void RBergomiST::copyToComplex(const Vector& x, fftw_complex* xc) {
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

void RBergomiST::copyToReal(Vector& x, const fftw_complex* xc) const {
	for (size_t i = 0; i < x.size(); ++i)
		x[i] = xc[i][0]; // real part
}

void RBergomiST::complexMult(const fftw_complex* x, const fftw_complex* y,
		fftw_complex* z) {
	for (size_t i = 0; i < nDFT; ++i)
		fftw_c_mult(x[i], y[i], z[i]);
}

void RBergomiST::fftw_c_mult(const fftw_complex a, const fftw_complex b,
		fftw_complex c) {
	c[0] = a[0] * b[0] - a[1] * b[1];
	c[1] = a[0] * b[1] + a[1] * b[0];
}

long RBergomiST::getM() const {
	return M;
}

void RBergomiST::setM(long m) {
	M = m;
}

int RBergomiST::getN() const {
	return N;
}

void RBergomiST::setN(int n) {
	N = n;
}

double RBergomiST::getXi() const {
	return par.xi();
}

double RBergomiST::intVdt(const Vector & v, double dt) const
{
	return dt * std::accumulate(v.begin(), v.end(), 0.0);
}

double RBergomiST::intRootVdW(const Vector & v, const Vector & W1, double sdt) const
{
	double IsvdW = 0.0;
	for(size_t i=0; i<v.size(); ++i)
		IsvdW += sqrt(v[i]) * sdt * W1[i];
	return IsvdW;
}

double RBergomiST::updatePayoff(long i, Vector& Wtilde,
		Vector& WtildeScaled, const Vector& W1, const Vector& W1perp, Vector& v){
	double dt = par.T(i) / N;
	double sdt = sqrt(dt);
	bool update = par.HTrigger(i);
	if (update)
		updateWtilde(Wtilde, W1, W1perp, par.H(i));
	// check if T has changed. If so, Wtilde and the time increment need re-scaling
	update = update || par.TTrigger(i);
	if (update) {
		scaleWtilde(WtildeScaled, Wtilde, par.T(i), par.H(i));
	}
	// check if eta has changed. If so, v needs to be updated
	update = update || par.etaTrigger(i);
	if (update)
		updateV(v, WtildeScaled, par.H(i), par.eta(i), dt);
	// now compute \int v_s ds, \int \sqrt{v_s} dW_s with W = W1
	double Ivdt = intVdt(v, dt);
	double IsvdW = intRootVdW(v, W1, sdt);
	// now compute the payoff by inserting properly into the BS formula
	double BS_vol = sqrt((1.0 - par.rho(i)*par.rho(i)) * Ivdt);
	double BS_spot = exp( - 0.5 * par.rho(i)*par.rho(i) * Ivdt + par.rho(i) * IsvdW );
	return BS_call_price(BS_spot, par.K(i), 1.0, BS_vol);
}

std::vector<Vector> RBergomiST::ComputePayoffRT(const std::vector<Vector> & W1Arr,
		const std::vector<Vector> & W1perpArr){
	// Check that input arrays have the right dimension
	if((W1Arr[0].size() != N) || (W1perpArr[0].size() != N)){
		std::cerr << "W1Arr[0].size() = " << W1Arr[0].size() << ", W1perpArr[0].size() = "
				<< W1perpArr[0].size() << ", but N = " << N << '\n';
		exit(17);
	}
	// The random vectors; the first 3 are independent, Z is a composite
	// Note that W1, W1perp, Wperp, Z correspond to UNNORMALIZED increments of Brownian motions,
	// i.e., are i.i.d. standard normal.
	Vector W1(N);
	Vector W1perp(N);
	Vector Wtilde(N);
	Vector WtildeScaled(N); // Wtilde scaled according to time
	Vector v(N);
	// Construct array for payoffs
	std::vector<Vector> payoffArr(W1Arr.size(), Vector(par.size()));

	// The big loop which needs to be parallelized in future
	for (int m = 0; m < W1Arr.size(); ++m) {
		// generate the fundamental Gaussians
		W1 = W1Arr[m];
		W1perp = W1perpArr[m];

		// now iterate through all parameters
		for (long i = 0; i < par.size(); ++i)
			payoffArr[m][i] = updatePayoff(i, Wtilde, WtildeScaled, W1, W1perp, v);
	}
	return payoffArr;
}

double RBergomiST::ComputePayoffRT_single(const Vector & W1, const Vector & W1perp){
	Vector Wtilde(N);
	Vector WtildeScaled(N); // Wtilde scaled according to time
	Vector v(N);
	double payoff = updatePayoff(0, Wtilde, WtildeScaled, W1, W1perp, v);
	return payoff;
}

