/*
 * RBergomiST.cpp
 *
 *  Created on: Nov 16, 2016
 *      Author: bayerc
 */

#include "RBergomi.h"

RBergomi::RBergomi() {
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
RBergomi::RBergomi(Vector x, Vector HIn, Vector e, Vector r, Vector t, Vector k,
		int NIn, long MIn,
		std::vector<uint64_t> seed = std::vector<uint64_t>(0)) {
	// Some safety tests/asserts?
	N = NIn;
	nDFT = 2 * N - 1;
	M = MIn;
	//numThreads = numThreadsIn;
	par = ParamTot(HIn, e, r, t, k, x);
	setGen(seed);
	dist = normDist(0.0, 1.0);
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
RBergomi::~RBergomi() {
	delete[] xC;
	delete[] xHat;
	delete[] yC;
	delete[] yHat;
	delete[] zC;
	delete[] zHat;
}

Result RBergomi::ComputePrice() {
	// The random vectors; the first 3 are independent, Z is a composite
	// Note that W1, W1perp, Wperp, Z correspond to UNNORMALIZED increments of Brownian motions,
	// i.e., are i.i.d. standard normal.
	Vector W1(N);
	Vector W1perp(N);
	Vector Wperp(N);
	Vector Wtilde(N);
	Vector WtildeScaled(N); // Wtilde scaled according to time
	Vector ZScaled(N);
	Vector v(N);
	Vector Z(N);
	double S; // maybe it is better to use a vector of S's corresponding to all different maturities!!
	// This would need a major re-organization of the code, including ParamTot...
	// The vectors of all combinations of parameter values and prices
	// vectors of prices and variances
	Vector price(par.size(), 0.0);
	Vector var(par.size(), 0.0);
	// other parameters used
	double dt; // time increment
	double sdt; // square root of time increment

	// The big loop which needs to be parallelized in future
	for (int m = 0; m < M; ++m) {
		// generate the fundamental Gaussians
		genGaussian(W1);
		genGaussian(W1perp);
		genGaussian(Wperp);

		// now iterate through all parameters
		for (long i = 0; i < par.size(); ++i) {
			updateWtilde(Wtilde, W1, W1perp, par.H(i));
			scaleWtilde(WtildeScaled, Wtilde, par.T(i), par.H(i));
			dt = par.T(i) / N;
			sdt = sqrt(dt);
			updateV(v, WtildeScaled, par.xi(i), par.H(i), par.eta(i), dt);

			updateZ(Z, W1, Wperp, par.rho(i));
			scaleZ(ZScaled, Z, sdt);
			S = updateS(v, ZScaled, dt);

			// now compute the payoff
			//double payoff = par.K(i) > 1.0 ? posPart(par.K(i) - S) : posPart(S - par.K(i));
			double payoff = posPart(S - par.K(i)); // call option
			price[i] += payoff;
			var[i] += payoff * payoff;
		}
	}

	// compute mean and variance
	scaleVector(price, 1.0 / double(M));
	scaleVector(var, 1.0 / double(M)); // = E[X^2]
	var = linearComb(1.0, var, -1.0, squareVector(price)); // = empirical var of price
	Vector stat = rootVector(var);
	scaleVector(stat, 1.0 / sqrt(double(M)));
	Vector iv(par.size());
	Result res { price, iv, stat, par, N, M, 0.0 };

	return res;
}

Result RBergomi::ComputeIV() {
	Result res = ComputePrice();
	// Compute the implied vol for each price saved in res
	for (long i = 0; i < res.par.size(); ++i) {
		res.iv[i] = IV_call(res.price[i], 1.0, res.par.K(i), res.par.T(i));
	}
	return res;
}

std::vector<Vector> RBergomi::ComputePayoffRT(const std::vector<Vector> & W1Arr,
		const std::vector<Vector> & W1perpArr) {
	// Check that input arrays have the right dimension
	if ((W1Arr[0].size() != N) || (W1perpArr[0].size() != N)) {
		std::cerr << "W1Arr[0].size() = " << W1Arr[0].size()
				<< ", W1perpArr[0].size() = " << W1perpArr[0].size()
				<< ", but N = " << N << '\n';
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
			payoffArr[m][i] = updatePayoff(i, Wtilde, WtildeScaled, W1, W1perp,
					v);
	}
	return payoffArr;
}

double RBergomi::ComputePayoffRT_single(const Vector & W1,
		const Vector & W1perp) {
	Vector Wtilde(N);
	Vector WtildeScaled(N); // Wtilde scaled according to time
	Vector v(N);
	double payoff = updatePayoff(0, Wtilde, WtildeScaled, W1, W1perp, v);
	return payoff;
}

Result RBergomi::ComputePriceRT() {
	// The random vectors; the first 3 are independent, Z is a composite
	// Note that W1, W1perp, Wperp, Z correspond to UNNORMALIZED increments of Brownian motions,
	// i.e., are i.i.d. standard normal.
	Vector W1(N);
	Vector W1perp(N);
	Vector Wtilde(N);
	Vector WtildeScaled(N); // Wtilde scaled according to time
	Vector v(N);
	Vector Z(N);
	// vectors of prices and variances
	Vector price(par.size(), 0.0);
	Vector var(par.size(), 0.0);

	// The big loop which needs to be parallelized in future
	for (int m = 0; m < M; ++m) {
		// generate the fundamental Gaussians
		genGaussian(W1);
		genGaussian(W1perp); // Wperp is  not needed!

		// now iterate through all parameters
		for (long i = 0; i < par.size(); ++i) {
			double payoff = updatePayoff(i, Wtilde, WtildeScaled, W1, W1perp,
					v);
			price[i] += payoff;
			var[i] += payoff * payoff;
		}
	}

	// compute mean and variance
	scaleVector(price, 1.0 / double(M));
	scaleVector(var, 1.0 / double(M)); // = E[X^2]
	var = linearComb(1.0, var, -1.0, squareVector(price)); // = empirical var of price
	Vector stat = rootVector(var);
	scaleVector(stat, 1.0 / sqrt(double(M)));
	Vector iv(par.size());
	Result res { price, iv, stat, par, N, M, 0.0 };

	return res;
}

Result RBergomi::ComputeIVRT() {
	Result res = ComputePriceRT();
	// Compute the implied vol for each price saved in res
	for (long i = 0; i < res.par.size(); ++i) {
		res.iv[i] = IV_call(res.price[i], 1.0, res.par.K(i), res.par.T(i));
	}
	return res;
}

void RBergomi::setGen(std::vector<uint64_t> seed) {
	std::seed_seq seeder(seed.begin(), seed.end());
	seeder.generate(seed.begin(), seed.end());
	gen.seed(seed[0]);
}

void RBergomi::genGaussian(Vector& X) {
	for (auto & x : X)
		x = dist(gen);
}

// Note that Wtilde[0] = 0!
void RBergomi::updateWtilde(Vector& Wtilde, const Vector& W1,
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

void RBergomi::scaleWtilde(Vector& WtildeScaled, const Vector& Wtilde, double T,
		double H) const {
	for (int i = 0; i < N; ++i)
		WtildeScaled[i] = pow(T, H) * Wtilde[i];
}

void RBergomi::scaleZ(Vector& ZScaled, const Vector& Z, double sdt) const {
	for (int i = 0; i < N; ++i)
		ZScaled[i] = sdt * Z[i];
}

void RBergomi::updateV(Vector& v, const Vector& WtildeScaled, double xi, double h,
		double e, double dt) const {
	v[0] = xi;
	for (int i = 1; i < N; ++i)
		v[i] = xi
				* exp(
						e * WtildeScaled[i - 1]
								- 0.5 * e * e * pow((i - 1) * dt, 2 * h));
}

void RBergomi::updateZ(Vector& Z, const Vector& W1, const Vector& Wperp,
		double r) const {
	Z = linearComb(r, W1, sqrt(1.0 - r * r), Wperp);
}

double RBergomi::updateS(const Vector& v, const Vector& ZScaled,
		double dt) const {
	double X = 0.0;
	for (size_t i = 0; i < v.size(); ++i)
		X += sqrt(v[i]) * ZScaled[i] - 0.5 * v[i] * dt;
	return exp(X); // Recall S_0 = 1.
}

void RBergomi::getGamma(Vector& Gamma, double H) const {
	double alpha = H - 0.5;
	Gamma[0] = 0.0;
	for (int i = 1; i < N; ++i)
		Gamma[i] = (pow(i + 1.0, alpha + 1.0) - pow(i, alpha + 1.0))
				/ (alpha + 1.0);
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

long RBergomi::getM() const {
	return M;
}

void RBergomi::setM(long m) {
	M = m;
}

int RBergomi::getN() const {
	return N;
}

void RBergomi::setN(int n) {
	N = n;
}

double RBergomi::intVdt(const Vector & v, double dt) const {
	return dt * std::accumulate(v.begin(), v.end(), 0.0);
}

double RBergomi::intRootVdW(const Vector & v, const Vector & W1,
		double sdt) const {
	double IsvdW = 0.0;
	for (size_t i = 0; i < v.size(); ++i)
		IsvdW += sqrt(v[i]) * sdt * W1[i];
	return IsvdW;
}

double RBergomi::updatePayoff(long i, Vector& Wtilde, Vector& WtildeScaled,
		const Vector& W1, const Vector& W1perp, Vector& v) {
	double dt = par.T(i) / N;
	double sdt = sqrt(dt);
	updateWtilde(Wtilde, W1, W1perp, par.H(i));
	scaleWtilde(WtildeScaled, Wtilde, par.T(i), par.H(i));
	updateV(v, WtildeScaled, par.xi(i), par.H(i), par.eta(i), dt);
	// now compute \int v_s ds, \int \sqrt{v_s} dW_s with W = W1
	double Ivdt = intVdt(v, dt);
	double IsvdW = intRootVdW(v, W1, sdt);
	// now compute the payoff by inserting properly into the BS formula
	double BS_vol = sqrt((1.0 - par.rho(i) * par.rho(i)) * Ivdt);
	double BS_spot = exp(
			-0.5 * par.rho(i) * par.rho(i) * Ivdt + par.rho(i) * IsvdW);
	return BS_call_price(BS_spot, par.K(i), 1.0, BS_vol);
}
