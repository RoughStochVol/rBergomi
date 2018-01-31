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
RBergomiST::~RBergomiST() {
	delete[] xC;
	delete[] xHat;
	delete[] yC;
	delete[] yHat;
	delete[] zC;
	delete[] zHat;
}

Result RBergomiST::ComputePrice() {
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
			// Note that each of the changes here forces all subsequent updates!
			// check if H has changed. If so, Wtilde needs to be updated (and, hence, everything else)
			bool update = par.HTrigger(i);
			if (update)
				updateWtilde(Wtilde, W1, W1perp, par.H(i));
			// check if T has changed. If so, Wtilde and the time increment need re-scaling
			update = update || par.TTrigger(i);
			if (update) {
				scaleWtilde(WtildeScaled, Wtilde, par.T(i), par.H(i));
				dt = par.T(i) / N;
				sdt = sqrt(dt);
			}
			// check if eta has changed. If so, v needs to be updated
			update = update || par.etaTrigger(i);
			if (update)
				updateV(v, WtildeScaled, par.H(i), par.eta(i), dt);
			// if rho has changed, then S needs to be re-computed
			update = update || par.rhoTrigger(i);
			if (update) {
				updateZ(Z, W1, Wperp, par.rho(i));
				scaleZ(ZScaled, Z, sdt);
				S = updateS(v, ZScaled, dt);
			}
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
	Result res { price, iv, par, stat, N, M, 1, 0.0 };

	return res;
}

Result RBergomiST::ComputeIV() {
	Result res = ComputePrice();
	// Compute the implied vol for each price saved in res
	for(long i=0; i<res.par.size(); ++i){
		res.iv[i] = IV_call(res.price[i], 1.0, res.par.K(i), res.par.T(i));
	}
	return res;
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


Result RBergomiST::ComputePriceRT() {
	// The random vectors; the first 3 are independent, Z is a composite
	// Note that W1, W1perp, Wperp, Z correspond to UNNORMALIZED increments of Brownian motions,
	// i.e., are i.i.d. standard normal.
	Vector W1(N);
	Vector W1perp(N);
	Vector Wtilde(N);
	Vector WtildeScaled(N); // Wtilde scaled according to time
	Vector v(N);
	Vector Z(N);
	double Ivdt, IsvdW; // \int v_s ds, \int \sqrt{v_s} dW_s; Here, W = W1
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
		genGaussian(W1perp); // Wperp is  not needed!

		// now iterate through all parameters
		for (long i = 0; i < par.size(); ++i) {
			double payoff = updatePayoff(i, Wtilde, WtildeScaled, W1, W1perp, v);
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
	Result res { price, iv, par, stat, N, M, 1, 0.0 };

	return res;
}

Result RBergomiST::ComputeIVRT() {
	Result res = ComputePriceRT();
	// Compute the implied vol for each price saved in res
	for(long i=0; i<res.par.size(); ++i){
		res.iv[i] = IV_call(res.price[i], 1.0, res.par.K(i), res.par.T(i));
	}
	return res;
}

void RBergomiST::setGen(std::vector<uint64_t> seed) {
	std::seed_seq seeder(seed.begin(), seed.end());
	seeder.generate(seed.begin(), seed.end());
	//for (int i = 0; i < numThreads; ++i) {
		//seeder = std::seed_seq{seed[i]};
		gen.seed(seed[0]);
		//gen[i] = MTGenerator(seeder);
	//}
}

void RBergomiST::genGaussian(Vector& X) {
	for(auto & x : X)
		x = dist(gen); // replace 0 by omp_get_thread_num() for multi-threaded code
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

void RBergomiST::updateZ(Vector& Z, const Vector& W1, const Vector& Wperp,
		double r) const {
	Z = linearComb(r, W1, sqrt(1.0 - r * r), Wperp);
}

double RBergomiST::updateS(const Vector& v, const Vector& ZScaled, double dt) const {
	double X = 0.0;
	for (size_t i = 0; i < v.size(); ++i)
		X += sqrt(v[i]) * ZScaled[i] - 0.5 * v[i] * dt;
	return exp(X); // Recall S_0 = 1.
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

int RBergomiST::getNumThreads() const {
	return 0;
}

double RBergomiST::getXi() const {
	return par.xi();
}

void RBergomiST::testScaleWtilde() {
	Vector W1(N);
	Vector Wtilde(N);
	Vector W1perp(N);
	Vector Wperp(N);
	Vector WtildeScaled(N); // Wtilde scaled according to time
	Vector v(N);
	Vector Z(N);
	Vector ZScaled(N);
	double S;
	double H_scalar = par.H(0);
	double T_scalar = par.T(0);
	double eta_scalar = par.eta(0);
	double rho_scalar = par.rho(0);

	std::cout << "H = " << H_scalar << ", T = " << T_scalar << ", eta = "
			<< eta_scalar << ", rho = " << rho_scalar << std::endl;

	// names of files for W1, Wtilde, W1perp, WtildeScaled
	std::string f_W1 = "./W1.txt";
	std::string f_Wtilde = "./Wtilde.txt";
	std::string f_W1perp = "./W1perp.txt";
	std::string f_Wperp = "./Wperp.txt";
	std::string f_WtildeScaled = "./WtildeScaled.txt";
	std::string f_v = "./v.txt";
	std::string f_Z = "./Z.txt";
	std::string f_S = "./S.txt";
	std::fstream file;
	// make sure that files are empty
	file.open(f_W1.c_str(), std::fstream::out | std::fstream::trunc);
	file << "";
	file.close();
	file.open(f_Wtilde.c_str(), std::fstream::out | std::fstream::trunc);
	file << "";
	file.close();
	file.open(f_W1perp.c_str(), std::fstream::out | std::fstream::trunc);
	file << "";
	file.close();
	file.open(f_WtildeScaled.c_str(), std::fstream::out | std::fstream::trunc);
	file << "";
	file.close();
	file.open(f_Wperp.c_str(), std::fstream::out | std::fstream::trunc);
	file << "";
	file.close();
	file.open(f_v.c_str(), std::fstream::out | std::fstream::trunc);
	file << "";
	file.close();
	file.open(f_Z.c_str(), std::fstream::out | std::fstream::trunc);
	file << "";
	file.close();
	file.open(f_S.c_str(), std::fstream::out | std::fstream::trunc);
	file << "";
	file.close();
	// The big loop which needs to be parallelized in future
	for (int m = 0; m < M; ++m) {
		// generate the fundamental Gaussians
		genGaussian(W1);
		genGaussian(W1perp);
		genGaussian(Wperp);
		updateWtilde(Wtilde, W1, W1perp, H_scalar);
		scaleWtilde(WtildeScaled, Wtilde, T_scalar, H_scalar);
		updateV(v, WtildeScaled, H_scalar, eta_scalar, T_scalar / N);
		updateZ(Z, W1, Wperp, rho_scalar);
		scaleZ(ZScaled, Z, sqrt(T_scalar / N));
		S = updateS(v, ZScaled, T_scalar / N);
		file.open(f_W1.c_str(), std::fstream::out | std::fstream::app);
		std::copy(W1.begin(), W1.end(),
				std::ostream_iterator<double>(file, " "));
		file << "\n";
		file.close();
		file.open(f_Wtilde.c_str(), std::fstream::out | std::fstream::app);
		std::copy(Wtilde.begin(), Wtilde.end(),
				std::ostream_iterator<double>(file, " "));
		file << "\n";
		file.close();
		file.open(f_W1perp.c_str(), std::fstream::out | std::fstream::app);
		std::copy(W1perp.begin(), W1perp.end(),
				std::ostream_iterator<double>(file, " "));
		file << "\n";
		file.close();
		file.open(f_WtildeScaled.c_str(), std::fstream::out | std::fstream::app);
		std::copy(WtildeScaled.begin(), WtildeScaled.end(),
				std::ostream_iterator<double>(file, " "));
		file << "\n";
		file.close();
		file.open(f_Wperp.c_str(), std::fstream::out | std::fstream::app);
		std::copy(Wperp.begin(), Wperp.end(),
				std::ostream_iterator<double>(file, " "));
		file << "\n";
		file.close();
		file.open(f_v.c_str(), std::fstream::out | std::fstream::app);
		std::copy(v.begin(), v.end(), std::ostream_iterator<double>(file, " "));
		file << "\n";
		file.close();
		file.open(f_Z.c_str(), std::fstream::out | std::fstream::app);
		std::copy(Z.begin(), Z.end(), std::ostream_iterator<double>(file, " "));
		file << "\n";
		file.close();
		file.open(f_S.c_str(), std::fstream::out | std::fstream::app);
		file << S;
		file << "\n";
		file.close();
	}
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

void RBergomiST::testConvolve() {
	// Convolve W1 and Gamma
	Vector Gamma(N);
	getGamma(Gamma, par.H(0));
	// Choose special W1:
	Vector W1(N);
	for (int i = 0; i < N; ++i)
		W1[i] = i;
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
	Vector Y2(N);
	copyToReal(Y2, zC);
	scaleVector(Y2, 1.0 / nDFT);
	std::cout << "Gamma = " << Gamma << std::endl;
	std::cout << "W1 = " << W1 << std::endl;
	std::cout << "Gamma * W1 = " << Y2 << std::endl;
}

void RBergomiST::testWtilde() {
	Vector W1 = { 0.0688183, 1.07636, -1.84502, -0.561138, 0.932903, 0.926146,
			-1.99161, 0.732015, -0.107891, -0.367463, 0.867313, 0.198056,
			-1.50203, -1.62697, 0.192327, -0.164254, -0.063662, -0.581115,
			-0.50949, 0.716903, -1.63944, -0.637341, 1.36474, 1.91016,
			0.912126, -0.0138699, 0.299401, -0.532777, 0.0128634, -0.334329,
			1.59338, 0.857188, -0.257853, 0.610727, -0.0761325, 0.657633,
			0.389386, -1.59859, 1.07015, 0.539286, -1.06533, -1.11202, 1.79673,
			-0.59522, 0.543368, 0.497823, -0.148346, 0.167303, 0.107486,
			-1.03891, -0.627193, 0.439227, 0.0449506, 0.0719449, -0.286124,
			0.338444, 0.21025, -0.953058, 0.420749, 0.129486, -0.568006,
			1.09033, 1.09979, -0.545845, 1.01928, 1.00078, -0.38289, 0.80968,
			1.01335, 0.0791322, -0.912119, 2.62517, -0.928324, 0.478441,
			-1.20624, -0.862473, -0.488587, -1.25914, -0.723156, 2.08981,
			0.591217, 1.29106, 0.810396, 0.133946, -0.760356, 0.463027,
			0.845751, 1.96196, -0.109956, -1.08388, 0.584277, 0.342109,
			0.288257, -1.69693, -0.858363, 0.399129, 0.0614502, -1.15613,
			0.691034, 0.593261, -1.07506, 1.34609, 0.660385, 0.698617, 2.03886,
			-1.39825, 0.935031, -0.889707, -1.80273, 0.523153, 0.165514,
			1.45271, 0.799404, 1.28545, 1.12914, -1.41436, 0.00587025,
			0.0437768, -0.422505, 1.2156, -0.687495, 0.966272, -0.0819264,
			1.19317, -0.881422, 0.698812, 2.3643, 0.569367, 0.0503847,
			-0.191175, -0.370943, -0.523355, 1.78402, -0.917506, 0.639609,
			-0.133216, -0.826925, 0.379126, 1.16025, -0.883843, -0.212405,
			0.0113117, 3.39823, 0.0168886, -0.604783, 0.395182, 0.0828049,
			1.74513, 0.187216, -0.567286, 0.445468, -1.14508, 1.90739, -1.0979,
			-1.28313, -0.785638, -1.14034, 0.0073991, -0.621622, -0.103495,
			-1.24885, -1.5407, -1.1973, -0.270162, 0.480444, -0.224886,
			-1.23382, -1.46449, 2.87728, -0.833233, -0.870657, -0.549865,
			0.741537, 0.853076, -0.8096, 0.67628, -0.217569, -0.730999,
			0.95267, 0.71366, 1.05699, 2.7518, -0.666633, -0.136098, 0.249606,
			-0.151262, -0.0319986, 0.476672, 1.05305, -1.05503, 0.106315,
			0.363139, -1.3216, -1.40817, 0.512966, 0.514397, -0.23803,
			-0.495276, -2.03592, 1.32942 };
	//(N);
	Vector W1perp = { 0.0688183, 1.07636, -1.84502, -0.561138, 0.932903,
			0.926146, -1.99161, 0.732015, -0.107891, -0.367463, 0.867313,
			0.198056, -1.50203, -1.62697, 0.192327, -0.164254, -0.063662,
			-0.581115, -0.50949, 0.716903, -1.63944, -0.637341, 1.36474,
			1.91016, 0.912126, -0.0138699, 0.299401, -0.532777, 0.0128634,
			-0.334329, 1.59338, 0.857188, -0.257853, 0.610727, -0.0761325,
			0.657633, 0.389386, -1.59859, 1.07015, 0.539286, -1.06533,
			-1.11202, 1.79673, -0.59522, 0.543368, 0.497823, -0.148346,
			0.167303, 0.107486, -1.03891, -0.627193, 0.439227, 0.0449506,
			0.0719449, -0.286124, 0.338444, 0.21025, -0.953058, 0.420749,
			0.129486, -0.568006, 1.09033, 1.09979, -0.545845, 1.01928, 1.00078,
			-0.38289, 0.80968, 1.01335, 0.0791322, -0.912119, 2.62517,
			-0.928324, 0.478441, -1.20624, -0.862473, -0.488587, -1.25914,
			-0.723156, 2.08981, 0.591217, 1.29106, 0.810396, 0.133946,
			-0.760356, 0.463027, 0.845751, 1.96196, -0.109956, -1.08388,
			0.584277, 0.342109, 0.288257, -1.69693, -0.858363, 0.399129,
			0.0614502, -1.15613, 0.691034, 0.593261, -1.07506, 1.34609,
			0.660385, 0.698617, 2.03886, -1.39825, 0.935031, -0.889707,
			-1.80273, 0.523153, 0.165514, 1.45271, 0.799404, 1.28545, 1.12914,
			-1.41436, 0.00587025, 0.0437768, -0.422505, 1.2156, -0.687495,
			0.966272, -0.0819264, 1.19317, -0.881422, 0.698812, 2.3643,
			0.569367, 0.0503847, -0.191175, -0.370943, -0.523355, 1.78402,
			-0.917506, 0.639609, -0.133216, -0.826925, 0.379126, 1.16025,
			-0.883843, -0.212405, 0.0113117, 3.39823, 0.0168886, -0.604783,
			0.395182, 0.0828049, 1.74513, 0.187216, -0.567286, 0.445468,
			-1.14508, 1.90739, -1.0979, -1.28313, -0.785638, -1.14034,
			0.0073991, -0.621622, -0.103495, -1.24885, -1.5407, -1.1973,
			-0.270162, 0.480444, -0.224886, -1.23382, -1.46449, 2.87728,
			-0.833233, -0.870657, -0.549865, 0.741537, 0.853076, -0.8096,
			0.67628, -0.217569, -0.730999, 0.95267, 0.71366, 1.05699, 2.7518,
			-0.666633, -0.136098, 0.249606, -0.151262, -0.0319986, 0.476672,
			1.05305, -1.05503, 0.106315, 0.363139, -1.3216, -1.40817, 0.512966,
			0.514397, -0.23803, -0.495276, -2.03592, 1.32942 };//(N);
	Vector Wtilde(N);
	Vector WtildeScaled(N); // Wtilde scaled according to time
	double H_scalar = par.H(0);
	updateWtilde(Wtilde, W1, W1perp, H_scalar);
	std::cout << "W1 = " << W1 << std::endl;
	std::cout << "W1perp = " << W1perp << std::endl;
	std::cout << "Wtilde = " << Wtilde << std::endl;
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
