/*
 * fftOMP.cpp
 *
 *  Created on: 8 Dec 2017
 *      Author: bayerc
 */


#include "fftOMP.h"

void copyToComplexNew(const Vector& x, fftw_complex* xc, size_t nDFT) {
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

void copyToRealNew(Vector& x, const fftw_complex* xc) {
	for (size_t i = 0; i < x.size(); ++i)
		x[i] = xc[i][0]; // real part
}

void complexMultNew(const fftw_complex* x, const fftw_complex* y,
		fftw_complex* z, size_t nDFT) {
	for (size_t i = 0; i < nDFT; ++i)
		fftw_c_multNew(x[i], y[i], z[i]);
}

void fftw_c_multNew(const fftw_complex a, const fftw_complex b,
		fftw_complex c) {
	c[0] = a[0] * b[0] - a[1] * b[1];
	c[1] = a[0] * b[1] + a[1] * b[0];
}

Array convolveST(const Array& X, const Array& Y){
	size_t M = X.size();
	size_t N = X[0].size();
	size_t nDFT = 2 * N - 1;
	fftw_complex * xC = new fftw_complex[nDFT];
	fftw_complex * xHat = new fftw_complex[nDFT];
	fftw_complex * yC = new fftw_complex[nDFT];
	fftw_complex * yHat = new fftw_complex[nDFT];
	fftw_complex * zC = new fftw_complex[nDFT];
	fftw_complex * zHat = new fftw_complex[nDFT];
	fftw_plan fPlanX = fftw_plan_dft_1d(nDFT, xC, xHat, FFTW_FORWARD,
			FFTW_ESTIMATE);
	fftw_plan fPlanY = fftw_plan_dft_1d(nDFT, yC, yHat, FFTW_FORWARD,
			FFTW_ESTIMATE);
	fftw_plan fPlanZ = fftw_plan_dft_1d(nDFT, zHat, zC, FFTW_BACKWARD,
			FFTW_ESTIMATE);
	Array Z(M, Vector(N));

	for(size_t m = 0; m < M; ++m){
		copyToComplexNew(X[m], xC, nDFT);
		copyToComplexNew(Y[m], yC, nDFT);
		// DFT both
		fftw_execute(fPlanX); // DFT saved in xHat[0]
		fftw_execute(fPlanY); // DFT saved in yHat[0]
		// multiply xHat and yHat and save in zHat
		complexMultNew(xHat, yHat, zHat, nDFT);
		// inverse DFT zHat
		fftw_execute(fPlanZ);
		// read out the real part, re-scale by 1/nDFT
		copyToRealNew(Z[m], zC);
		scaleVector(Z[m], 1.0 / nDFT);
	}
	return Z;
}

// This version seeems to work very well!
Array convolveMT(const Array& X, const Array& Y, int numThreads){
	size_t M = X.size();
	size_t N = X[0].size();
	size_t nDFT = 2 * N - 1;
	fftData fft(nDFT, numThreads);
	Array Z(M, Vector(N));

	// Enforce that OMP use numThreads threads
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(numThreads);

	#pragma omp parallel for shared(X, Y, M, N, nDFT, fft, Z)
	for(size_t m = 0; m < M; ++m)
	{
		copyToComplexNew(X[m], fft.xC[omp_get_thread_num()], nDFT);
		copyToComplexNew(Y[m], fft.yC[omp_get_thread_num()], nDFT);
		// DFT both
		fftw_execute(fft.fPlanX[omp_get_thread_num()]); // DFT saved in xHat[0]
		fftw_execute(fft.fPlanY[omp_get_thread_num()]); // DFT saved in yHat[0]
		// multiply xHat and yHat and save in zHat
		complexMultNew(fft.xHat[omp_get_thread_num()], fft.yHat[omp_get_thread_num()], fft.zHat[omp_get_thread_num()], nDFT);
		// inverse DFT zHat
		fftw_execute(fft.fPlanZ[omp_get_thread_num()]);
		// read out the real part, re-scale by 1/nDFT
		copyToRealNew(Z[m], fft.zC[omp_get_thread_num()]);
		scaleVector(Z[m], 1.0 / nDFT);
	}
	return Z;
}

Array sampleArrayCos(size_t M, size_t N){
	Array Z(M, Vector(N));
	for(size_t m=0; m<M; ++m){
		for(size_t n=0; n<N; ++n)
			Z[m][n] = cos(n*m);
	}
	return Z;
}

Array sampleArraySin(size_t M, size_t N){
	Array Z(M, Vector(N));
	for(size_t m=0; m<M; ++m){
		for(size_t n=0; n<N; ++n)
			Z[m][n] = sin(n*m);
	}
	return Z;
}

std::vector<bool> compareArrays(const Array& X, const Array& Y){
	std::vector<bool> comp(X.size());
	for(size_t i = 0; i < comp.size(); ++i)
		comp[i] = std::equal(X[i].begin(), X[i].end(), Y[i].begin(), pred);
	return comp;
}

bool pred(double x, double y){
	return fabs(x - y) < 0.0000001;
}
