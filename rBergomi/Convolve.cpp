/*
 * Convolve.cpp
 *
 *  Created on: 31 Jan 2018
 *      Author: bayerc
 */

#include "Convolve.h"

Convolve::Convolve(int n) {
	N = n;
}

ConvolveFFTW::ConvolveFFTW() : Convolve() {
	nDFT = 0;
	xC = 0;
	xHat = 0;
	yC = 0;
	yHat = 0;
	zC = 0;
	zHat = 0;
	fPlanX = 0;
	fPlanY = 0;
	fPlanZ = 0;
}

ConvolveFFTW::ConvolveFFTW(int N) : Convolve(N) {
	nDFT = 2 * N - 1;
	xC = fftw_alloc_complex(nDFT);
	xHat = fftw_alloc_complex(nDFT);
	yC = fftw_alloc_complex(nDFT);
	yHat = fftw_alloc_complex(nDFT);
	zC = fftw_alloc_complex(nDFT);
	zHat = fftw_alloc_complex(nDFT);
	fPlanX = fftw_plan_dft_1d(nDFT, xC, xHat, FFTW_FORWARD,
					FFTW_PATIENT);
	fPlanY = fftw_plan_dft_1d(nDFT, yC, yHat, FFTW_FORWARD,
						FFTW_PATIENT);
	fPlanZ = fftw_plan_dft_1d(nDFT, zHat, zC, FFTW_BACKWARD,
			FFTW_PATIENT);
}

ConvolveFFTW::ConvolveFFTW(ConvolveFFTW& conv) {
	N = conv.N;
	nDFT = conv.nDFT;
	xC = fftw_alloc_complex(nDFT);
	xHat = fftw_alloc_complex(nDFT);
	yC = fftw_alloc_complex(nDFT);
	yHat = fftw_alloc_complex(nDFT);
	zC = fftw_alloc_complex(nDFT);
	zHat = fftw_alloc_complex(nDFT);
	fPlanX = fftw_plan_dft_1d(nDFT, xC, xHat, FFTW_FORWARD,
					FFTW_PATIENT);
	fPlanY = fftw_plan_dft_1d(nDFT, yC, yHat, FFTW_FORWARD,
						FFTW_PATIENT);
	fPlanZ = fftw_plan_dft_1d(nDFT, zHat, zC, FFTW_BACKWARD,
			FFTW_PATIENT);
}

ConvolveFFTW& ConvolveFFTW::operator =(const ConvolveFFTW& conv) {
	this->N = conv.N;
	this->nDFT = conv.nDFT;
	this->xC = fftw_alloc_complex(nDFT);
	this->xHat = fftw_alloc_complex(nDFT);
	this->yC = fftw_alloc_complex(nDFT);
	this->yHat = fftw_alloc_complex(nDFT);
	this->zC = fftw_alloc_complex(nDFT);
	this->zHat = fftw_alloc_complex(nDFT);
	this->fPlanX = fftw_plan_dft_1d(nDFT, xC, xHat, FFTW_FORWARD,
					FFTW_PATIENT);
	this->fPlanY = fftw_plan_dft_1d(nDFT, yC, yHat, FFTW_FORWARD,
						FFTW_PATIENT);
	this->fPlanZ = fftw_plan_dft_1d(nDFT, zHat, zC, FFTW_BACKWARD,
			FFTW_PATIENT);
	return *this;
}

ConvolveFFTW::~ConvolveFFTW() {
	fftw_free(xC);
	fftw_free(xHat);
	fftw_free(yC);
	fftw_free(yHat);
	fftw_free(zC);
	fftw_free(zHat);
	fftw_destroy_plan(fPlanX);
	fftw_destroy_plan(fPlanY);
	fftw_destroy_plan(fPlanZ);
}

void ConvolveFFTW::copyToComplex(const Vector & x, fftw_complex *xc) {
	for (int i = 0; i < N; ++i) {
		xc[i][0] = x[i]; // real part
		xc[i][1] = 0.0; // imaginary part
	}
	// fill up with 0es
	for (int i = N; i < nDFT; ++i) {
		xc[i][0] = 0.0; // real part
		xc[i][1] = 0.0; // imaginary part
	}
}

void ConvolveFFTW::fftw_c_mult(const fftw_complex a, const fftw_complex b, fftw_complex c) {
	c[0] = a[0] * b[0] - a[1] * b[1];
	c[1] = a[0] * b[1] + a[1] * b[0];
}

void ConvolveFFTW::complexMult(const fftw_complex *x, const fftw_complex *y,
		fftw_complex *z) {
	for (int i = 0; i < nDFT; ++i)
		fftw_c_mult(x[i], y[i], z[i]);
}

void ConvolveFFTW::copyToReal(Vector & x, const fftw_complex *xc) {
	for (int i = 0; i < N; ++i)
		x[i] = xc[i][0]; // real part
}

void ConvolveFFTW::run(const Vector& x, const Vector& y, Vector& z) {
	copyToComplex(x, xC);
	copyToComplex(y, yC);
	// DFT both
	fftw_execute(fPlanX); // DFT saved in xHat
	fftw_execute(fPlanY); // DFT saved in yHat[0]
	// multiply xHat and yHat and save in zHat
	complexMult(xHat, yHat, zHat);
	// inverse DFT zHat
	fftw_execute(fPlanZ);
	// read out the real part, re-scale by 1/nDFT
	copyToReal(z, zC);
	scaleVector(z, 1.0 / nDFT);
}

void ConvolveFFTW::operator ()(const Vector& x, const Vector& y, Vector& z) {
	this->run(x, y, z);
}

ConvolveGSL::ConvolveGSL() : Convolve() {
	N = 0;
	nDFT = 0;
	wavetable = 0;
	workspace = 0;
}

ConvolveGSL::ConvolveGSL(int N) : Convolve(N) {
	nDFT = 2 * N - 1;
	xC = Vector(2*nDFT);
	yC = Vector(2*nDFT);
	zC = Vector(2*nDFT);
	wavetable = gsl_fft_complex_wavetable_alloc(nDFT);
	workspace = gsl_fft_complex_workspace_alloc(nDFT);
}

ConvolveGSL::ConvolveGSL(ConvolveGSL& conv){
	N = conv.N;
	nDFT = conv.nDFT;
	xC = Vector(2*nDFT);
	yC = Vector(2*nDFT);
	zC = Vector(2*nDFT);
	wavetable = gsl_fft_complex_wavetable_alloc(nDFT);
	workspace = gsl_fft_complex_workspace_alloc(nDFT);
}

ConvolveGSL& ConvolveGSL::operator=(const ConvolveGSL& conv){
	this->N = conv.N;
	this->nDFT = conv.nDFT;
	this->xC = conv.xC;
	this->yC = conv.yC;
	this->zC = conv.zC;
	this->wavetable = gsl_fft_complex_wavetable_alloc(nDFT);
	this->workspace = gsl_fft_complex_workspace_alloc(nDFT);
	return *this;
}

ConvolveGSL::~ConvolveGSL(){
	gsl_fft_complex_wavetable_free(wavetable);
	gsl_fft_complex_workspace_free(workspace);
}

void ConvolveGSL::copyToComplex(const Vector& x, Vector& xc) {
	for(int i=0; i<N; ++i){
		xc[2*i] = x[i];
		xc[2*i+1] = 0.0;
	}
	for(int i=N; i<nDFT; ++i){
		xc[2*i] = 0.0;
		xc[2*i+1] = 0.0;
	}
}

void ConvolveGSL::complexMult(const Vector& xc, const Vector& yc, Vector& zc) {
	for(int i=0; i<nDFT; ++i){
		zc[2*i] = xc[2*i]*yc[2*i] - xc[2*i+1]*yc[2*i+1];
		zc[2*i+1] = xc[2*i]*yc[2*i+1] + xc[2*i+1]*yc[2*i];
	}
}

void ConvolveGSL::copyToReal(Vector& x, const Vector& xc) {
	for(int i=0; i<N; ++i)
		x[i] = xc[2*i];
}

void ConvolveGSL::run(const Vector& x, const Vector& y, Vector& z) {
	copyToComplex(x, xC);
	copyToComplex(y, yC);
	// xC = DFT(xC), yC = DFT(yC)
	gsl_fft_complex_forward(&(xC[0]), 1, nDFT, wavetable, workspace);
	gsl_fft_complex_forward(&(yC[0]), 1, nDFT, wavetable, workspace);
	// zC = xC * yC
	complexMult(xC, yC, zC);
	// inverse DFT zC
	gsl_fft_complex_backward(&(zC[0]), 1, nDFT, wavetable, workspace);
	// read out the real part, re-scale by 1/nDFT
	copyToReal(z, zC);
	scaleVector(z, 1.0 / nDFT);
}

void ConvolveGSL::operator ()(const Vector& x, const Vector& y, Vector& z) {
	this->run(x, y, z);
}
