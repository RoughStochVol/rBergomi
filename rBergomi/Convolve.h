/*
 * Convolve.h
 *
 * Class for computing discrete convolutions as used for the hybrid scheme.
 * There will be an abstract base class Convolve, and (at least) one derived class
 * ConvolveFFTW, which implements convolutions based on FFTW.
 *
 * Note that the class cannot assumed to be thread-save, inheriting the problems for fftw.
 * Note further that each dimension N requires a new object for convolution.
 *
 *  Created on: 31 Jan 2018
 *      Author: bayerc
 */

#ifndef CONVOLVE_H_
#define CONVOLVE_H_

#include <vector>
#include <fftw3.h>
#include <gsl/gsl_fft_complex.h>
#include "aux.h"

typedef std::vector<double> Vector;

class Convolve {
protected:
	int N; // size of vectors to be convolved
public:
	Convolve(){N = 0;}
	Convolve(int n);
	virtual ~Convolve(){};
	virtual void run(const Vector& x, const Vector& y, Vector& z) = 0; // z = x * y
	virtual void operator()(const Vector& x, const Vector& y, Vector& z) = 0; // z = x * y
};

class ConvolveFFTW : public Convolve {
private:
	int nDFT;
	fftw_complex *xC; // complex array (of size nDFT = 2*N-1)
	fftw_complex *xHat;
	fftw_complex *yC;
	fftw_complex *yHat;
	fftw_complex *zC;
	fftw_complex *zHat;
	fftw_plan fPlanX; // encodes xHat = DFT(xC)
	fftw_plan fPlanY; // yHat = DFT(yC)
	fftw_plan fPlanZ; // zC = IDFT(zC)
	void copyToComplex(const Vector & x, fftw_complex *xc);
	void fftw_c_mult(const fftw_complex a, const fftw_complex b, fftw_complex c);
	void copyToReal(Vector & x, const fftw_complex *xc);
	void complexMult(const fftw_complex *x, const fftw_complex *y,
			fftw_complex *z);
public:
	ConvolveFFTW();
	ConvolveFFTW(int N);
	ConvolveFFTW(ConvolveFFTW& conv);
	ConvolveFFTW& operator=(const ConvolveFFTW& conv);
	~ConvolveFFTW();
	void run(const Vector& x, const Vector& y, Vector& z); // z = x * y
	void operator()(const Vector& x, const Vector& y, Vector& z); // z = x * y
};

class ConvolveGSL : public Convolve {
private:
	int nDFT;
	Vector xC; // used for the FFT; As it is complex, the size is 2*nDFT, with even indeces for real, odd ones for imaginary parts
	Vector yC;
	Vector zC;
	gsl_fft_complex_wavetable* wavetable;
	gsl_fft_complex_workspace* workspace;
	void copyToComplex(const Vector & x, Vector& xc);
	void complexMult(const Vector& xc, const Vector& yc, Vector& zc);
	void copyToReal(Vector& x, const Vector& xc);
public:
	ConvolveGSL();
	ConvolveGSL(int N);
	ConvolveGSL(ConvolveGSL& conv);
	ConvolveGSL& operator=(const ConvolveGSL& conv);
	~ConvolveGSL();
	void run(const Vector& x, const Vector& y, Vector& z); // z = x * y
	void operator()(const Vector& x, const Vector& y, Vector& z); // z = x * y
};

#endif /* CONVOLVE_H_ */
