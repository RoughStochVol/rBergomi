/*
 * RfBm.h
 *
 * Generate samples from W1 and Wtilde jointly as provided in compute_Wtilde in rBergomi.cpp.
 *
 * n denotes the number N in rBergomi.cpp, i.e., the timestep is 1/n.
 *
 *  Created on: 30 Jan 2018
 *      Author: bayerc
 */

#ifndef RFBM_H_
#define RFBM_H_

#include <gsl/gsl_sf_hyperg.h>
#include "RNorm.h"

class RfBm {
private:
	int n;
	double H;
	double gamma;
	RNorm * rnorm;
	std::vector<Vector> L; // Cholesky factorization of the covariance matrix.
	double A(int i, int j); // covariance matrix
	double G(double x); // function G as in Benjamin's write-up
public:
	RfBm();
	RfBm(int nI, double HI, RNorm * rnormI);
	void generate(Vector& W1, Vector& Wtilde);
	void operator()(Vector& W1, Vector& Wtilde);
};

#endif /* RFBM_H_ */
