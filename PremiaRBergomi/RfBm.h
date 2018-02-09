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
	double A(int i, int j) const; // covariance matrix
	double G(double x) const; // function G as in Benjamin's write-up
	// individual parts of covariance matrix
	double cBm(int i, int j) const; // covariance of scaled Brownian increments
	double cfBm(int i, int j) const; // covariance between fBm at time (i+1)/n and (j+1)/n
	double cMixed(int i, int j) const; // covariance between \sqrt{n} (W_{(i+1)/n} - W_{i/n}) and \tilde{W}_{(j+1)/n}
public:
	RfBm();
	RfBm(int nI, double HI, RNorm * rnormI);
	void generate(Vector& W1, Vector& Wtilde);
	void operator()(Vector& W1, Vector& Wtilde);
	std::vector<Vector> GetL() const {return L;}
	std::vector<Vector> GetA() const;
};

#endif /* RFBM_H_ */
