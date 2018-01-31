/*
 * RfBm.cpp
 *
 *  Created on: 30 Jan 2018
 *      Author: bayerc
 */

#include "RfBm.h"

RfBm::RfBm() {
	n = 0;
	H = 0.0;
	gamma = 0.0;
	rnorm = 0;
	L = std::vector<Vector>(0);
}

double RfBm::A(int i, int j) {
	double ret = 0.0;
	if ((i < n) && (j < n)) {
		if (i == j)
			ret = 1.0;
	} else if ((i < n) && (j >= n)) {
		if (i <= (j - n))
			ret = 2.0 * sqrt(2.0 * H) * pow(double(j - n) + 1.0, 0.5 + H)
					* pow(double(n), -H) / (1.0 + 2.0 * H);
	} else if ((i >= n) && (j < n)) {
		if (j <= (i - n))
			ret = 2.0 * sqrt(2.0 * H) * pow(double(i - n) + 1.0, 0.5 + H)
					* pow(double(n), -H) / (1.0 + 2.0 * H);
	} else {
		if (i < j)
			ret = pow(double(i + 1) / double(n), 2.0 * H)
					* G(double(j + 1) / double(i + 1));
		else if (i > j)
			ret = pow(double(j + 1) / double(n), 2.0 * H)
					* G(double(i + 1) / double(j + 1));
		else if (i == j)
			ret = pow(double(j + 1) / double(n), 2.0 * H);
	}
	return ret;
}

double RfBm::G(double x) {
	return 2.0 * H
			* (pow(x, -gamma) / (1.0 - gamma)
					+ gamma * pow(x, -(1.0 + gamma))
							* gsl_sf_hyperg_2F1(1.0, 1.0 + gamma, 3.0 - gamma,
									1.0 / x) / ((1.0 - gamma) * (2.0 - gamma)));
}

RfBm::RfBm(int nI, double HI, RNorm* rnormI) {
	n = nI;
	H = HI;
	gamma = 0.5 - H;
	rnorm = rnormI;
	L = std::vector<Vector>(2 * n, Vector(2 * n, 0.0)); // fill with 0s first.
	// Now fill the matrix L.
	for (int k = 0; k < 2 * n; ++k) {
		double temp = 0.0;
		for (int j = 0; j < k; ++j)
			temp += L[k][j] * L[k][j];
		L[k][k] = sqrt(A(k, k) - temp);
		for (int i = k + 1; i < n; ++i) {
			double temp = 0.0;
			for (int j = 0; j < k; ++j)
				temp += L[i][j] * L[k][j];
			L[i][k] = (A(i, k) - temp) / L[k][k];
		}
	}
}

void RfBm::generate(Vector& W1, Vector& Wtilde) {
	Vector x(2 * n);
	(*rnorm)(x);
	for (int i = 0; i < n; ++i) {
		W1[i] = 0.0;
		for (int j = 0; j <= i; ++j)
			W1[i] += L[i][j] * x[j];
	}
	for (int i = 0; i < n; ++i) {
		Wtilde[i] = 0.0;
		for (int j = 0; j <= (i + n); ++j)
			Wtilde[i] += L[i + n][j] * x[j];
	}
}

void RfBm::operator ()(Vector& W1, Vector& Wtilde) {
	this->generate(W1, Wtilde);
}
