/*
 * RBergomi.cpp
 *
 *  Created on: 12 Jun 2018
 *      Author: bayerc
 */

#include "RBergomi.h"

RBergomi::RBergomi() {
	rep = 0;
	N = 0;
	ds = 0.0;
	r = 0.0;
	H = 0.0;
	eta = 0.0;
	rho = 0.0;
	rhobar = 0.0;
	S0 = 0.0;
	xi0 = [](double t) {return 0.0;};
}

RBergomi::RBergomi(HierarchicalRep* hier, double s0, double R, double Eta,
		double Rho, std::function<double(double)> f) {
	rep = hier;
	N = rep->getN();
	ds = rep->getDs();
	S0 = s0;
	r = R;
	H = rep->getH();
	eta = Eta;
	rho = Rho;
	rhobar = sqrt(1.0 - rho * rho);
	xi0 = f;
}

Vector RBergomi::generate() {

	return std::get<0>(generate_SV());
}

void RBergomi::operator()(Vector& S, Vector& X) {
	S = generate();
	X = rep->getX();
}

std::pair<Vector, Vector> RBergomi::generate_SV() {
	// generate random numbers
	rep->generate();
	Vector Wtilde = rep->Wtilde();
	Vector dB1 = rep->dB(0);
	Vector dB2 = rep->dB(1);
	Vector S(N + 1);
	Vector X(N + 1); // log-price

	// define the renormalisation constant; WHY DOES IT HAVE NO EFFECT???
	const double Cepsilon = sqrt(2 * H) / ((H + 0.5) * (H + 1.5)) *
			pow(2.0, log2(rep->getN() / rep->getK()) * (0.5 - H));
//	const double Cepsilon = sqrt(2 * H) / ((H + 0.5) * (H + 1.5)) *
//				pow(2.0, rep->getN() / rep->getK() * (0.5 - H));
//	const double Cepsilon = 0.5;
	// initialize state variables
	X[0] = log(S0);
	Vector v(N); // do not need last value, but logic is easier this way.
	v[0] = xi0(0.0);
	// Maybe I should not use the left-point rule here!!!
	for (int n = 0; n < N - 1; ++n) {
		v[n + 1] = xi0((n + 1) * ds)
				* exp(
						eta * Wtilde[n + 1]
								- 0.5 * eta * eta * pow((n + 1) * ds, 2.0 * H));
		double vMid = 0.5 * (v[n] + v[n+1]);
		X[n + 1] = X[n] + (r - 0.5 * vMid) * ds
				+ sqrt(vMid) * (rho * dB1[n] + rhobar * dB2[n])
				- 0.5 * Cepsilon * eta * rho * sqrt(vMid) * ds;
	}
	// BUG! SHOULD REALLY USE MIDPOINT RULE AT LAST TIME STEP, TOO.
	X[N] = X[N-1] + (r - 0.5 * v[N - 1]) * ds
			+ sqrt(v[N - 1]) * (rho * dB1[N - 1] + rhobar * dB2[N - 1]) - 0.5 * Cepsilon * eta * rho * sqrt(v[N-1]) * ds;
	// S = exp(X)
	std::transform(X.begin(), X.end(), S.begin(),
			[](double x) -> double {return exp(x);});

	return std::pair<Vector, Vector> { S, v };
}
