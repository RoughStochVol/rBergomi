/*
 * aux.cpp
 *
 *  Created on: 31 May 2018
 *      Author: bayerc
 */

#include "aux.h"

double G(double x, double gamma) {
	double ret;
	if (x > 1.0)
		ret = (1.0 - 2.0 * gamma) * pow(x, -gamma)
				* gsl_sf_hyperg_2F1(gamma, 1.0, 2.0 - gamma, 1.0 / x)
				/ (1.0 - gamma);
	else
		ret = 1.0;
	return ret;
}

double covWtilde(double t, double s, double H){
	double tMin = fmin(t, s);
	double tMax = fmax(t, s);
	return pow(tMin, 2.0*H) * G(tMax / tMin, 0.5 - H);
}

double covCross(double t, double s, double H, int i){
	double ret;
	if(i != 0)
		ret = 0.0;
	else
		ret = sqrt(2.0*H) * (pow(t, H + 0.5) - pow(t - fmin(t,s), H+0.5)) / (H + 0.5);
	return ret;
}
