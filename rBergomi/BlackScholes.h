/*
 * BlackScholes.h
 *
 *  Created on: Apr 13, 2017
 *      Author: bayerc
 */

#pragma once

#include<cmath>
#include <boost/math/tools/roots.hpp>

double Phi(double value);
double phi(double value);
double BS_call_price(double S0, double K, double tau, double sigma, double r = 0.0);
double BS_put_price(double S0, double K, double tau, double sigma, double r = 0.0);
double BS_Vega(double S0, double K, double tau, double sigma, double r = 0.0);
double IV_call(double price, double S0, double K, double tau, double r = 0.0);
double IV_put(double price, double S0, double K, double tau, double r = 0.0);
