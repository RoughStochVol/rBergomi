/*
 * BlackScholes.cpp
 *
 * Auxiliary functions for computing Black-Scholes prices and implied vols.
 * 
 * In the next step, try to export these function from the library
 *
 *  Created on: Apr 13, 2017
 *      Author: bayerc
 */

#include<Rcpp.h>
#include "BlackScholes.h"

using namespace Rcpp;
using boost::math::tools::eps_tolerance;
using boost::math::tools::bracket_and_solve_root;

// normal CDF
double Phi(double value)
{
   return 0.5 * erfc(-value * M_SQRT1_2);
}

// BS call price.
// [[Rcpp::export]]
double BS_call_price(double S0, double K, double tau, double sigma, double r){
	double d1 = (log(S0/K) + (r+0.5*sigma*sigma)*tau)/(sigma*sqrt(tau));
	double d2 = d1 - sigma*sqrt(tau);
	return Phi(d1)*S0 - Phi(d2)*K*exp(-r*tau);
}

// BS put price.
// [[Rcpp::export]]
double BS_put_price(double S0, double K, double tau, double sigma, double r){
	double d1 = (log(S0/K) + (r+0.5*sigma*sigma)*tau)/(sigma*sqrt(tau));
	double d2 = d1 - sigma*sqrt(tau);
	return -Phi(-d1)*S0 + Phi(-d2)*K*exp(-r*tau);
}

// implied vol
// [[Rcpp::export]]
double IV_call(double price, double S0, double K, double tau, double r){
	// if ITM, apply put-call-parity and compute IV on the put option
	if(S0 > K){
		double put_price = price - S0 + exp(-r*tau) * K;
		return IV_put(put_price, S0, K, tau, r);
	}
	// Apply the bracketing root solver from boost.

	// first construct the objective function; we use a lambda
	auto objFunc = [&](double sigma) {return price - BS_call_price(S0, K, tau, sigma, r);};

	// set initial guess
	double guess = 0.2;

	// set various numerical parameters for the root finder
	double factor = 2.0;// size of steps for searching??
	const boost::uintmax_t maxit = 20;            // Limit to maximum iterations.
	boost::uintmax_t it = maxit;                  // Initally our chosen max iterations, but updated with actual.
	bool is_rising = false;                        // Is the objective function increasing? (Otherwise: decreasing)
	int get_digits = 11;                 // tolerance = 2^(1-get_digits)
	eps_tolerance<double> tol(get_digits);             // Set the tolerance.

	// run the solver
	std::pair<double, double> res = bracket_and_solve_root(objFunc, guess, factor, is_rising, tol, it);

	// return the middle value in the interval
	return res.first + (res.second - res.first)/2;
}

// [[Rcpp::export]]
double IV_put(double price, double S0, double K, double tau, double r){
	// if ITM, apply put-call-parity and compute IV on the call option
	if(S0 < K){
		double call_price = price + S0 - exp(-r*tau) * K;
		return IV_call(call_price, S0, K, tau, r);
	}
	// Apply the bracketing root solver from boost.

	// first construct the objective function; we use a lambda
	auto objFunc = [&](double sigma) {return price - BS_put_price(S0, K, tau, sigma, r);};

	// set initial guess
	double guess = 0.2;

	// set various numerical parameters for the root finder
	double factor = 2.0;// size of steps for searching??
	const boost::uintmax_t maxit = 20;            // Limit to maximum iterations.
	boost::uintmax_t it = maxit;                  // Initally our chosen max iterations, but updated with actual.
	bool is_rising = false;                       // Is the objective function increasing? (Otherwise: decreasing)
	int get_digits = 11;                 //tolerance = 2^(1-get_digits)
	eps_tolerance<double> tol(get_digits);             // Set the tolerance.

	// run the solver
	std::pair<double, double> res = bracket_and_solve_root(objFunc, guess, factor, is_rising, tol, it);

	// return the middle value in the interval
	return res.first + (res.second - res.first)/2;
}

// [[Rcpp::export]]
List rcpp_hello_world() {
  
  CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
  NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
  List z            = List::create( x, y ) ;
  
  return z ;
}

// [[Rcpp::export]]
double stupid_plus(double a, double b){
  return a+b;
}

// [[Rcpp::export]]
double stupid_one(){
  return 1.0;
}

// [[Rcpp::export]]
List stupid_two(){
  
  List z = List::create(2.0);
  return z;
}