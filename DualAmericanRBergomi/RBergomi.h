/*
 * RBergomi.h
 *
 * Simulation of rBergomi trajectories based on hierarchical representation.
 *
 * NEED TO USE STRATONOVICH CORRECTION BETWEEN S AND V!!! (Using the regularity structures paper.)
 *
 *  Created on: 12 Jun 2018
 *      Author: bayerc
 */

#ifndef RBERGOMI_H_
#define RBERGOMI_H_

#include <functional>
#include "HierarchicalRep.h"

class RBergomi {
private:
	HierarchicalRep* rep;
	int N; // # of timesteps
	double ds; // timestep
	double S0; // spot price
	double r; // interest rate respectively negative dividend rate
	double H; // Hurst index
	double eta; // vol of vol
	double rho; // correlation
	double rhobar; // = sqrt(1 - rho^2)
	std::function<double(double)> xi0; // forward variance curve
public:
	RBergomi();
	RBergomi(HierarchicalRep* hier, double s0, double R, double Eta, double Rho,
			std::function<double(double)> f);
	virtual ~RBergomi(){};
	// Return one simulated path of the rBergomi process
	// Note: this changes the inner state of rep, as it generates a new sample X.
	Vector generate();
	// For testing purposes: return both S and V
	std::pair<Vector, Vector> generate_SV();
	// Return one simulated path of S and the corresponding vector X of standard normals
	// used in the hierarchical representation underlying the Brownian motions.
	// Note: this changes the inner state of rep, as it generates a new sample X.
	void operator()(Vector& S, Vector& X);
};

#endif /* RBERGOMI_H_ */
