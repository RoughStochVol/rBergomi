/*
 * ParamTot.h
 *
 * All possible combinations of (H, eta, rho) and (T, K), saved in order of H and then eta.
 *
 * It may be better to re-order: H ordered first, then T, then eta, then rho, then K.
 *
 * ParamTotUnordered makes no assumption on oredering of parameter: instead,
 * no combinations of parameter vectors are considered.
 *
 *  Created on: Nov 17, 2016
 *      Author: bayerc
 */

#ifndef PARAMTOT_H_
#define PARAMTOT_H_

#include<vector>
#include<cmath>
#include<iostream>
#include<cstdlib>

typedef std::vector<double> Vector;

class ParamTot{
private:
	Vector HV;
	Vector etaV;
	Vector rhoV;
	Vector TV;
	Vector KV;
	Vector xiV;
public:
	ParamTot();
	ParamTot(Vector h, Vector e, Vector r, Vector t, Vector k, Vector x);
	inline double H(long I) const {return HV[I];}
	inline double eta(long I) const {return etaV[I];}
	inline double rho(long I) const {return rhoV[I];}
	inline double T(long I) const {return TV[I];}
	inline double K(long I) const {return KV[I];}
	inline long size() const {return HV.size();}
	inline double xi(long I) const {return xiV[I];}
};

#endif /* PARAMTOT_H_ */
