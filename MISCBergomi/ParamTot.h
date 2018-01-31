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

class ParamTot {
private:
	double eps;
	Vector HV;
	Vector etaV;
	Vector rhoV;
	Vector TV;
	Vector KV;
	double xiV;
	long Hsize; // size of H, eta and rho
	long Tsize; // size of K and T
	long TotSize; // Hsize*Tsize
	inline long toHindex(long I) const { return I / Tsize;}
	inline long toTindex(long I) const { return I % Tsize;}
public:
	ParamTot();
	ParamTot(Vector h, Vector e, Vector r, Vector t, Vector k, double x);
	inline double H(long I) const {return HV[toHindex(I)];}
	inline double eta(long I) const {return etaV[toHindex(I)];}
	inline double rho(long I) const {return rhoV[toHindex(I)];}
	inline double T(long I) const {return TV[toTindex(I)];}
	inline double K(long I) const {return KV[toTindex(I)];}
	inline long size() const {return TotSize;}
	inline double xi() const {return xiV;}
	// Is the current H value different from the one before (i.e., for I-1)?
	inline bool HTrigger(long I) const {return (I == 0) || (fabs(H(I) - H(I-1)) > eps);}
	inline bool etaTrigger(long I) const {return (I == 0) || (fabs(eta(I) - eta(I-1)) > eps);}
	inline bool rhoTrigger(long I) const {return (I == 0) || (fabs(rho(I) - rho(I-1)) > eps);}
	inline bool TTrigger(long I) const {return (I == 0) || (fabs(T(I) - T(I-1)) > eps);}
	// note that K is always supposed to be triggered!
};

class ParamTotUnordered{
private:
	Vector HV;
	Vector etaV;
	Vector rhoV;
	Vector TV;
	Vector KV;
	double xiV;
public:
	ParamTotUnordered();
	ParamTotUnordered(Vector h, Vector e, Vector r, Vector t, Vector k, double x);
	inline double H(long I) const {return HV[I];}
	inline double eta(long I) const {return etaV[I];}
	inline double rho(long I) const {return rhoV[I];}
	inline double T(long I) const {return TV[I];}
	inline double K(long I) const {return KV[I];}
	inline long size() const {return HV.size();}
	inline double xi() const {return xiV;}
};

#endif /* PARAMTOT_H_ */
