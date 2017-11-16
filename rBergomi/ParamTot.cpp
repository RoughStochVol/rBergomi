/*
 * ParamTot.cpp
 *
 *  Created on: Nov 17, 2016
 *      Author: bayerc
 */

#include "ParamTot.h"

ParamTot::ParamTot() {
	Tsize = 0;
	Hsize = 0;
	TotSize = 0;
	eps = 0.0;
	xiV = 0.0;
}

ParamTot::ParamTot(Vector h, Vector e, Vector r, Vector t, Vector k, double x) {
	HV = h;
	etaV = e;
	rhoV = r;
	TV = t;
	KV = k;
	Hsize = HV.size();
	Tsize = TV.size();
	TotSize = Hsize * Tsize;
	eps = 0.0000000001;
	xiV = x;
}

ParamTotUnordered::ParamTotUnordered() {
	xi = 0.0;
}

ParamTotUnordered::ParamTotUnordered(Vector h, Vector e, Vector r, Vector t, Vector k, double x) {
	// Check that all input vector have the same size.
	if((h.size() != e.size()) || (h.size() != r.size()) || (h.size() != t.size()) || (h.size() != k.size())){
		std::cerr << "Input vectors do not have the same sizes.\n";
		exit(47);
	}
	H = h;
	eta = e;
	rho = r;
	T = t;
	K = k;
	xi = x;
}
