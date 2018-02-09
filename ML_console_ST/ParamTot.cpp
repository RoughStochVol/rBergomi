/*
 * ParamTot.cpp
 *
 *  Created on: Nov 17, 2016
 *      Author: bayerc
 */

#include "ParamTot.h"

ParamTot::ParamTot() {};

ParamTot::ParamTot(Vector h, Vector e, Vector r, Vector t, Vector k, Vector x) {
	// Check that all input vector have the same size.
	if((h.size() != e.size()) || (h.size() != r.size()) || (h.size() != t.size()) || (h.size() != k.size())){
		std::cerr << "Input vectors do not have the same sizes.\n";
		exit(47);
	}
	HV = h;
	etaV = e;
	rhoV = r;
	TV = t;
	KV = k;
	xiV = x;
}
