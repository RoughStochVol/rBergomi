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
