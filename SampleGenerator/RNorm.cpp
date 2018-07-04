/*
 * RNorm.cpp
 *
 *  Created on: 26 Jan 2018
 *      Author: bayerc
 */

#include "RNorm.h"

RNorm::RNorm(uint64_t seed) {
	gen.seed(seed);
	dist = normDist(0.0, 1.0);
}

void RNorm::generate(Vector& x) {
	for(auto& v : x)
		v = dist(gen);
}

void RNorm::operator ()(Vector& x) {
	this->generate(x);
}
