/*
 * HierarchicalRep.cpp
 *
 *  Created on: 8 May 2018
 *      Author: bayerc
 */

#include "HierarchicalRep.h"

HierarchicalRep::HierarchicalRep() {
	H = 0.0;
	T = 0.0;
	K = 0;
	dt = 0.0;
	N = 0;
	ds = 0.0;
	d = 0;
	I = 0;
	X = Vector(0);
	setGen(std::vector<uint64_t>(0));
	dist = normDist(0.0, 1.0);
}

HierarchicalRep::HierarchicalRep(double h, double t, int k, int n, int D, long i,
		std::vector<uint64_t> seed = std::vector<uint64_t>(0)) {
	H = h;
	T = t;
	K = k;
	dt = T/K;
	N = n;
	ds = T/N;
	d = D;
	I = i;
	X = Vector(I, 0.0);
	setGen(seed);
	dist = normDist(0.0, 1.0);
}

void HierarchicalRep::setGen(std::vector<uint64_t> seed) {
	std::seed_seq seeder(seed.begin(), seed.end());
	seeder.generate(seed.begin(), seed.end());
	gen.seed(seed[0]);
}

void HierarchicalRep::generate() {
	for (auto & x : X)
		x = dist(gen);
}

Vector HierarchicalRep::getX() const {
	return X;
}
