/*
 * RNorm.h
 *
 *  Created on: 26 Jan 2018
 *      Author: bayerc
 */

#ifndef RNORM_H_
#define RNORM_H_

#include<vector>
#include<random>

typedef std::vector<double> Vector;
typedef std::mt19937_64 MTGenerator;
typedef std::normal_distribution<double> normDist;

class RNorm {
private:
	MTGenerator gen;
	normDist dist;
public:
	RNorm(uint64_t seed = 0);
	void generate(Vector& x);
	void operator()(Vector& x);
};

#endif /* RNORM_H_ */
