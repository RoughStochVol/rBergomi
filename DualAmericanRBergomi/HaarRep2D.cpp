/*
 * HaarRep2D.cpp
 *
 *  Created on: 8 May 2018
 *      Author: bayerc
 */

#include "HaarRep2D.h"
/* \alpha = (l, j, i, k)
 * We have: \iota = i + (K - k - 1)*d + K * d * int(2^(l-1) + j)
 * i = \iota % d; \iota_1 := (\iota - i)/d
 * k = K - (\iota_1 % K) - 1; \iota_2 := (\iota_1 - (K - k - 1)) / K = int(2^(l-1) + j)
 * if \iota_2 == 0 then l = j = 0;
 * else l = int(log_2(\iota_2)) + 1, j = \iota_2 - 2^(l-1)
 * */

int HaarRep2D::alpha2iota(MultiIndex alpha) const {
	return alpha[2] + (K - alpha[3] - 1) * d + K * d * static_cast<int>(pow(2, alpha[0] - 1) + alpha[1]);
}

MultiIndex HaarRep2D::iota2alpha(int iota) const {
	MultiIndex alpha;
	alpha[2] = iota % d;
	int iota1 = (iota - alpha[2]) / d;
	alpha[3] = K - (iota1 % K) - 1;
	int iota2 = (iota1 - (K - alpha[3] - 1)) / K;
	if(iota2 == 0){
		alpha[0] = 0;
		alpha[1] = 0;
	}
	else{
		alpha[0] = static_cast<int>(log2(iota2)) + 1;
		alpha[1] = iota2 - pow(2, alpha[0]);
	}
	return alpha;
}

void HaarRep2D::getWtildeKer() {
}

void HaarRep2D::getdBKer() {
}

HaarRep2D::HaarRep2D() : HierarchicalRep(){
}

HaarRep2D::HaarRep2D(double h, double t, int k, int n,
		std::vector<uint64_t> seed) : HierarchicalRep(h, t, k, n, 2, n/2, seed){
	// First check that k and n have the right properties.
	this->n = ilog2(I / (2 * K) - 1);
	if(K * 2 * (pow(2, n) + 1) != I){
		std::cerr << "N = " << N << " and K = " << K << " do not provide a full Haar basis." << std::endl;
		exit(5);
	}
	// construct auxiliary vectors
	getWtildeKer();
	getdBKer();
}

Vector HaarRep2D::Wtilde() const {
}

Vector HaarRep2D::dB(int j) const {
}
