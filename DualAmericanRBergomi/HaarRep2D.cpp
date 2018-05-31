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
		alpha[1] = iota2 - pow(2, alpha[0] - 1);
	}
	return alpha;
}

void HaarRep2D::getWtildeKer() {
	WtildeKer = std::vector<Vector>(N, Vector(I));
	for(int n = 0; n < N; ++n){
		for(int iota = 0; iota < I; ++iota){
			MultiIndex alpha = iota2alpha(iota);
			if(alpha[2] == 0){
				WtildeKer[n][iota] = PsiTilde(n * ds, alpha[0], alpha[1]);
			}
			else{
				WtildeKer[n][iota] = 0.0;
			}
		}
	}
}

void HaarRep2D::getdBKer() {
	dB1Ker = std::vector<Vector>(N, Vector(I));
	dB2Ker = std::vector<Vector>(N, Vector(I));
	for(int n = 0; n < N; ++n){
		for(int iota = 0; iota < I; ++iota){
//			MultiIndex alpha = iota2alpha(iota);
//			if(alpha[2] == 0){
//				dB2Ker[n][iota] = 0.0;
//				dB1Ker[n][iota] = Psi((n+1) * ds, alpha[0], alpha[1]) - Psi(n * ds, alpha[0], alpha[1]);
//			}
//			else{
//				dB1Ker[n][iota] = 0.0;
//				dB2Ker[n][iota] = Psi((n+1) * ds, alpha[0], alpha[1]) - Psi(n * ds, alpha[0], alpha[1]);
//			}
			dB1Ker[n][iota] = Phi((n+1) * ds, iota, 0) - Phi(n * ds, iota, 0);
			dB2Ker[n][iota] = Phi((n+1) * ds, iota, 1) - Phi(n * ds, iota, 1);
		}
	}
}

double HaarRep2D::Psi(double t, int l, int j) const{
	double ret;
	if(l == 0)
		ret = t;
	else{
		double tMin = pow(2.0, -(l-1)) * j;
		double tMid = pow(2.0, -(l-1)) * (j + 0.5);
		double tMax = pow(2.0, -(l-1)) * (j + 1.0);
		ret = pow(2.0, 0.5*(l-1)) * (2.0 * fmin(tMid, t) - fmin(tMin, t) - fmin(tMax, t));
	}
	return ret;
}

double HaarRep2D::Phi(double t, int iota, int i) const{
	const double eps = 0.0000001;
	MultiIndex alpha = iota2alpha(iota);
	double ret;
	// return 0 if iota acts in the other dimension or t is not contained in the active exercise interval.
	if((alpha[2] != i) || (alpha[3] * dt > t + eps)) // bug!!!
		ret = 0.0;
	else if((alpha[3] + 1) * dt < t  - eps){
		// Take value of Psi at the end of the exercise interval.
		ret = Psi(1.0, alpha[0], alpha[1]) * sdt;
	}
	else{
		// compute the re-scaled and shifted version of Psi
		double tScaled = (t - alpha[3] * dt) / dt;
		ret = Psi(tScaled, alpha[0], alpha[1]) * sdt;
	}
	return ret;
}

double HaarRep2D::PsiTilde(double t, int l, int j) const{
	double ret;
	if(l == 0)
		ret = sqrt(2*H/T) / (H + 0.5) * pow(t, H + 0.5); // to be changed
	else{
		double tMin = pow(2.0, -(l-1)) * j * T;
		double tMid = pow(2.0, -(l-1)) * (j + 0.5) * T;
		double tMax = pow(2.0, -(l-1)) * (j + 1.0) * T;
		double factor = sqrt(2*H/T) / (H + 0.5) * pow(2.0, 0.5*(l-1));
		if(t < tMin)
			ret = 0.0;
		else if(t < tMid)
			ret = factor * pow(t - tMin, H + 0.5);
		else if(t < tMax)
			ret = factor * (pow(t - tMin, H + 0.5) - 2.0 * pow(t - tMid, H + 0.5));
		else
			ret = factor * (pow(t - tMin, H + 0.5) - 2.0 * pow(t - tMid, H + 0.5) +
					pow(t - tMax, H + 0.5));
	}
	return ret;
}

HaarRep2D::HaarRep2D() : HierarchicalRep(){
}

HaarRep2D::HaarRep2D(double h, double t, int k, int n,
		std::vector<uint64_t> seed) : HierarchicalRep(h, t, k, n, 2, n*2, seed){
	// First check that k and n have the right properties.
	this->n = ilog2(I / (2 * K));
	sdt = sqrt(dt);
	if(K * 2 * pow(2, this->n) != I){
		std::cerr << "N = " << N << " and K = " << K << " give n = " << this->n << " and do not provide a full Haar basis." << std::endl;
		exit(5);
	}
	// construct auxiliary vectors
	getWtildeKer();
	getdBKer();
}

Vector HaarRep2D::Wtilde() const {
	Vector W(N);
	for(int n = 0; n < N; ++n)
		W[n] = std::inner_product(WtildeKer[n].begin(), WtildeKer[n].end(), X.begin(), 0.0);
	return W;
}

Vector HaarRep2D::dB(int j) const {
	Vector dB(N);
	const std::vector<Vector>& dBKer = (j == 0 ? dB1Ker : dB2Ker);
	for(int n = 0; n < N; ++n)
		dB[n] = std::inner_product(dBKer[n].begin(), dBKer[n].end(), X.begin(), 0.0);
	return dB;
}
