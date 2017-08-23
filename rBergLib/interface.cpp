/*
 * interface.cpp
 *
 *  Created on: May 13, 2017
 *      Author: bayerc
 */

#include"interface.h"


ResultSimple ComputeIV_singleThread(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K, int N, long  M, std::vector<uint64_t> seed)
{
	RBergomiST rBergomi(xi, H, eta, rho, T, K, N, M, seed);
	Result res = rBergomi.ComputeIV();
	return Result2ResultSimple(res);
}

ResultSimple ComputeIVRT_multiThread(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K, int N, long  M, int numThreads, std::vector<uint64_t> seed)
{
	Result res = ComputeIVRTMT(xi, H, eta, rho, T, K, N, M, numThreads, seed);
	return Result2ResultSimple(res);
}

ResultSimple ComputeIV_multiThread(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K, int N, long  M, int numThreads, std::vector<uint64_t> seed)
{
	Result res = ComputeIVMT(xi, H, eta, rho, T, K, N, M, numThreads, seed);
	return Result2ResultSimple(res);
}

ResultSimple ComputeIVRT_singleThread(double xi, Vector H, Vector eta, Vector rho, Vector T, Vector K, int N, long  M, std::vector<uint64_t> seed)
{
	RBergomiST rBergomi(xi, H, eta, rho, T, K, N, M, seed);
	Result res = rBergomi.ComputeIVRT();
	return Result2ResultSimple(res);
}

ResultSimple Result2ResultSimple(const Result & r)
{
	Vector HTot(r.price.size(), 0.0);
	Vector rhoTot(r.price.size(), 0.0);
	Vector etaTot(r.price.size(), 0.0);
	Vector TTot(r.price.size(), 0.0);
	Vector KTot(r.price.size(), 0.0);
	for(size_t i=0; i<HTot.size(); ++i){
		HTot[i] = r.par.H(i);
		rhoTot[i] = r.par.rho(i);
		etaTot[i] = r.par.eta(i);
		TTot[i] = r.par.T(i);
		KTot[i] = r.par.K(i);
	}
	return ResultSimple{r.price, r.iv, HTot, rhoTot, etaTot, TTot, KTot, r.par.xi(), r.stat, r.N, r.M, r.numThreads};
}

// test this
double stupidS(std::vector<double> x){
	double res = 0.0;
	for(size_t i=0; i<x.size(); ++i)
		res += x[i];
	return res;
}

ParSimple::ParSimple()
{
	n_H = 0;
	H = new double[n_H];
}

ParSimple::ParSimple(double *h, size_t n_h)
{
	n_H = n_h;
	H = new double[n_H];
	std::copy(h, h+n_H, H);
}

ParSimple::~ParSimple()
{
	delete[] H;
}

Vector ParSimple::getH() const{
	return Vector(H, H+n_H);
}

ResultSimple stupid(ParSimple x){
	ResultSimple res;
	res.H = x.getH();
	return res;
}
