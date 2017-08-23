/*
 * cinterface.cpp
 *
 *  Created on: May 26, 2017
 *      Author: bayerc
 */

#include"cinterface.h"

void computeIV_tensor(double xi, double* H, double* eta, double* rho, int par_size,
		double* T, double* K, int op_size, int N, long M, uint64_t* seed, int seed_size,
		int num_threads, double* price, double* iv, double* stat, double* H_tot, double* eta_tot,
		double* rho_tot, double* T_tot, double* K_tot, int out_size){
	// check that enough memory is available for output
	if(out_size != par_size * op_size){
		std::cerr << "Output arrays have wrong size:\n"
				<< "par_size * op_size = " << par_size * op_size
				<< " != out_size.\n";
		exit(17);
	}
	// turn C-style arrays into vectors
	Vector Hv(H, H + par_size);
	Vector etav(eta, eta + par_size);
	Vector rhov(rho, rho + par_size);
	Vector Tv(T, T + op_size);
	Vector Kv(K, K + op_size);
	std::vector<uint64_t> seedv(seed, seed + seed_size);

	// compute prices
	Result res;
	if(num_threads == 1){
		RBergomiST rBergomi(xi, Hv, etav, rhov, Tv, Kv, N, M, seedv);
		res = rBergomi.ComputeIV();
	}
	else if(num_threads > 1){
		res = ComputeIVMT(xi, Hv, etav, rhov, Tv, Kv, N, M, num_threads, seedv);
	}
	else{
		std::cerr << "Number of threads has to be a positive integer,\n"
				<< "obtained num_threads = " << num_threads << ".\n";
		exit(21);
	}

	// copy results into provided arrays
	Result2Array(res, price, iv, stat, H_tot, rho_tot, eta_tot, T_tot, K_tot);
}


void computeIVRT_tensor(double xi, double* H, double* eta, double* rho, int par_size,
		double* T, double* K, int op_size, int N, long M, uint64_t* seed, int seed_size,
		int num_threads, double* price, double* iv, double* stat, double* H_tot, double* eta_tot,
		double* rho_tot, double* T_tot, double* K_tot, int out_size){
	// check that enough memory is available for output
	if(out_size != par_size * op_size){
		std::cerr << "Output arrays have wrong size:\n"
				<< "par_size * op_size = " << par_size * op_size
				<< " != out_size.\n";
		exit(17);
	}
	// turn C-style arrays into vectors
	Vector Hv(H, H + par_size);
	Vector etav(eta, eta + par_size);
	Vector rhov(rho, rho + par_size);
	Vector Tv(T, T + op_size);
	Vector Kv(K, K + op_size);
	std::vector<uint64_t> seedv(seed, seed + seed_size);

	// compute prices
	Result res;
	if(num_threads == 1){
		RBergomiST rBergomi(xi, Hv, etav, rhov, Tv, Kv, N, M, seedv);
		res = rBergomi.ComputeIVRT();
	}
	else if(num_threads > 1){
		res = ComputeIVRTMT(xi, Hv, etav, rhov, Tv, Kv, N, M, num_threads, seedv);
	}
	else{
		std::cerr << "Number of threads has to be a positive integer,\n"
				<< "obtained num_threads = " << num_threads << ".\n";
		exit(21);
	}

	// copy results into provided arrays
	Result2Array(res, price, iv, stat, H_tot, eta_tot, rho_tot, T_tot, K_tot);
}


void Result2Array(const Result& res, double* price, double* iv, double* stat, double* H_tot, double* eta_tot,
		double* rho_tot, double* T_tot, double* K_tot){
	std::copy(res.price.begin(), res.price.end(), price);
	std::copy(res.iv.begin(), res.iv.end(), iv);
	std::copy(res.stat.begin(), res.stat.end(), stat);
	for(size_t i = 0; i < res.price.size(); ++i){
		H_tot[i] = res.par.H(i);
		rho_tot[i] = res.par.rho(i);
		eta_tot[i] = res.par.eta(i);
		T_tot[i] = res.par.T(i);
		K_tot[i] = res.par.K(i);
	}
}

void computeIVRT_tensor_quadrature(double xi, double *H, double *eta, double *rho, int par_size, double *T,
		double *K, int op_size, int N, long  M, double *Z, int num_threads, double *payoff)
{

	// turn C-style arrays into vectors
	Vector Hv(H, H + par_size);
	Vector etav(eta, eta + par_size);
	Vector rhov(rho, rho + par_size);
	Vector Tv(T, T + op_size);
	Vector Kv(K, K + op_size);

	// construct the arrays W1Arr and W1perpArr
	std::vector<Vector> W1Arr(M, Vector(N));
	std::vector<Vector> W1perpArr(M, Vector(N));
	for(size_t i=0; i<M; ++i){
		W1Arr[i] = Vector(Z + 2*i*N, Z + (2*i+1)*N);
		W1perpArr[i] = Vector(Z + (2*i+1)*N, Z + (2*i+2)*N);
	}

	// compute payoffs
	std::vector<Vector> payoffArr = ComputePayoffRTsamples(xi, Hv, etav, rhov, Tv, Kv,
			num_threads, W1Arr, W1perpArr);

	size_t out_size = payoffArr[0].size();
	// copy payoffArr into payoff
	for(size_t m =0; m<M; ++m){
		for(size_t i=0; i<out_size; ++i)
			payoff[m*out_size + i] = payoffArr[m][i];
	}
}


