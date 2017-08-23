/*
 * cinterface.h
 *
 * C-style interface to rBergomi.
 *
 *  Created on: May 26, 2017
 *      Author: bayerc
 */

#pragma once

/*
 * computeIV_tensor
 *
 * Compute implied vols using MC simulation based on the hybrid scheme in the rBergomi model.
 *
 * Input:
 * xi ... the constant forward variance curve
 * H, eta, rho ... arrays of size par_size contianing the model parameters. For good performance,
 * 				   rho should change fastest and H slowest if these parameters are the output
 * 				   of a tensorization.
 * T, K ... arrays of size op_size, contianing the maturities and strikes of the options to be computed.
 * 			Note that option prices are computed for all possible combinations of option parameters
 * 			(T,K) and model parameters (H, eta, rho). In other words, for every model parameter triple,
 * 			the full vol surface is computed.
 * N, M ... number of steps and trajectories for the Euler-MC scheme.
 * seed ... array of size seed_size containing the seed for the RNG.
 * num_threads ... number of threads. If num_threads == 1, the single-threaded code is called.
 *
 * Output:
 * price ... array of size out_size containing the call prices. Note that out_size should be equal to
 * 			 op_size * par_size.
 * iv ... implied volatilities of size out_size.
 * H_tot,..., K_tot ... parameter array of size out_size. price[i] corresponds to the call option price
 * 						with model parameters H_tot[i], rho_tot[i], eta_tot[i], and maturity T_tot[i]
 * 						and strike price K_tot[i].
 *
 * computeIVRT_tensor does the same but uses the Romano-Touzi trick.
 */

extern "C"{
void computeIV_tensor(double xi, double* H, double* eta, double* rho, int par_size,
		double* T, double* K, int op_size, int N, long M, uint64_t* seed, int seed_size,
		int num_threads, double* price, double* iv, double* stat, double* H_tot, double* eta_tot,
		double* rho_tot, double* T_tot, double* K_tot, int out_size);

void computeIVRT_tensor(double xi, double* H, double* eta, double* rho, int par_size,
		double* T, double* K, int op_size, int N, long M, uint64_t* seed, int seed_size,
		int num_threads, double* price, double* iv, double* stat, double* H_tot, double* eta_tot,
		double* rho_tot, double* T_tot, double* K_tot, int out_size);
}
