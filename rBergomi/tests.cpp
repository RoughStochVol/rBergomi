/*
 * tests.cpp
 *
 * Unit tests for RBergomi and the multi-threaded version.
 *
 *  Created on: Apr 13, 2017
 *      Author: bayerc
 */

//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include"catch.hpp"
#include "RBergomi.h"
#include "rBergomiMT.h"
#include "BlackScholes.h"

TEST_CASE( "Test the single and multi-threaded pricing routines", "[pricing]" )
{
	const double xi = 0.07;
	const Vector H { 0.07, 0.1 };
	Vector eta { 2.2, 2.4 };
	Vector rho { -0.9, -0.85 };
	Vector T { 0.5, 1.0 };
	Vector K { 0.8, 1.0 };
	Vector truePrice { 0.2149403, 0.08461882, 0.2157838, 0.07758564 };
	long steps = 256;
	long M = 40 * static_cast<long> (pow(static_cast<double> (steps), 2.0));
	std::vector<uint64_t> seed { 123, 452, 567, 248, 9436, 675, 194, 6702 };
	RBergomiST rBergomi(xi, H, eta, rho, T, K, steps, M, seed);
	const int numThreads = 8;

	const double epsilon = 0.005;

	SECTION("Single-threaded pricing:") {
		// Check that the price is within epsilon from the confidence interval around the true price
		//Result res = rBergomi.ComputePrice();
		Result res = rBergomi.ComputeIV();
		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(res.price[i] - truePrice[i]) < epsilon + 2*res.stat[i]);

		Vector trueVol(truePrice.size());
		trueVol[0] = IV_call(truePrice[0], 1.0, K[0], T[0]);
		trueVol[1] = IV_call(truePrice[1], 1.0, K[1], T[1]);
		trueVol[2] = IV_call(truePrice[2], 1.0, K[0], T[0]);
		trueVol[3] = IV_call(truePrice[3], 1.0, K[1], T[1]);

		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(res.iv[i] - trueVol[i]) < epsilon + 2*res.stat[i]);
	}

	SECTION("Multi-threaded pricing:") {
		Result res = ComputeIVMT(xi, H, eta, rho, T, K, steps, M, numThreads,
				seed);

		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(res.price[i] - truePrice[i]) < epsilon + 2*res.stat[i]);

		Vector trueVol(truePrice.size());
		trueVol[0] = IV_call(truePrice[0], 1.0, K[0], T[0]);
		trueVol[1] = IV_call(truePrice[1], 1.0, K[1], T[1]);
		trueVol[2] = IV_call(truePrice[2], 1.0, K[0], T[0]);
		trueVol[3] = IV_call(truePrice[3], 1.0, K[1], T[1]);

		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(res.iv[i] - trueVol[i]) < epsilon + 2*res.stat[i]);
	}


	SECTION("Multi-threaded pricing with Romano-Touzi trick:") {
		Result res = ComputeIVRTMT(xi, H, eta, rho, T, K, steps, M, numThreads,
				seed);

		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(res.price[i] - truePrice[i]) < epsilon + 2*res.stat[i]);

		Vector trueVol(truePrice.size());
		trueVol[0] = IV_call(truePrice[0], 1.0, K[0], T[0]);
		trueVol[1] = IV_call(truePrice[1], 1.0, K[1], T[1]);
		trueVol[2] = IV_call(truePrice[2], 1.0, K[0], T[0]);
		trueVol[3] = IV_call(truePrice[3], 1.0, K[1], T[1]);

		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(res.iv[i] - trueVol[i]) < epsilon + 2*res.stat[i]);
	}

	SECTION("Multi-threaded pricing with Romano-Touzi trick based on externally provided samples:"){
		// generate the normals
		std::vector<Vector> W1Arr(M, Vector(steps));
		std::vector<Vector> W1perpArr(M, Vector(steps));
		RNG rng(1, seed);
		for(long i = 0; i<M; ++i){
			genGaussianMT(W1Arr[i], rng, 0);
			genGaussianMT(W1perpArr[i], rng, 0);
		}
		// now call the function to compute all the samples.
		std::vector<Vector> payoffArr = ComputePayoffRTsamples(xi, H, eta, rho, T, K,
				numThreads, W1Arr, W1perpArr);

		// now compute the prices and statistical error estimates
		size_t L = payoffArr[0].size();
		Vector price(L, 0.0);
		Vector stat(L, 0.0);
		for(long i=0; i<M; ++i){
			for(size_t j = 0; j<L; ++j){
				price[j] += payoffArr[i][j];
				stat[j] += payoffArr[i][j] * payoffArr[i][j];
			}
		}
		for(size_t j = 0; j<L; ++j){
			price[j] = price[j] / static_cast<double>(M);
			stat[j] = sqrt(stat[j] / static_cast<double>(M) - price[j]*price[j]) / sqrt(static_cast<double>(M));
		}
		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(price[i] - truePrice[i]) < epsilon + 2*stat[i]);
	}

	SECTION("Multi-threaded pricing with Romano-Touzi trick based on externally provided samples in hierarchical representation:"){
		// generate the normals
		std::vector<Vector> Z1Arr(M, Vector(steps));
		std::vector<Vector> Z1perpArr(M, Vector(steps));
		std::vector<Vector> W1Arr(M, Vector(steps));
		std::vector<Vector> W1perpArr(M, Vector(steps));
		RNG rng(1, seed);
		for(long i = 0; i<M; ++i){
			genGaussianMT(Z1Arr[i], rng, 0);
			genGaussianMT(Z1perpArr[i], rng, 0);
			hierarchical2increments(Z1Arr[i], W1Arr[i]);
			hierarchical2increments(Z1perpArr[i], W1perpArr[i]);
		}
		// now call the function to compute all the samples.
		std::vector<Vector> payoffArr = ComputePayoffRTsamples(xi, H, eta, rho, T, K,
				numThreads, W1Arr, W1perpArr);

		// now compute the prices and statistical error estimates
		size_t L = payoffArr[0].size();
		Vector price(L, 0.0);
		Vector stat(L, 0.0);
		for(long i=0; i<M; ++i){
			for(size_t j = 0; j<L; ++j){
				price[j] += payoffArr[i][j];
				stat[j] += payoffArr[i][j] * payoffArr[i][j];
			}
		}
		for(size_t j = 0; j<L; ++j){
			price[j] = price[j] / static_cast<double>(M);
			stat[j] = sqrt(stat[j] / static_cast<double>(M) - price[j]*price[j]) / sqrt(static_cast<double>(M));
		}
		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(price[i] - truePrice[i]) < epsilon + 2*stat[i]);
	}

	/*
	SECTION("Basic testing of hierarchical representation:"){
		// generate the normals
		std::vector<Vector> Z1Arr(M, Vector(steps));
		std::vector<Vector> Z1perpArr(M, Vector(steps));
		std::vector<Vector> W1Arr(M, Vector(steps));
		std::vector<Vector> W1perpArr(M, Vector(steps));
		RNG rng(1, seed);
		for(long i = 0; i<M; ++i){
			genGaussianMT(Z1Arr[i], rng, 0);
			genGaussianMT(Z1perpArr[i], rng, 0);
			hierarchical2increments(Z1Arr[i], W1Arr[i]);
			hierarchical2increments(Z1perpArr[i], W1perpArr[i]);
		}
		// Compute very simple statistics: sample variances of the entries should be one.
		Vector varZ1 = sampleVar(Z1Arr);
		Vector varW1 = sampleVar(W1Arr);


		// Check that variances are close to 1
		double max_diff_Z1 = 0.0;
		double max_diff_W1 = 0.0;
		for(size_t i = 0; i < varZ1.size(); ++i){
			double diff1 = fabs(varZ1[i] - 1.0);
			if(diff1 > max_diff_Z1)
				max_diff_Z1 = diff1;
			double diff2 = fabs(varW1[i] - 1.0);
			if(diff2 > max_diff_W1)
				max_diff_W1 = diff2;
		}

		REQUIRE( max_diff_Z1 < 2.0 * epsilon);
		REQUIRE( max_diff_W1 < 2.0 * epsilon);
	}

	SECTION("Basic testing of hierarchical representation 2:"){
		// generate the normals
		Vector var = hierarchical2incrementsVar(log2(steps));
		// print variances
		//std::cout << "Var before hierarchical transform:\n" << varZ1 << std::endl;
		std::cout << "Theoretical var after hierarchical transform:\n" << var << std::endl;

		// Strange: the actual sample variance does not work, but the theoretical one does. Is there some problem with the input???
		// Or maybe with the computation of the "index" in hierarchical2increments?

		// Try explicit formula
		auto dW = hierarchical2string(log2(steps));
		// plot out the last entry
		std::cout << "dW[" << steps - 1 << "]:\n" << dW[steps-1] << std::endl;

		// Check that variances are close to 1

	}
	*/
}
