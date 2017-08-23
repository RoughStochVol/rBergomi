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
#include "interface.h"

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
	long M = 10 * static_cast<long> (pow(static_cast<double> (steps), 2.0));
	std::vector<uint64_t> seed { 123, 452, 567, 248, 9436, 675, 194, 6702 };
	const int numThreads = 8;

	const double epsilon = 0.005;

	SECTION("Single-threaded pricing:") {
		// Check that the price is within epsilon from the confidence interval around the true price
		//Result res = rBergomi.ComputePrice();
		ResultSimple res = ComputeIV_singleThread( xi, H, eta, rho, T, K, steps, M, seed);
		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(res.price[i] - truePrice[i]) < epsilon + res.stat[i]);

		Vector trueVol(truePrice.size());
		trueVol[0] = IV_call(truePrice[0], 1.0, K[0], T[0]);
		trueVol[1] = IV_call(truePrice[1], 1.0, K[1], T[1]);
		trueVol[2] = IV_call(truePrice[2], 1.0, K[0], T[0]);
		trueVol[3] = IV_call(truePrice[3], 1.0, K[1], T[1]);

		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(res.iv[i] - trueVol[i]) < epsilon);
	}

	SECTION("Multi-threaded pricing:") {
		ResultSimple res = ComputeIV_multiThread(xi, H, eta, rho, T, K, steps, M, numThreads,
				seed);

		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(res.price[i] - truePrice[i]) < epsilon + res.stat[i]);

		Vector trueVol(truePrice.size());
		trueVol[0] = IV_call(truePrice[0], 1.0, K[0], T[0]);
		trueVol[1] = IV_call(truePrice[1], 1.0, K[1], T[1]);
		trueVol[2] = IV_call(truePrice[2], 1.0, K[0], T[0]);
		trueVol[3] = IV_call(truePrice[3], 1.0, K[1], T[1]);

		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(res.iv[i] - trueVol[i]) < epsilon);
	}

	SECTION("Multi-threaded pricing with Romano-Touzi trick:") {
		ResultSimple res = ComputeIVRT_multiThread(xi, H, eta, rho, T, K, steps, M, numThreads,
				seed);

		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(res.price[i] - truePrice[i]) < epsilon + res.stat[i]);

		Vector trueVol(truePrice.size());
		trueVol[0] = IV_call(truePrice[0], 1.0, K[0], T[0]);
		trueVol[1] = IV_call(truePrice[1], 1.0, K[1], T[1]);
		trueVol[2] = IV_call(truePrice[2], 1.0, K[0], T[0]);
		trueVol[3] = IV_call(truePrice[3], 1.0, K[1], T[1]);

		for (int i = 0; i < 4; ++i)
			REQUIRE(fabs(res.iv[i] - trueVol[i]) < epsilon);
	}

}
