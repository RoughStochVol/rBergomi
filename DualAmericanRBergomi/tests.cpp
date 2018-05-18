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
#include "HaarRep2D.h"

TEST_CASE( "Test the hierarchical representation in the Haar sense.", "[hierarchical]" )
{
	const double H = 0.07;
	const double T = 0.8;
	const int K = 3; // # of exercise dates
	const int N = K * 4; // # of grid points on fine grid
	std::vector<uint64_t> seed {1, 2};

	HaarRep2D haar(H, T, K, N, seed);

	SECTION("Test the multi-index notation:") {

	}

}
