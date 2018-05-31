/*
 * tests.cpp
 *
 * Unit tests for RBergomi and the multi-threaded version.
 *
 *  Created on: Apr 13, 2017
 *      Author: bayerc
 */

//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <algorithm>
#include"catch.hpp"
#include "HaarRep2D.h"

TEST_CASE( "Test the hierarchical representation in the Haar sense.", "[hierarchical]" ) {
	const double H = 0.07;
	const double T = 0.8;
	const int K = 3; // # of exercise dates
	const int N = K * 4; // # of grid points on fine grid
	std::vector<uint64_t> seed { 1, 2 };
	const double confidenceFactor = 4.5; // corresponds to a confidence set at level 3e-06

	HaarRep2D haar(H, T, K, N, seed);

	SECTION("Test the multi-index notation:"){
		MultiIndex alpha, beta;
		const int testSize = 13;
		std::vector<MultiIndex> ref(testSize);
		ref[0] = MultiIndex {0, 0, 0, 2};
		ref[1] = MultiIndex {0, 0, 1, 2};
		ref[2] = MultiIndex {0, 0, 0, 1};
		ref[3] = MultiIndex {0, 0, 1, 1};
		ref[4] = MultiIndex {0, 0, 0, 0};
		ref[5] = MultiIndex {0, 0, 1, 0};
		ref[6] = MultiIndex {1, 0, 0, 2};
		ref[7] = MultiIndex {1, 0, 1, 2};
		ref[8] = MultiIndex {1, 0, 0, 1};
		ref[9] = MultiIndex {1, 0, 1, 1};
		ref[10] = MultiIndex {1, 0, 0, 0};
		ref[11] = MultiIndex {1, 0, 1, 0};
		ref[12] = MultiIndex {2, 0, 0, 2};
		int kappa;
		for(int iota = 0; iota < testSize; ++iota) {
			alpha = ref[iota];
			beta = haar.iota2alpha(iota);
			kappa = haar.alpha2iota(alpha);
			//std::cout << "alpha = " << alpha << ", beta = " << beta << std::endl;
			REQUIRE(alpha == beta);
			REQUIRE(iota == kappa);
		}
	}

	SECTION("Test the integrated Haar basis functions:"){
		// evaluate each Haar function for the different values of l on the grid.

		// Compare against reference Haar values.
		Vector haarRef {0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1290994, 0.2581989, 0.3872983, 0.5163978, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1290994, 0.2581989, 0.3872983, 0.5163978, 0, 0,
			0, 0, 0, 0.1290994, 0.2581989, 0.3872983, 0.5163978, 0.5163978, 0.5163978, 0.5163978, 0.5163978, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1290994, 0.2581989, 0.3872983, 0.5163978,
			0.5163978, 0.5163978, 0.5163978, 0.5163978, 0, 0.1290994, 0.2581989, 0.3872983, 0.5163978, 0.5163978, 0.5163978,
			0.5163978, 0.5163978, 0.5163978, 0.5163978, 0.5163978, 0.5163978, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1290994, 0.2581989, 0.3872983, 0.5163978, 0.5163978, 0.5163978, 0.5163978, 0.5163978,
			0.5163978, 0.5163978, 0.5163978, 0.5163978, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1290994, 0.2581989, 0.1290994, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1290994, 0.2581989,
			0.1290994, 0, 0, 0, 0, 0, 0, 0.1290994, 0.2581989, 0.1290994, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1290994, 0.2581989, 0.1290994, 0, 0, 0, 0, 0, 0, 0.1290994, 0.2581989,
			0.1290994, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0.1290994, 0.2581989, 0.1290994, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1825742, 0, -3.252607e-18, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1825742, 0,
			-3.252607e-18, 0, 0, 0, 0, 0, 0, 0.1825742, 0, -9.486769e-18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1825742, 0, -9.486769e-18, 0, 0, 0, 0, 0, 0, 0.1825742, 0, -9.486769e-18,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1825742, 0,
			-9.486769e-18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1825742, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1825742, 0, 0, 0, 0, 0, 0, 0, 0, 0.1825742,
			8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1825742, 8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17,
			0, 0, 0, 0.1825742, 8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17,
			8.326673e-17, 8.326673e-17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1825742,
			8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17, 8.326673e-17};

		// compute values
		const Vector tGrid {0, 0.06666667, 0.1333333, 0.2, 0.2666667, 0.3333333, 0.4, 0.4666667, 0.5333333, 0.6,
			0.6666667, 0.7333333, 0.8};
		Vector haarVal(haarRef.size());
		int count = 0;
		const int iotaMax = 2 * N;
		for(int iota = 0; iota < iotaMax; ++iota){
			for(size_t ti=0; ti < tGrid.size(); ++ti){
				haarVal[count] = haar.Phi(tGrid[ti], iota, 0);
				count++;
			}
			for(size_t ti=0; ti < tGrid.size(); ++ti){
				haarVal[count] = haar.Phi(tGrid[ti], iota, 1);
				count++;
			}
		}

//		std::vector<int> haarDiff;
//		for(size_t i=0; i<haarVal.size(); ++i)
//			if(fabs(haarVal[i] - haarRef[i]) > 0.001)
//				haarDiff.push_back(i);
//
////		std::cout << "haarVal = " << haarVal << std::endl;
////		std::cout << "haarRef = " << haarRef << std::endl;
//		std::cout << "haarDiff = " << haarDiff << std::endl;
//		std::cout << "haarVal[108] = " << haarVal[108] << ", haarRef[108] = " << haarRef[108] << std::endl;
//		std::cout << "iota = 4, t = 0.2666667, i = 0, Phi = " << haar.Phi(0.2666667, 4, 0) << std::endl;
		// compare
		for(size_t i = 0; i < haarVal.size(); ++i)
			REQUIRE(fabs(haarVal[i] - haarRef[i]) < 0.000001);
	}

	SECTION("Test mean and covariance functions:"){
		std::vector<Vector> X(2);
		X[0] = Vector{1.0, 2.0};
		X[1] = Vector{3.0, 4.0};
		std::vector<Vector> Y(2);
		Y[0] = Vector{2.0, 1.0};
		Y[1] = Vector{4.0, 3.0};
		Vector meanX{2.0, 3.0};
		std::vector<Vector> meanXY(2);
		meanXY[0] = Vector{7.0, 5.0};
		meanXY[1] = Vector{10.0, 7.0};
		auto mX = mean(X);
		auto mXY = covariance(X, Y);
		REQUIRE(fabs(meanX[0] - mX[0]) < 0.00001);
		REQUIRE(fabs(meanX[1] - mX[1]) < 0.00001);
		REQUIRE(fabs(meanXY[0][0] - mXY[0][0]) < 0.00001);
		REQUIRE(fabs(meanXY[0][1] - mXY[0][1]) < 0.00001);
		REQUIRE(fabs(meanXY[1][0] - mXY[1][0]) < 0.00001);
		REQUIRE(fabs(meanXY[1][1] - mXY[1][1]) < 0.00001);
	}

	SECTION("Test generation of standard normal vector X:"){
		const int M = 40000;
		const int I = N*2;
		std::vector<Vector> X(M, Vector(I));
		// check the internal value of I
		haar.generate();
		auto temp = haar.getX();
		REQUIRE(temp.size() == I);

		for(int m=0; m<M; ++m){
			haar.generate();
			X[m] = haar.getX();
		}

		// check first and second moments
		auto meanX = mean(X);
		auto covX = covariance(X, X);
		// subtract 1 at the diagonal
		for(int i = 0; i<I; ++i)
			covX[i][i] -= 1.0;
		for(auto& x : meanX)
			REQUIRE(fabs(x) < confidenceFactor / sqrt(M));
		for(auto& v : covX)
			for(auto& w : v)
				REQUIRE(fabs(w) < confidenceFactor / sqrt(M));
	}

	SECTION("Test generation of Brownian motion:"){
		const int M = 40000; // number of samples
		std::vector<Vector> dB1(M, Vector(N));
		std::vector<Vector> dB2(M, Vector(N));
		for(int m = 0; m < M; ++m){
			haar.generate();
			dB1[m] = haar.dB(0);
			dB2[m] = haar.dB(1);
		}

		// first compute means
		auto mean1 = mean(dB1);
		auto mean2 = mean(dB2);
		for(int n=0; n<N; ++n){
			REQUIRE(fabs(mean1[n]) < confidenceFactor/sqrt(M));
			REQUIRE(fabs(mean2[n]) < confidenceFactor/sqrt(M));
		}

		// Next compute covariances of the individual terms
		auto cov1 = covariance(dB1, dB1);
		auto cov2 = covariance(dB2, dB2);
		auto cov12 = covariance(dB1, dB2);

		// Compare with true values
		std::cout << "ds = " << T / N << ", threshold = " << confidenceFactor/sqrt(M) << std::endl;
		for(int n1=0; n1<N; ++n1){
			for(int n2=0; n2<N; ++n2){
				double trueVar = (n1 == n2) ? T / N : 0.0;
//				if((fabs(cov1[n1][n2] - trueVar) > confidenceFactor/sqrt(M)) || (fabs(cov2[n1][n2] - trueVar) > confidenceFactor/sqrt(M))
//						|| (fabs(cov12[n1][n2] - 0.0) > confidenceFactor/sqrt(M))){
//				std::cout << "cov(dB1[" << n1 << "], dB1[" << n2 << "]) = " << cov1[n1][n2]
//						<< "\ncov(dB2[" << n1 << "], dB2[" << n2 << "]) = " << cov2[n1][n2]
//						<< "\ncov(dB1[" << n1 << "], dB2[" << n2 << "]) = " << cov12[n1][n2]
//						<< "\ntrueVar = " << trueVar << "\n========================================"
//																						 << std::endl;
//				}
				REQUIRE(fabs(cov1[n1][n2] - trueVar) < confidenceFactor/sqrt(M));
				REQUIRE(fabs(cov2[n1][n2] - trueVar) < confidenceFactor/sqrt(M));
				REQUIRE(fabs(cov12[n1][n2] - 0.0) < confidenceFactor/sqrt(M));
			}
		}
	}

}
