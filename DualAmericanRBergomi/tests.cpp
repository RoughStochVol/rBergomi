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
#include "RBergomi.h"

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

		// Compare against reference Haar values (from R).
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
		//std::cout << "ds = " << T / N << ", threshold = " << confidenceFactor/sqrt(M) << std::endl;
		for(int n1=0; n1<N; ++n1){
			for(int n2=0; n2<N; ++n2){
				double trueVar = (n1 == n2) ? T / N : 0.0;
				REQUIRE(fabs(cov1[n1][n2] - trueVar) < confidenceFactor/sqrt(M));
				REQUIRE(fabs(cov2[n1][n2] - trueVar) < confidenceFactor/sqrt(M));
				REQUIRE(fabs(cov12[n1][n2] - 0.0) < confidenceFactor/sqrt(M));
			}
		}
	}

	SECTION("Test std::partial_sum:"){
		Vector x {0.1, 0.2, 0.1};
		Vector xSum(x.size()+1, 0.0);
		std::partial_sum(x.begin(), x.end(), xSum.begin() + 1);
		Vector xSumRef { 0.0, 0.1, 0.3, 0.4};
		for(size_t i=0; i<xSum.size(); ++i)
			REQUIRE(fabs(xSum[i] - xSumRef[i]) < 0.000001);
	}

	SECTION("Test the integral of the fractional kernel against the Haar basis:"){
		// Compare against reference Haar values from R. Note that these values seem to be inaccurate!
		Vector PhiTildeRef {0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2715394, 0.4031066, 0.507916, 0.5984212, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0.2715394, 0.4031066, 0.507916, 0.5984212, 0.4080486, 0.3509067, 0.3153471, 0.2899492, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2715394, 0.4031066, 0.507916, 0.5984212, 0.4080486, 0.3509067, 0.3153471, 0.2899492,
			0.270472, 0.2548511, 0.2419257, 0.2309801, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2715394,
			0.4031066, -0.03516271, -0.2077919, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2715394, 0.4031066, -0.03516271,
			-0.2077919, -0.06470468, -0.03972256, -0.02799681, -0.02123496, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2715394,
			0.4031066, -0.03516271, -0.2077919, -0.06470468, -0.03972256, -0.02799681, -0.02123496, -0.01687833, -0.01386313,
			-0.01166793, -0.01000785, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3840146, -0.1979505,
			-0.03784122, -0.02022923, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3840146, -0.1979505, -0.03784122,
			-0.02022923, -0.01320655, -0.009533861, -0.007319203, -0.005858502, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3840146,
			-0.1979505, -0.03784122, -0.02022923, -0.01320655, -0.009533861, -0.007319203, -0.005858502, -0.004833354, -0.004080223,
			-0.003507141, -0.003058736, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3840146, -0.1979505,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3840146, -0.1979505, -0.03784122, -0.02022923, -0.01320655,
			-0.009533861, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3840146, -0.1979505, -0.03784122, -0.02022923, -0.01320655,
			-0.009533861, -0.007319203, -0.005858502, -0.004833354, -0.004080223, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};


		//  Compare against reference Haar values (from Mathematica). These values essentially coincide with the C++ ones, up to order
		// e-07 at least.
		Vector PhiTildeMath {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.271539, 0.403107, 0.507916,
			0.598421, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
			0., 0., 0.0000695414, 0.271539, 0.403107, 0.507916, 0.598421,
			0.408049, 0.350907, 0.315347, 0.289949, 0., 0., 0., 0., 0., 0., 0.,
			0., 0., 0., 0., 0., 0., 0., 0.271539, 0.403107, 0.507916, 0.598352,
			0.408049, 0.350907, 0.315347, 0.289949, 0.270472, 0.254851, 0.241926,
			0.23098, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
			0., 0., 0., 0., 0., 0., 0., 0.271539, 0.402968, -0.0351626,
			-0.207792, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
			0., 0., 0., 0.0000695414, 0.271539, 0.403107, -0.0351628, -0.207792,
			-0.0647047, -0.0397226, -0.0279968, -0.021235, 0., 0., 0., 0., 0.,
			0., 0., 0., 0., 0., 0., 0., 0., 0., 0.271539, 0.403107, -0.0351627,
			-0.207722, -0.0647047, -0.0397226, -0.0279968, -0.021235, -0.0168783,
			-0.0138631, -0.0116679, -0.0100079, 0., 0., 0., 0., 0., 0., 0., 0.,
			0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.384015,
			-0.197852, -0.0378412, -0.0202292, 0., 0., 0., 0., 0., 0., 0., 0.,
			0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0000983463, 0.384015, -0.19795,
			-0.0378412, -0.0202292, -0.0132065, -0.00953386, -0.0073192,
			-0.0058585, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
			0.383962, -0.19795, -0.0378412, -0.0202292, -0.0132066, -0.00953386,
			-0.0073192, -0.0058585, -0.00483335, -0.00408022, -0.00350714,
			-0.00305874, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
			0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0000983459, 0.384015, -0.19795,
			0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
			0., 0., 0., 0.383818, -0.19795, -0.0378412, -0.0202292, -0.0132066,
			-0.00953386, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
			0., 0., 0.384015, -0.197852, -0.0378412, -0.0202292, -0.0132065,
			-0.00953386, -0.0073192, -0.0058585, -0.00483335, -0.00408022, 0.,
			0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

		// compute values
		const Vector tGrid {0, 0.06666667, 0.1333333, 0.2, 0.2666667, 0.3333333, 0.4, 0.4666667, 0.5333333, 0.6,
			0.6666667, 0.7333333, 0.8};
		Vector PhiTildeVal(PhiTildeRef.size());
		int count = 0;
		const int iotaMax = 2 * N;
		for(int iota = 0; iota < iotaMax; ++iota){
			for(size_t ti=0; ti < tGrid.size(); ++ti){
				PhiTildeVal[count] = haar.PhiTilde(tGrid[ti], iota);
				count++;
			}
		}

		for(size_t i = 0; i < PhiTildeVal.size(); ++i){
//			if((fabs(PhiTildeVal[i] - PhiTildeRef[i]) > 0.000001) || (fabs(PhiTildeVal[i] - PhiTildeMath[i]) > 0.000001) ||
//					(fabs(PhiTildeMath[i] - PhiTildeRef[i]) > 0.000001))
//				std::cout << "i = " << i << ", diff(C++, R) = " << fabs(PhiTildeVal[i] - PhiTildeRef[i])
//				<< ", diff(C++,Mathematica) = " << fabs(PhiTildeVal[i] - PhiTildeMath[i])
//				<< ", diff(R,Mathematica) = " << fabs(PhiTildeRef[i] - PhiTildeMath[i])
//				<< std::endl;
			REQUIRE(fabs(PhiTildeVal[i] - PhiTildeMath[i]) < 0.00001);
		}

	}

	SECTION("Test generation of fBm together with Bm:"){
		const int M = 40000; // number of samples
		std::vector<Vector> dB1(M, Vector(N));
		std::vector<Vector> dB2(M, Vector(N));
		std::vector<Vector> Wtilde(M, Vector(N));
		// compute B1 and B2 as well; includes the Bi[0], but not the final term.
		std::vector<Vector> B1(M, Vector(N, 0.0));
		std::vector<Vector> B2(M, Vector(N, 0.0));
		for(int m = 0; m < M; ++m){
			haar.generate();
			dB1[m] = haar.dB(0);
			dB2[m] = haar.dB(1);
			Wtilde[m] = haar.Wtilde();
			std::partial_sum(dB1[m].begin(), dB1[m].end() - 1, B1[m].begin() + 1);
			std::partial_sum(dB2[m].begin(), dB2[m].end() - 1, B2[m].begin() + 1);
		}

		// compute the mean and covariance matrix of Wtilde and the covariance between Wtilde and Bi
		auto meanW = mean(Wtilde);
		auto covW = covariance(Wtilde, Wtilde);
		auto covWB1 = covariance(Wtilde, B1);
		auto covWB2 = covariance(Wtilde, B2);

		std::vector<Vector> covWMath {{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0.294894,
				  0.142913, 0.113842, 0.098304, 0.0881606, 0.0808381, 0.0752166,
				  0.070717, 0.0670048, 0.0638711, 0.0611774, 0.0588282}, {0.,
				  0.142913, 0.364174, 0.198061, 0.16147, 0.141018, 0.127328, 0.117282,
				   0.10948, 0.103182, 0.0979515, 0.0935127, 0.0896808}, {0., 0.113842,
				   0.198061, 0.408114, 0.236004, 0.195498, 0.17222, 0.15636, 0.144577,
				   0.135343, 0.127835, 0.121565, 0.116219}, {0., 0.098304, 0.16147,
				  0.236004, 0.440728, 0.265427, 0.222449, 0.197294, 0.179933,
				  0.166912, 0.156633, 0.148227, 0.141173}, {0., 0.0881606, 0.141018,
				  0.195498, 0.265427, 0.467231, 0.289551, 0.244924, 0.218428,
				  0.199958, 0.186002, 0.174918, 0.16581}, {0., 0.0808381, 0.127328,
				  0.17222, 0.222449, 0.289551, 0.489388, 0.310167, 0.264307, 0.236793,
				   0.217465, 0.20277, 0.191042}, {0., 0.0752166, 0.117282, 0.15636,
				  0.197294, 0.244924, 0.310167, 0.508419, 0.328241, 0.281402,
				  0.253086, 0.233069, 0.217774}, {0., 0.070717, 0.10948, 0.144577,
				  0.179933, 0.218428, 0.264307, 0.328241, 0.525526, 0.344268,
				  0.296709, 0.267751, 0.247172}, {0., 0.0670048, 0.103182, 0.135343,
				  0.166912, 0.199958, 0.236793, 0.281402, 0.344268, 0.540748,
				  0.358778, 0.310608, 0.281116}, {0., 0.0638711, 0.0979515, 0.127835,
				  0.156633, 0.186002, 0.217465, 0.253086, 0.296709, 0.358778,
				  0.554429, 0.372066, 0.323355}, {0., 0.0611774, 0.0935127, 0.121565,
				  0.148227, 0.174918, 0.20277, 0.233069, 0.267751, 0.310608, 0.372066,
				   0.56727, 0.384229}, {0., 0.0588282, 0.0896808, 0.116219, 0.141173,
				  0.16581, 0.191042, 0.217774, 0.247172, 0.281116, 0.323355, 0.384229,
				   0.579004}};

		std::vector<Vector> covMixedMath {{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0.140213,
		  0.140222, 0.140222, 0.140222, 0.140222, 0.140222, 0.140222,
		  0.140222, 0.140222, 0.140222, 0.140222, 0.140222}, {0., 0.067941,
		  0.208163, 0.208163, 0.208163, 0.208163, 0.208163, 0.208163,
		  0.208163, 0.208163, 0.208163, 0.208163, 0.208163}, {0., 0.0541234,
		  0.122064, 0.262287, 0.262287, 0.262287, 0.262287, 0.262287,
		  0.262287, 0.262287, 0.262287, 0.262287, 0.262287}, {0., 0.0467367,
		  0.10086, 0.168801, 0.308987, 0.309023, 0.309023, 0.309023, 0.309023,
		   0.309023, 0.309023, 0.309023, 0.309023}, {0., 0.0419143, 0.088651,
		  0.142774, 0.210715, 0.350938, 0.350938, 0.350938, 0.350938,
		  0.350938, 0.350938, 0.350938, 0.350938}, {0., 0.0384331, 0.0803474,
		  0.127084, 0.181207, 0.249148, 0.389371, 0.389371, 0.389371,
		  0.389371, 0.389371, 0.389371, 0.389371}, {0., 0.0357605, 0.0741935,
		  0.116108, 0.162845, 0.216968, 0.284909, 0.425095, 0.425131,
		  0.425131, 0.425131, 0.425131, 0.425131}, {0., 0.0336212, 0.0693817,
		  0.107815, 0.149729, 0.196466, 0.250589, 0.31853, 0.458752, 0.458752,
		   0.458752, 0.458752, 0.458752}, {0., 0.0318564, 0.0654776, 0.101238,
		   0.139671, 0.181585, 0.228322, 0.282446, 0.350386, 0.490609,
		  0.490609, 0.490609, 0.490609}, {0., 0.0303665, 0.0622228, 0.0958441,
		   0.131605, 0.170038, 0.211952, 0.258689, 0.312812, 0.380753,
		  0.520939, 0.520975, 0.520975}, {0., 0.0290858, 0.0594523, 0.0913087,
		   0.12493, 0.16069, 0.199123, 0.241038, 0.287774, 0.341898, 0.409839,
		   0.550061, 0.550061}, {0., 0.0279689, 0.0570548, 0.0874213,
		  0.119278, 0.152899, 0.188659, 0.227092, 0.269007, 0.315743,
		  0.369867, 0.437808, 0.57803}};

		// meanW must be 0
		for(auto& m : meanW)
			REQUIRE(fabs(m) < confidenceFactor/sqrt(M));

		// NOTE THAT THERE IS A BIAS HERE!!!
		// Check covariance
		const double ds = T/N;
		for(size_t n1=0; n1<covW.size(); ++n1){
			for(size_t n2=0; n2<covW[0].size(); ++n2){
				double temp = covWMath[n1][n2];//covWtilde(n1 * ds, n2 * ds, H);
				REQUIRE(fabs(covW[n1][n2] - temp) < confidenceFactor/sqrt(M));
				temp = covMixedMath[n1][n2];//covCross(n1 * ds, n2 * ds, H, 0);
				REQUIRE(fabs(covWB1[n1][n2] - temp) < confidenceFactor/sqrt(M));
				temp = covCross(n1 * ds, n2 * ds, H, 1);
				REQUIRE(fabs(covWB2[n1][n2] - temp) < confidenceFactor/sqrt(M));
			}
		}
	}


}

TEST_CASE( "Test the rBergomi simulation with H = 0.5.", "[Bergomi]" ) {
	// model parameters
	const double H = 0.5; // Sanity check: test standard Bm!
	const double rho = -0.9;
	const double eta = 1.9;
	const std::function<double(double)> xi0 = [](double t) -> double {return 0.04;};
	const double S0 = 1.0;
	const double r = 0.0;

	// option parameters
	const int K = 3;
	const double T = 0.8;
	const double strike = 0.95;
	const double EuropRefPrice = 0.09104141;

	// numerical parameters
	const int N = 4 * K;
	std::vector<uint64_t> seed { 1, 2 };
	const double confidenceFactor = 4.5; // corresponds to a confidence set at level 3e-06
	const int M = 100000;

	// initialize classes
	HaarRep2D haar(H, T, K, N, seed);
	RBergomi rBergomi(&haar, S0, r, eta, rho, xi0);


	SECTION("Test the vectorized exponential function:"){
		Vector x{1.0, 3.4, 1.7, -0.9};
		Vector y(x.size());
		for(size_t i=0; i<x.size(); ++i)
			y[i] = exp(x[i]);
		Vector z(x.size());
		expVector(x, z);
		for(size_t i=0; i<x.size(); ++i)
			REQUIRE(fabs(y[i] - z[i]) < 0.0000001);
	}

	SECTION("Test that S is a martingale:"){
		std::vector<Vector> S(M, Vector(N+1));
		for(int m=0; m<M; ++m)
			S[m] = rBergomi.generate();
		auto mS = mean(S);
		//std::cout << "mean(S) = " << mS << std::endl;
		for(int n=0; n<N+1; ++n)
			REQUIRE(fabs(mS[n] - S0) < confidenceFactor/sqrt(M));
	}

	SECTION("Test for the expectation of v:"){
		std::vector<Vector> v(M, Vector(N, 0.0));
		for(int m=0; m<M; ++m)
			v[m] = std::get<1>(rBergomi.generate_SV());
		auto mv = mean(v);
		for(auto vv : mv)
			REQUIRE(fabs(vv - xi0(0.0)) < confidenceFactor/sqrt(M));// note that even the theoretical values are far off!
	}

	SECTION("Test the European option price:"){
		Vector payoffs(M);
		for(long m=0; m<M; ++m)
			payoffs[m] = rBergomi.generate().back(); // the asset price at the terminal time
		// now compute the payoffs
		std::transform(payoffs.begin(), payoffs.end(), payoffs.begin(),
				[strike](double s){return s > strike ? s - strike : 0.0;});
		// compute mean and variance of payoffs
		double price = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / M;
		double stat = sqrt(std::inner_product(payoffs.begin(), payoffs.end(), payoffs.begin(), 0.0) / M -
						price * price) / sqrt(M);

		// From other experiments, expect the error to be around 10 %
		// Note: this is a problematic condition!!
		double error = fabs(price - EuropRefPrice);
		REQUIRE(error < 0.01);
	}
}

TEST_CASE( "Test the rBergomi simulation with H = 0.3.", "[rBergomi]" ) {
	// model parameters
	const double H = 0.3; // Sanity check: test standard Bm!
	const double rho = -0.9;
	const double eta = 2.5;
	const std::function<double(double)> xi0 = [](double t) -> double {return 0.04;};
	const double S0 = 1.0;
	const double r = 0.0;

	// option parameters
	const int K = 3;
	const double T = 0.8;

	// numerical parameters
	const int N = 4 * K;
	std::vector<uint64_t> seed { 1, 2 };
	const double confidenceFactor = 4.5; // corresponds to a confidence set at level 3e-06
	const int M = 100000;

	// initialize classes
	HaarRep2D haar(H, T, K, N, seed);
	RBergomi rBergomi(&haar, S0, r, eta, rho, xi0);


	SECTION("Test the vectorized exponential function:"){
		Vector x{1.0, 3.4, 1.7, -0.9};
		Vector y(x.size());
		for(size_t i=0; i<x.size(); ++i)
			y[i] = exp(x[i]);
		Vector z(x.size());
		expVector(x, z);
		for(size_t i=0; i<x.size(); ++i)
			REQUIRE(fabs(y[i] - z[i]) < 0.0000001);
	}

	SECTION("Test that S is a martingale:"){
		std::vector<Vector> S(M, Vector(N+1));
		for(int m=0; m<M; ++m)
			S[m] = rBergomi.generate();
		auto mS = mean(S);
		//std::cout << "mean(S) = " << mS << std::endl;
		for(int n=0; n<N+1; ++n)
			REQUIRE(fabs(mS[n] - S0) < confidenceFactor/sqrt(M));
	}

	SECTION("Test for the expectation of v:"){
		std::vector<Vector> v(M, Vector(N, 0.0));
		for(int m=0; m<M; ++m)
			v[m] = std::get<1>(rBergomi.generate_SV());
		auto mv = mean(v);
		for(auto vv : mv)
			REQUIRE(fabs(vv - xi0(0.0)) < confidenceFactor/sqrt(M));// note that even the theoretical values are far off!
	}
}

// Note:
// For both H = 0.5 and H = 0.3, the martingality of S is satisfied with similar precision (i.e.,
// mean(S) - 1 is of similar size for both Hs.
