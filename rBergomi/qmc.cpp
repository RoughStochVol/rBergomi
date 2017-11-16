/*
 * qmc.cpp
 *
 *  Created on: Nov 3, 2017
 *      Author: bayerc
 */

#include"qmc.h"

std::vector<Vector> normalQMC(long M, int N){
	// Check N <= 1111
	if(N > maxDim){
		std::cerr << "The dimension N = " << N << " is larger than the maximum dimension " << maxDim
				<< " for Sobol' numbers\n";
		exit(24);
	}
	// Allocate memory.
	std::vector<Vector> X(M, Vector(N));
	for(long seed = 0; seed < M; ++seed)
		normalQMC_sample(X[seed], N, seed + 1); // start with 1 to avoid problems at inversion.

	return X;
}

void normalQMC_sample(Vector& Z, int N, long long i){
	i8_sobol(N, &i, &(Z[0])); // uniform number. Now transform by normal quantile.
	// Check that none of the uniform numbers is 0.
	if(minVec(Z) < 0.000000000000001){
		std::cerr << "Z = " << Z << " hits the singularity of the normal quantile function.\n"
				<< "Z was generated with i = " << i << ".\n";
		long long ii = 1;
		i8_sobol(N, &ii, &(Z[0]));
		std::cerr << "Z = " << Z << " hits the singularity of the normal quantile function.\n"
						<< "Z was generated with ii = " << ii << ".\n";
		exit(33);
	}
	Phi_inv(Z); // done!
}
