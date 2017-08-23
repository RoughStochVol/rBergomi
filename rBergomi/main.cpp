// TODO:
// -> Arbitrary input vectors
// -> Combination of console-input and .txt in put for use from other environments without formal
//	  interface.
// Idea for R integration: Construct C++ library and then call it from Rcpp
// Other idea: go back to old .C-interface: just use an interface based on pointers

// Console and .txt-file input:
// call in the following form
// rBergomi N M num_threads out.txt
// N ... number of timesteps
// M ... number of samples
// num_threads ... number of threads
// path ... path for input and output-files
// out.txt ... name of the output-text-file.
// Further input variables:
// H ... read from local file .H.txt
// eta ... read from local file .eta.txt
// rho ... read from local file .rho.txt
// xi ... read from local file .xi.txt
// T ... read from local file .T.txt
// K ... read from local file .K.txt
// Those files are supposed to only contain the corresponding numbers.

#include <chrono>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <random>
#include "RBergomi.h"
#include "rBergomiMT.h"
#include "BlackScholes.h"
#include "interface.h"
#include "aux.h"

typedef std::vector<double> Vector;

// There are other clocks, but this is usually the one you want.
// It corresponds to CLOCK_MONOTONIC at the syscall level.
using Clock = std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

int main(int argc, char* argv[]){
	if(argc != 6){
		std::cerr << "Error: Wrong number of arguments.\nUsage: rBergomi N M num_threads path out_file\n";
		exit(5);
	}
	// Convert the console input
	long N = atoi(argv[1]);
	long M = atoi(argv[2]);
	int num_threads = atoi(argv[3]);
	std::string path_name(argv[4]);
	std::string fname(argv[5]);

	// Double check: output of console arguments as recorded
	std::cerr << "N = " << N << ", M = " << M << ", num_threads = " << num_threads
			<< ",\noutfile = " << fname << std::endl;

	// read in the further input files.
	const std::string hName = path_name + std::string(".H.txt");
	Vector H = file2vector<double>(hName);
	const std::string etaName = path_name + std::string(".eta.txt");
	Vector eta = file2vector<double>(etaName);
	const std::string rhoName = path_name + std::string(".rho.txt");
	Vector rho = file2vector<double>(rhoName);
	const std::string TName = path_name + std::string(".T.txt");
	Vector T = file2vector<double>(TName);
	const std::string KName = path_name + std::string(".K.txt");
	Vector K = file2vector<double>(KName);
	const std::string xiName = path_name + std::string(".xi.txt");
	Vector xi = file2vector<double>(xiName);

	// Double check that all the inputs have positive size
	if((H.size() <=0) || (eta.size() <=0) || (rho.size() <=0) ||
			(T.size() <=0) || (K.size() <=0) || (xi.size() <=0)){
		std::cerr << "One or more parameter vectors have size 0.\n";
		exit(17);
	}

	// generate seed;
	std::vector<uint64_t> seed(num_threads);
	std::random_device rd;
	for(size_t i=0; i < seed.size(); ++i)
		seed[i] = rd();

	// run the code
	ResultSimple res;
	auto start = Clock::now();
	if(num_threads > 1){
		res = ComputeIVRT_multiThread(xi[0], H, eta, rho, T, K, N, M, num_threads, seed);
	}
	else{
		res = ComputeIVRT_singleThread(xi[0], H, eta, rho, T, K, N, M, seed);
	}
	auto end = Clock::now();
	auto diff = duration_cast<milliseconds>(end - start);

	// print result in output file
	std::string fname_full = path_name + fname;
	std::ofstream file;
	file.open(fname_full.c_str(), std::ofstream::out | std::ofstream::trunc);
	if(!file.is_open()){
		std::cerr << "/nError while opening file " << fname << "./n";
		exit(1);
	}
	file << "H eta rho T K price iv stat\n";
	// set precision
	file << std::setprecision(10);
	for(size_t i=0; i<res.H.size(); ++i)
		file << res.H[i] << " " << res.eta[i] << " " << res.rho[i]
		     << " " << res.T[i] << " " << res.K[i] << " " << res.price[i]
		     << " " << res.iv[i] << " " << res.stat[i] << "\n";
	file.close();

	std::cout << "Time elapsed: " << diff.count() << "ms.\n";

	return 0;
}

//int main() {
//	const double xi = 0.07;
//	const Vector H { 0.07, 0.1 };
//	Vector eta { 2.2, 2.4 };
//	Vector rho { -0.9, -0.85 };
//	Vector T { 0.5, 1.0 };
//	Vector K { 0.8, 1.0 };
//	Vector truePrice { 0.2149403, 0.08461882, 0.2157838, 0.07758564 };
//	long steps = 256;
//	long M = 10 * static_cast<long> (pow(static_cast<double> (steps), 2.0));
//	std::vector<uint64_t> seed { 123, 452, 567, 248, 9436, 675, 194, 6702 };
//	RBergomiST rBergomi(xi, H, eta, rho, T, K, steps, M, seed);
//	const int numThreads = 8; // numThreads = 1 exactly gives the same results as single-threaded version!
//
//	const double epsilon = 0.005;
//
//	// generate the normals
//	std::vector<Vector> W1Arr(M, Vector(steps));
//	std::vector<Vector> W1perpArr(M, Vector(steps));
//	RNG rng(1, seed);
//	for(long i = 0; i<M; ++i){
//		genGaussianMT(W1Arr[i], rng, 0);
//		genGaussianMT(W1perpArr[i], rng, 0);
//	}
//
//	//std::cout << "First entries of W1 and W1perp:\n"
//	//		<< W1Arr[0] << "\n" << W1perpArr[0] << "\n";
//
//	if(numThreads < 1)
//		std::cout << "Single-threaded version.\n";
//	else
//		std::cout << "Multi-threaded version. Number of threads = " << numThreads << "\n";
//
//	// now call the function to compute all the samples.
//	std::vector<Vector> payoffArr;
//	if(numThreads < 1)
//		payoffArr = ComputePayoffRTsamples_ST(xi, H, eta, rho, T, K, W1Arr, W1perpArr);
//	else
//		payoffArr = ComputePayoffRTsamples(xi, H, eta, rho, T, K, numThreads, W1Arr, W1perpArr);
//
//	// now compute the prices and statistical error estimates
//	size_t L = payoffArr[0].size();
//	//std::cout << "Dimensions of payoffArr = (" << payoffArr.size() << ", " << payoffArr[0].size() << ")\n";
//	//std::cout << "First entry of payoffArr = (" << payoffArr[0][0]
//	         // << ", " << payoffArr[0][1] << ", " << payoffArr[0][2]<< ", " << payoffArr[0][3]
//	//          << ")\n";
//	Vector price(L, 0.0);
//	Vector stat(L, 0.0);
//	//std::cout << "Length of price = " << price.size() << ", first entry = " << price[0] << "\n";
//	for(size_t j = 0; j<L; ++j){
//		for(long i=0; i<M; ++i){
//			price[j] += payoffArr[i][j];
//			stat[j] += payoffArr[i][j] * payoffArr[i][j];
//		}
//		price[j] = price[j] / static_cast<double>(M);
//		stat[j] = sqrt(stat[j] / static_cast<double>(M) - price[j]*price[j]) / sqrt(static_cast<double>(M));
//	}
//	for (int i = 0; i < price.size(); ++i)
//		std::cout << fabs(price[i] - truePrice[i]) << " < " << epsilon + 2*stat[i] << std::endl;
//}
