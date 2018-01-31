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
//
// If no input is given, use default values for debugging, profiling,...

#include <chrono>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include "RBergomi.h"
//#include "BlackScholes.h"
#include "aux.h"

typedef std::vector<double> Vector;

// There are other clocks, but this is usually the one you want.
// It corresponds to CLOCK_MONOTONIC at the syscall level.
using Clock = std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

int main(int argc, char* argv[]) {
	/*
	if ((argc != 6) && (argc != 1)) {
		std::cerr
				<< "Error: Wrong number of arguments.\nUsage: rBergomi N M num_threads path out_file\n";
		exit(5);
	}

	long N, M;
	int num_threads;
	Vector H, eta, rho, T, K, xi;
	std::string path_name, fname;

	if (argc == 1) {
		N = 100;
		M = 100000;
		num_threads = 8;
		H = Vector { 0.05, 0.2 };
		eta = Vector { 1.0, 3.0 };
		rho = Vector { -0.98, -0.8 };
		T = Vector { 0.05, 2.0 };
		K = Vector { 1.0, 1.3 };
		xi = Vector { 0.04 };
		path_name = std::string("./");
		fname = std::string("out.txt");
	} else {
		// Convert the console input
		N = atoi(argv[1]);
		M = atoi(argv[2]);
		num_threads = atoi(argv[3]);
		path_name = std::string(argv[4]);
		fname = std::string(argv[5]);

		// Double check: output of console arguments as recorded
		std::cerr << "N = " << N << ", M = " << M << ", num_threads = "
				<< num_threads << ",\noutfile = " << fname << std::endl;

		// read in the further input files.
		const std::string hName = path_name + std::string(".H.txt");
		H = file2vector<double> (hName);
		const std::string etaName = path_name + std::string(".eta.txt");
		eta = file2vector<double> (etaName);
		const std::string rhoName = path_name + std::string(".rho.txt");
		rho = file2vector<double> (rhoName);
		const std::string TName = path_name + std::string(".T.txt");
		T = file2vector<double> (TName);
		const std::string KName = path_name + std::string(".K.txt");
		K = file2vector<double> (KName);
		const std::string xiName = path_name + std::string(".xi.txt");
		xi = file2vector<double> (xiName);
	}
	// Double check that all the inputs have positive size
	if ((H.size() <= 0) || (eta.size() <= 0) || (rho.size() <= 0) || (T.size()
			<= 0) || (K.size() <= 0) || (xi.size() <= 0)) {
		std::cerr << "One or more parameter vectors have size 0.\n";
		exit(17);
	}

	// generate seed;
	std::vector<uint64_t> seed(num_threads);
	std::random_device rd;
	for (size_t i = 0; i < seed.size(); ++i)
		seed[i] = rd();

	// run the code
	ResultSimple res;
	auto start = Clock::now();
	if (num_threads > 1) {
		res = ComputeIVRT_multiThread(xi[0], H, eta, rho, T, K, N, M,
				num_threads, seed);
	} else {
		res = ComputeIVRT_singleThread(xi[0], H, eta, rho, T, K, N, M, seed);
	}
	auto end = Clock::now();
	auto diff = duration_cast<milliseconds> (end - start);

	// print result in output file
	std::string fname_full = path_name + fname;
	std::ofstream file;
	file.open(fname_full.c_str(), std::ofstream::out | std::ofstream::trunc);
	if (!file.is_open()) {
		std::cerr << "/nError while opening file " << fname << "./n";
		exit(1);
	}
	file << "H eta rho T K price iv stat\n";
	// set precision
	file << std::setprecision(10);
	for (size_t i = 0; i < res.H.size(); ++i)
		file << res.H[i] << " " << res.eta[i] << " " << res.rho[i] << " "
				<< res.T[i] << " " << res.K[i] << " " << res.price[i] << " "
				<< res.iv[i] << " " << res.stat[i] << "\n";
	file.close();

	std::cout << "Time elapsed: " << diff.count() << "ms.\n";
*/
	return 0;
}
