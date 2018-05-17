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
// out.txt ... name of the output-text-file
// in_name ... root name for input texts: for instance, H is read from file .in_text_H.txt
// method ... 's' for standard MC, 't' for "turbo-charged" MC
// Further input variables:
// H ... read from local file .in_text_H.txt
// eta ... read from local file .in_text_eta.txt
// rho ... read from local file .in_text_rho.txt
// xi ... read from local file .in_text_xi.txt
// T ... read from local file .in_text_T.txt
// K ... read from local file .in_text_K.txt
// Those files are supposed to only contain the corresponding numbers.

#include <chrono>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <random>
#include <fenv.h>
#include "RBergomi.h"
//#include "BlackScholes.h"
//#include "aux.h"

typedef std::vector<double> Vector;

// There are other clocks, but this is usually the one you want.
// It corresponds to CLOCK_MONOTONIC at the syscall level.
using Clock = std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

int main(int argc, char* argv[]) {
	feenableexcept(FE_INVALID);// | FE_OVERFLOW);

	if ((argc != 7) && (argc != 1)) {
		std::cerr
				<< "Error: Wrong number of arguments.\nUsage: rBergomi N M path out_file in_name method\n";
		exit(5);
	}

	long N, M;
	Vector H, eta, rho, T, K, xi;
	std::string path_name, in_name, out_name;
	char method;

	if (argc == 1) {
		N = 100;
		M = 100000;
		H = Vector { 0.05, 0.2 };
		eta = Vector { 1.0, 3.0 };
		rho = Vector { -0.98, -0.8 };
		T = Vector { 0.05, 2.0 };
		K = Vector { 1.0, 1.3 };
		xi = Vector { 0.04, 0.04 };
		method = 't';
	} else {
		// Convert the console input
		N = atoi(argv[1]);
		M = atoi(argv[2]);
		path_name = std::string(argv[3]);
		out_name = std::string(argv[4]);
		in_name = std::string(argv[5]);
		method = *(argv[6]);

		// Double check: output of console arguments as recorded
		std::cerr << "N = " << N << ", M = " << M << ",\npath = " << path_name << ",\noutfile = " << out_name
				<< "\nmethod = " << method
				<< std::endl;

		// read in the further input files.
		const std::string hName = path_name + std::string(".") + in_name
				+ std::string("H.txt");
		H = file2vector<double>(hName);
		const std::string etaName = path_name + std::string(".") + in_name
				+ std::string("eta.txt");
		eta = file2vector<double>(etaName);
		const std::string rhoName = path_name + std::string(".") + in_name
				+ std::string("rho.txt");
		rho = file2vector<double>(rhoName);
		const std::string TName = path_name + std::string(".") + in_name
				+ std::string("T.txt");
		T = file2vector<double>(TName);
		const std::string KName = path_name + std::string(".") + in_name
				+ std::string("K.txt");
		K = file2vector<double>(KName);
		const std::string xiName = path_name + std::string(".") + in_name
				+ std::string("xi.txt");
		xi = file2vector<double>(xiName);
	}

	// Double check that all the inputs have positive size
	if ((H.size() <= 0) || (eta.size() <= 0) || (rho.size() <= 0)
			|| (T.size() <= 0) || (K.size() <= 0) || (xi.size() <= 0)) {
		std::cerr << "One or more parameter vectors have size 0.\n";
		exit(17);
	}
	// Double check that all inputs have the same length
	if ((eta.size() != H.size()) || (rho.size() != H.size())
			|| (T.size() != H.size()) || (K.size() != H.size())
			|| (xi.size() != H.size())) {
		std::cerr << "The parameter arrays are not equal in size.\n";
		exit(18);
	}

	// generate seed;
	std::vector<uint64_t> seed(1);
	std::random_device rd;
	for (size_t i = 0; i < seed.size(); ++i)
		seed[i] = rd();

	// bugfix: output the seed
	std::cout << "Seed = " << seed[0] << ", K[0] = " << K[0] << std::endl;

	// run the code
	auto start = Clock::now();
	RBergomi rberg(xi, H, eta, rho, T, K, N, M, seed);
	Result res;
	if (method == 's')
		res = rberg.ComputeIVRT();
	else if (method == 't')
		res = rberg.ComputeIVRTVarRed();
	// bugfix
	std::cerr << "Computation completed." << std::endl;
	auto end = Clock::now();
	auto diff = duration_cast<milliseconds>(end - start);

	// output results
	if (argc > 1) {
		// print result in output file
		std::string out_name_full = path_name + out_name;
		std::ofstream file;
		file.open(out_name_full.c_str(),
				std::ofstream::out | std::ofstream::trunc);
		// bugfix
		std::cerr << "Writing results in file " << out_name_full << std::endl;
		if (!file.is_open()) {
			std::cerr << "/nError while opening file " << out_name << "./n";
			exit(1);
		}
		file << "xi H eta rho T K price iv stat\n";
		// set precision
		file << std::setprecision(10);
		for (size_t i = 0; i < H.size(); ++i)
			file << res.par.xi(i) << " " << res.par.H(i) << " "
					<< res.par.eta(i) << " " << res.par.rho(i) << " "
					<< res.par.T(i) << " " << res.par.K(i) << " "
					<< res.price[i] << " " << res.iv[i] << " " << res.stat[i]
					<< "\n";
		file.close();
	} else {
		for (size_t i = 0; i < H.size(); ++i)
			std::cout << res.par.xi(i) << " " << res.par.H(i) << " "
					<< res.par.eta(i) << " " << res.par.rho(i) << " "
					<< res.par.T(i) << " " << res.par.K(i) << " "
					<< res.price[i] << " " << res.iv[i] << " " << res.stat[i]
					<< "\n";
	}

	std::cout << "Time elapsed: " << diff.count() << "ms.\n";

	return 0;
}
