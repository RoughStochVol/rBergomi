// Auxiliary functions

#pragma once

#include<iostream>
#include<fstream>
#include<vector>
#include<iterator>
#include<string>
#include<sstream>
#include<cstdlib>
#include<cmath>
#include<algorithm>
#include<numeric>
#include<fftw3.h>
#include"ParamTot.h"

template <typename T>
inline std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  if ( !v.empty() ) {
    out << '(';
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b)";
  }
  return out;
}

inline std::ostream& operator<< (std::ostream & out, const ParamTot& par) {
    out << "\n# H eta rho T K";
    for(long i=0; i<par.size(); ++i)
    	out << "\n" << i << " " << par.H(i) << " " << par.eta(i)
			<< " " << par.rho(i) << " " << par.T(i) << " "
			<< par.K(i);
    out << std::endl;
    return out;
}

template<typename T>
inline T posPart(T x){
	return x > 0 ? x : 0;
}

template<typename T>
inline void scaleVector(std::vector<T>& x, T s){
	for(auto& v : x)
		v = v*s;
}

// x = x + y
template<typename T>
inline void addTo(std::vector<T>& x, const std::vector<T>& y){
	for(size_t i = 0; i<x.size(); ++i)
		x[i] = x[i] + y[i];
}

// z = a*x+b*y
template<typename T>
inline std::vector<T> linearComb(T a, const std::vector<T>& x, T b, const std::vector<T>& y){
	std::vector<T> z(x.size());
	for(size_t i=0; i<z.size(); ++i)
		z[i] = a*x[i] + b*y[i];
	return z;
}

template<typename T>
inline std::vector<T> squareVector(const std::vector<T>& x){
	std::vector<T> y = x;
	for(auto& v : y)
		v = v*v;
	return y;
}

template<typename T>
inline std::vector<T> rootVector(const std::vector<T>& x){
	std::vector<T> y = x;
	for(auto& v : y)
		v = sqrt(v);
	return y;
}

// read file and save in vector of given type
template<typename T>
inline std::vector<T> file2vector(std::string fname){
	std::ifstream file(fname.c_str(), std::ios::in);
	if(!file.is_open()){
		std::cerr << "\nError while opening file " << fname << ".\n";
		exit(27);
	}
	T num;
	std::vector<T> vec;
	while(file >> num)
		vec.push_back(num);

	file.close();
	return vec;
}

// Hierarchical representation of Gaussian random variables.
// This function transforms a IID standard normal Vector Z
// (of dimension M = 2^m for simplicity (for some m)) into another
// IID standard normal vector dW (of the same size). If we construct
// a Brownian motion W by letting W(t_{i+1}) - W(t_i) = \sqrt{t_{i+1} - t_i} dW[i],
// then this corresponds to the Brownian bridge construction of
// W in terms of the original random vector Z. I.e., we have
// W(t_M) = \sqrt{t_M} Z[0], ...
inline void hierarchical2increments(const Vector& Z, Vector& dW){
	// read out m from the input Z
	int m = log2(Z.size());
	// Check that size of Z is = 2^m
	if(Z.size() != exp2(m)){
		std::cerr << "Error in hierarchical2increments.\nSize of input vector Z is not an exact power of 2.\n";
		exit(15);
	}
	for(int i=0; i<exp2(m); ++i){
		dW[i] = exp2(-0.5*m) * Z[0];
		for(int n=1; n<=m; ++n){
			int i2nm = floor(i * exp2(n - m));
			int iota = (i2nm + 1) % 2; // = 1 if i2nm is even and 0 else
			//int i1 = exp2(n-1) + floor(i * exp2(n-m));
			//int i2 = i1 - 1;
			int weight = 2 * iota - 1; // 1 if i2nm is even and -1 else
			int index = exp2(n-1) + (i2nm + iota + 1) / 2 - 1; // this has to be an integer! Is that a dangerous assumption?
			dW[i] += exp2(0.5*(n-m-1)) * weight * Z[index];
		}
	}
}

inline std::vector<std::string> hierarchical2string(int m){
	std::vector<std::string> dW(exp2(m), "");
	std::stringstream ss;
	for(int i=0; i<exp2(m); ++i){
		ss.str("");
		ss << exp2(-0.5*m) << " Z[0]";
		for(int n=1; n<=m; ++n){
			int i2nm = floor(i * exp2(n - m));
			int iota = (i2nm + 1) % 2; // = 1 if i2nm is even and 0 else
			//int i1 = exp2(n-1) + floor(i * exp2(n-m));
			//int i2 = i1 - 1;
			int weight = 2 * iota - 1; // 1 if i2nm is even and -1 else
			int index = exp2(n-1) + (i2nm + iota + 1) / 2 - 1; // "-1" due to C style array indexing
			ss << " + " << exp2(0.5*(n-m-1)) * weight << " Z[" << index << "]";
		}
		dW[i] = ss.str();
	}
	return dW;
}

// For debugging purposes: compute the "expected" variances of increments given
inline Vector hierarchical2incrementsVar(int m){
	Vector var(exp2(m), 0.0);
	for(int i=0; i<exp2(m); ++i){
		var[i] = pow(exp2(-0.5*m), 2);
		for(int n=1; n<=m; ++n){
			int i2nm = floor(i * exp2(n - m));
			int iota = (i2nm + 1) % 2; // = 1 if i2nm is even and 0 else
			//int i1 = exp2(n-1) + floor(i * exp2(n-m));
			//int i2 = i1 - 1;
			int weight = 2 * iota - 1; // 1 if i2nm is even and -1 else
			int index = exp2(n-1) + (i2nm + iota + 1) / 2 - 1; // this has to be an integer! Is that a dangerous assumption?
			var[i] += pow(exp2(0.5*(n-m-1)) * weight, 2);
		}
	}
	return var;
}

// Compute sample variances
inline Vector sampleVar(std::vector<Vector> Z){
	size_t M = Z.size();
	size_t N = Z[0].size();
	Vector sample_mean(N, 0.0);
	Vector sample_var(N, 0.0);
	for(size_t i=0; i<M; ++i){
		for(size_t j = 0; j<N; ++j){
			sample_mean[j] += Z[i][j];
			sample_var[j] += Z[i][j] * Z[i][j];
		}
	}
	// compute sample mean and sample var.
	scaleVector(sample_mean, 1.0/M);
	scaleVector(sample_var, 1.0/M);
	sample_var = linearComb(1.0, sample_var, -1.0, squareVector(sample_mean));
	return sample_var;
}

// Put even indices of z into x1 and odd ones into x2
inline void breakZ(const Vector& z, Vector& x1, Vector& x2){
	for(size_t i = 0; i < x1.size(); ++i){
		x1[i] = z[2*i];
		x2[i] = z[2*i+1];
	}
}

// return the minimal value in a vector
template <typename T>
inline T minVec(const std::vector<T>& x){
	auto iter = std::min_element(x.begin(), x.end());
	return *iter;
}

// A series of "small" deterministic vectors
inline void debugFillVector(Vector& x, int m){
	for(size_t i=0; i<x.size(); ++i)
		x[i] = cos(double(i)*double(m));
}

// Compute covariance between two samples
inline double covariance(const Vector& x, const Vector& y){
	if(x.size() != y.size()){
		std::cerr << "Computing covariances for samples of different size.\n";
		exit(34);
	}
	double xMean = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
	double yMean = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
	Vector xDiff(x.size()); // x - xMean
	Vector yDiff(y.size());
	std::transform(x.begin(), x.end(), xDiff.begin(), [xMean](double x) { return x - xMean; });
	std::transform(y.begin(), y.end(), yDiff.begin(), [yMean](double x) { return x - yMean; });
	double covariance = std::inner_product(xDiff.begin(), xDiff.end(), yDiff.begin(), 0.0) / (x.size() - 1);
	return covariance;
}
