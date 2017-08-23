// Auxiliary functions

#pragma once

#include<iostream>
#include<vector>
#include<iterator>
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
