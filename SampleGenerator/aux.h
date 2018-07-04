// Auxiliary functions

#pragma once

#include<vector>

typedef std::vector<double> Vector;

template<typename T>
inline void scaleVector(std::vector<T>& x, T s){
	for(auto& v : x)
		v = v*s;
}

// z = a*x+b*y
template<typename T>
inline std::vector<T> linearComb(T a, const std::vector<T>& x, T b, const std::vector<T>& y){
	std::vector<T> z(x.size());
	for(size_t i=0; i<z.size(); ++i)
		z[i] = a*x[i] + b*y[i];
	return z;
}
