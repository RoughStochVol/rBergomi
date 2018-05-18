/*
 * aux.h
 *
 *  Created on: 16 May 2018
 *      Author: bayerc
 */

#include <cmath>


// integer logarithm for base 2
inline int ilog2(int n){
	int ret = 0;
	if(n > 0){
		while (n >>= 1)
			++ret;
	} else
		ret = -1;
	return ret;
}

// double logarithm for base 2
inline double log2(double x){
	return log(x) / log(2);
}
