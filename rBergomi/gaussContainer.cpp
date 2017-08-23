/*
 * gaussContainer.cpp
 *
 *  Created on: Aug 9, 2017
 *      Author: bayerc
 */

#include"gaussContainer.h"

Vector gaussContainer::operator()()
{
	if(i >= Z.size()){
		std::cerr << "The number of provided random numbers was too small.\n";
		std::cerr << "i = " << i << "Z.size = " << Z.size() << "\n";
		std::abort();
	}
	i = i + 1;
	return Z[i-1];

}

gaussContainer::gaussContainer(std::vector<Vector> X)
{
	Z = X;
	i = 0;
}

gaussContainer::gaussContainer()
{
	Z = std::vector<Vector>(0);
	i = 0;
}

std::vector<gaussContainer> containerize(std::vector<Vector> X, size_t k){
	size_t l = X.size() / k;
	std::vector<gaussContainer> cont(k);
	for(size_t i  = 0; i < cont.size()-1; ++i)
		cont[i] = gaussContainer(std::vector<Vector>(&(X[i*l]), &(X[i*l]) + l));
	cont[cont.size()-1] = gaussContainer(std::vector<Vector>(&(X[(cont.size() - 1)*l]), &(X[X.size()-1]) + 1));
	return cont;
}

void printGCSizes(const std::vector<gaussContainer>& x){
	std::cout << "Size of containers = ";
	for(size_t i=0; i<x.size(); ++i)
		std::cout << x[i].size() << " ";
	std::cout << "\n";
}
