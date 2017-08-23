/*
 * gaussContainer.h
 *
 *  Created on: Aug 9, 2017
 *      Author: bayerc
 */

#ifndef GAUSSCONTAINER_H_
#define GAUSSCONTAINER_H_

#include<vector>
#include<iostream>
#include<cstdlib>

typedef std::vector<double> Vector;

// class containing Gaussian random numbers to be successively returned
class gaussContainer{
private:
	std::vector<Vector> Z;
	long i; // what is the current location
public:
	gaussContainer(std::vector<Vector> X);
	gaussContainer();
	Vector operator()();
	size_t size() const {return Z.size();}
};

// divide a vector<Vector> into k gaussContainers
std::vector<gaussContainer> containerize(std::vector<Vector> X, size_t k);

void printGCSizes(const std::vector<gaussContainer>& x);

#endif
