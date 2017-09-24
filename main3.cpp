#include <iostream>
#include "CPUTimer.h"
#include <typeinfo>
#include <cxxabi.h>

#include "ADEngine.h"
#include "HestonAnalyticFormula.h"

#ifndef MAX_PATHS
#define MAX_PATHS 100000
#endif

#ifndef VECTOR_SIZE
#define VECTOR_SIZE 1
#endif

using namespace std;

int main()
{
	string stString;

	for(int i=0;i<NUM_VARIABLES;i++)
	{
		values[i] = 0;
		sensitivities[i] = 0;
	}

	//HestonFormula hestonFormula;


	// St
	Variable<__COUNTER__> kappa;
	Variable<__COUNTER__> theta;
	Variable<__COUNTER__> sigma;
	Variable<__COUNTER__> corr;
	Variable<__COUNTER__> Vt;
	Variable<__COUNTER__> rate;
	Variable<__COUNTER__> T;
	Variable<__COUNTER__> St;
	Variable<__COUNTER__> strike;

	corr = 0;
	St = 100;
	rate = 0.005;
	kappa = 0.4;
	theta = 0.1;
	Vt = 0.02;
	sigma = 0.3;
	strike = 50;
	T = 1;

	//HestonFormula<double> hs;
	//double res;
	ForLoop<200> f;
	auto res2 = kappa*theta;
	auto res = f.add(corr, res2);
	//HestonCallQuad(kappa, theta, sigma, corr, Vt, rate, T, St, strike);// func(kappa);
	//hs.HestonCallQuad(kappa, theta, sigma, corr, Vt, rate, T, St, strike);
	//auto f = HestonCallQuad(kappa, theta, sigma, corr, Vt, rate, T, St, strike);


#define ADJOINT
	int status;
	stString = abi::__cxa_demangle(typeid(res).name(), 0, 0, &status);

	std::size_t found;
	do{
		found = stString.find("class");
		if(found > stString.size())
			break;
		stString.replace(found, 5, "");
	}while(found<stString.size() && found!=std::string::npos);

	do{
		found = stString.find(" ");
		if(found > stString.size())
			break;
		stString.replace(found, 1, "");
	}while(found<stString.size() && found!=std::string::npos);


	cout << stString.c_str() << endl << endl;


	CPUTimer cpuTimer;
	double executionTime;
	return 0;
};
