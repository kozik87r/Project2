#include <iostream>
#include "CPUTimer.h"
#include <typeinfo>

#include "ADEngine.h"

#include <mkl.h>
#include <mkl_vsl.h>

#ifndef MAX_PATHS
#define MAX_PATHS 100000
#endif

#ifndef VECTOR_SIZE
#define VECTOR_SIZE 1
#endif

using namespace std;

int main()
{
	int numberOfSteps = 100;
	string volString;
	string stString;
	DATA_TYPE optionPrice = 0;
	DATA_TYPE squareOptionPrice = 0;
	DATA_TYPE T = 1.0;
	DATA_TYPE strike = 50;

	for(int i=0;i<NUM_VARIABLES;i++)
	{
		values[i] = 0;
		sensitivities[i] = 0;
	}

	// St
	Variable<0> St;
	Variable<1> rate;
	Variable<2> q;
	Variable<3> half;
	Variable<4> Vtt;
	Variable<5> dt;
	Variable<6> Zs;

	// Vt
	Variable<7> Vt;
	Variable<8> kappa;
	Variable<9> theta;
	Variable<10> zero;
	Variable<11> sigma;
	Variable<12> Zv;

	St = 100;
	rate = 0.005;
	q = 0;
	half = 0.5;
	Vtt = 0.02;
	dt = T/numberOfSteps;
	Zs = 1;
	kappa = 0.4;
	theta = 0.1;
	Vt = 0.02;
	sigma = 0.3;
	zero = 0;

	double rho = 0.3;
	double xmu = 0.0;
	double var = 1.0;
	double* randZ1 = new double[VECTOR_SIZE];
	double* randZ2 = new double[VECTOR_SIZE];
	double* randZ3 = new double[VECTOR_SIZE];

	#define ADJOINT


	auto nextSt = St * exp(((rate - q) - half * Vtt) * dt + sqrt(Vtt * dt) * Zs);
	auto nextVt = Vt + kappa*(theta-max2(Vt, zero))*dt + sigma * sqrt(max2(Vt, zero)*dt) * Zv;

	stString = typeid(nextSt).name();
	volString = typeid(nextVt).name();

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

	do{
		found = volString.find("class");
		if(found > volString.size())
			break;
		volString.replace(found, 5, "");
	}while(found<volString.size() && found!=std::string::npos);

	do{
		found = volString.find(" ");
		if(found > volString.size())
			break;
		volString.replace(found, 1, "");
	}while(found<volString.size() && found!=std::string::npos);

	cout << stString.c_str() << endl << endl;
	cout << volString.c_str() << endl << endl;

	int stepSize = 10000;


	CPUTimer cpuTimer;
	double executionTime;

	for(int pathsNr = 10000; pathsNr <= MAX_PATHS; pathsNr += stepSize)
	{
		if(pathsNr>=50000)
			stepSize = 50000;
		cpuTimer.startTime();

		for(int pathId=0;pathId<pathsNr;pathId++)
		{
			VSLStreamStatePtr  stream0, stream1, stream2, stream3;
			vslNewStream(&stream0, VSL_BRNG_MCG59, pathId);
			vslCopyStream(&stream1, stream0);  // Create a second instance of the same randome stream
			vslCopyStream(&stream2, stream0);  // Create a second instance of the same randome stream
			vslCopyStream(&stream3, stream0);  // Create a second instance of the same randome stream
			vslLeapfrogStream( stream0, 0, 4); // The master stream will have 2 partitions.  0 takes
			vslLeapfrogStream( stream1, 1, 4); // the first, 1 the second
			vslLeapfrogStream( stream2, 2, 4); // the first, 1 the second
			vslLeapfrogStream( stream3, 3, 4); // the first, 1 the second

			St = 100;
			Vt = 0.02;
			Vtt = 0.02;
			for(int stepId = 0; stepId < numberOfSteps; stepId++)
			{
				vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream0, VECTOR_SIZE, randZ1, xmu, var);
				vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream1, VECTOR_SIZE, randZ2, xmu, var);

				Zs.value = randZ1[0];
				Zv.value = rho * Zs.value + sqrt(1-rho*rho)*randZ1[0];

				nextVt.evaluate();

#if defined(ADJOINT)
				nextVt.setDiff();
				nextVt.diff();
#endif

				nextSt.evaluate();

#if defined(ADJOINT)
				nextSt.setDiff();
				nextSt.diff();
#endif

				St.value = nextSt.value;
				Vtt.value = Vt.value = nextVt.value;

				vslSkipAheadStream(stream0, VECTOR_SIZE);
				vslSkipAheadStream(stream1, VECTOR_SIZE);
				vslSkipAheadStream(stream2, VECTOR_SIZE);
			}

			optionPrice +=  exp(-rate.value*T)*max(St.value-strike, 0);
			squareOptionPrice += optionPrice;
			vslDeleteStream(&stream0);
			vslDeleteStream(&stream1);
			vslDeleteStream(&stream2);
			vslDeleteStream(&stream3);
		}
		executionTime =  cpuTimer.calculateElapsedTime();
		cout << pathsNr << "\t" << executionTime << endl;
	}


	Mul<Variable<0>,Exp<Add<Mul<Sub<Sub<Variable<1>,Variable<2>>,Mul<Variable<3>,Variable<4>>>,\
		Variable<5>>,Mul<Sqrt<Mul<Variable<4>,Variable<5>>>,Variable<6>>>>> nextSt2;

	Add<Add<Variable<7>,Mul<Mul<Variable<8>,Sub<Variable<9>,Max<Variable<7>,Variable<10>>>>,\
		Variable<5>>>,Mul<Mul<Variable<11>,Sqrt<Mul<Max<Variable<7>,Variable<10>>,Variable<5>>>>,Variable<12>>> nextVt2;


	stepSize = 10000;
	optionPrice = 0;
	squareOptionPrice = 0;

	for(int pathsNr = 10000; pathsNr <= MAX_PATHS; pathsNr += stepSize)
	{
		if(pathsNr>=50000)
			stepSize = 50000;
		cpuTimer.startTime();


		for(int pathId=0;pathId<pathsNr;pathId++)
		{
			VSLStreamStatePtr  stream0, stream1, stream2, stream3;
			vslNewStream(&stream0, VSL_BRNG_MCG59, pathId);
			vslCopyStream(&stream1, stream0);  // Create a second instance of the same randome stream
			vslCopyStream(&stream2, stream0);  // Create a second instance of the same randome stream
			vslCopyStream(&stream3, stream0);  // Create a second instance of the same randome stream
			vslLeapfrogStream( stream0, 0, 4); // The master stream will have 2 partitions.  0 takes
			vslLeapfrogStream( stream1, 1, 4); // the first, 1 the second
			vslLeapfrogStream( stream2, 2, 4); // the first, 1 the second
			vslLeapfrogStream( stream3, 3, 4); // the first, 1 the second

			St = 100;
			Vt = 0.02;
			Vtt = 0.02;
			for(int stepId = 0; stepId < numberOfSteps; stepId++)
			{



				vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream0, VECTOR_SIZE, randZ1, xmu, var);
				vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream1, VECTOR_SIZE, randZ2, xmu, var);

				Zs.value = randZ1[0];
				Zv.value = rho * Zs.value + sqrt(1-rho*rho)*randZ1[0];

				nextVt2.evaluate();

#if defined(ADJOINT)
				nextVt2.setDiff();
				nextVt2.diff();
#endif
				nextSt2.evaluate();

#if defined(ADJOINT)
				nextSt2.setDiff();
				nextSt2.diff();
#endif

				St.value = nextSt2.value;
				Vtt.value = Vt.value = nextVt2.value;

				vslSkipAheadStream(stream0, VECTOR_SIZE);
				vslSkipAheadStream(stream1, VECTOR_SIZE);
				vslSkipAheadStream(stream2, VECTOR_SIZE);
			}
			optionPrice +=  exp(-rate.value*T)*max(St.value-strike, 0);
			squareOptionPrice += optionPrice;

			vslDeleteStream(&stream0);
			vslDeleteStream(&stream1);
			vslDeleteStream(&stream2);
			vslDeleteStream(&stream3);

		}

		executionTime =  cpuTimer.calculateElapsedTime();
		cout << pathsNr << "\t" << executionTime << endl;
	}
	return 0;
};




/*


for(int i=0;i<7;i++)
{
cout << "sensitivity[" << i << "]:=" << sensitivities[i] << endl;
}
// the template argument index refers to the index in the array of input variables


//class Add<class Add<class Mul<class Variable<0>,class Variable<1> >,class Variab
//	le<1> >,class Mul<class Mul<class Variable<2>,class Variable<1> >,class Variable
//	<0> > >


for(int i=0;i<NUM_VARIABLES;i++)
{
//	values[i] = 0;
sensitivities[i] = 0;
}

for(int k=0;k<1000;k++)
{
Mul<Variable<0>,Exp<Add<Mul<Sub<Sub<Variable<1>,Variable<2>>,Mul<Variable<3>, \
Variable<4>>>,Variable<5>>,Mul<Sqrt<Mul<Variable<4>,Variable<5>>>,Variable<6>>>>> nextSt2;

nextSt2.evaluate();

//	cout << "---" << endl;
//	cout << "Value: " << f5.evaluate() << endl;
nextSt2.setDiff();
nextSt2.diff();
}

for(int i=0;i<7;i++)
{
cout << "sensitivity[" << i << "]:=" << sensitivities[i] << endl;
}

//	cout << " " << f3.evaluate();
//	f3.diff();
return 0;
};*/
