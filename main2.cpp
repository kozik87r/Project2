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
	DATA_TYPE St;
	DATA_TYPE rate;
	DATA_TYPE q;
	DATA_TYPE half;
	DATA_TYPE Vtt;
	DATA_TYPE dt;
	DATA_TYPE Zs;

	// Vt
	DATA_TYPE Vt;
	DATA_TYPE kappa;
	DATA_TYPE theta;
	DATA_TYPE zero;
	DATA_TYPE sigma;
	DATA_TYPE Zv;

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


	DATA_TYPE nextSt; //= St * exp(((rate - q) - half * Vtt) * dt + sqrt(Vtt * dt) * Zs);
	DATA_TYPE nextVt; //= Vt + kappa*(theta-max2(Vt, zero))*dt + sigma * sqrt(max2(Vt, zero)*dt) * Zv;

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
	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10;
	double df10, df9, dSt, df8, df4, df7, df6, df5, dVtt, df3, ddt, df1, df2, drate, dq;
	double v1, v2, v3, v4, v5, v6, v7, v8;
	double dv8, dv7, dv6, dVt, dv4, dsigma, dv5, dv1, dddt, dv3, dkappa, dv2, dtheta;

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

				Zs = randZ1[0];
				Zv = rho * Zs + sqrt(1-rho*rho)*randZ1[0];

#ifndef ADJOINT
				nextSt = St * exp(((rate - q) - half * Vtt) * dt + sqrt(Vtt * dt) * Zs);
				nextVt = Vt + kappa*(theta-max(Vt, zero))*dt + sigma * sqrt(max(Vt, zero)*dt) * Zv;

#else

				
				f1 = f2 = f3 = f4 = f5 = f6 = f7 = f8 = f9 = f10 = 0;

				f1 = rate -q;
				f2 = half * Vtt;
				f3 = f1 + f2;
				f4 = f3*dt;
				f5 = Vtt * dt;
				f6 = sqrt(f5);
				f7 = f6*Zs;
				f8 = f4+f7;
				f9 = exp(f8);
				f10 = St * f9;
				nextSt = f10;
				
				df10 = df9 = dSt = df8 = df4 = df7 = df6 = df5 = dVtt = dt = df3 = ddt = df1 = df2 = drate = dq = 0;

				df10 += 1;
				df9 += St*df10;
				dSt += St*df10;
				df8 += exp(f8)*df9;
				df4 += df8;
				df7 += df8;
				df6 += Zs*df7;
				df5 += (0.5/sqrt(f5))*df6;
				dVtt += dt * df5;
				dt += Vtt * df5;
				df3 += dt * df4;
				ddt += f3*df4;
				df1 += df3;
				df2 += df3;
				dVtt += half * df2;
				drate += df1;
				dq = -df1;

			
				
				v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = 0;

				v1 = max(Vt, zero);
				v2 = theta - v1;
				v3 = kappa * v2;
				v4 = v3 * dt;
				v5 = v1*dt;
				v6 = sigma * v5 * Zv;
				v7 = Vt + v4;
				v8 = v7 + v6;
				nextVt = v8;
				
				dv8 = dv7 = dv6 = dVt = dv4 = dsigma = dv5 = dv1 = dddt = dv3 = dkappa = dv2 = dtheta = 0;
				dv8 += 1;
				dv7 += 1;
				dv6 += 1;
				dVt += 1;
				dv4 += 1;
				dsigma += v5*Zv*dv6;
				dv5 += sigma *Zv * dv6;
				dv1 += dt*dv5;
				dddt += v1*dv5;
				dv3 += dt*dv4;
				dddt += v3*dv4;
				dkappa += v2*dv3;
				dv2 += kappa*dv3;
				dtheta = dv2;
				dv1 += -dv2;




#endif

				//	nextSt.evaluate();


				St = nextSt;
				Vtt = Vt = nextVt;

				vslSkipAheadStream(stream0, VECTOR_SIZE);
				vslSkipAheadStream(stream1, VECTOR_SIZE);
				vslSkipAheadStream(stream2, VECTOR_SIZE);
			}

			optionPrice +=  exp(-rate*T)*max(St-strike, 0);
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

				Zs = randZ1[0];
				Zv = rho * Zs + sqrt(1-rho*rho)*randZ1[0];

				nextSt = St * exp(((rate - q) - half * Vtt) * dt + sqrt(Vtt * dt) * Zs);
				nextVt = Vt + kappa*(theta-max(Vt, zero))*dt + sigma * sqrt(max(Vt, zero)*dt) * Zv;
				//	nextSt.evaluate();


				St = nextSt;
				Vtt = Vt = nextVt;

				vslSkipAheadStream(stream0, VECTOR_SIZE);
				vslSkipAheadStream(stream1, VECTOR_SIZE);
				vslSkipAheadStream(stream2, VECTOR_SIZE);
			}
			optionPrice +=  exp(-rate*T)*max(St-strike, 0);
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
