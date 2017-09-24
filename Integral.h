#ifndef INTEGRAL_H
#define INTEGRAL_H

template <class Type>
class HestonFormula;

template <class Type>
class Integral
{
public:
	Integral();
	~Integral() {};
	Type GetIntegrationResult(HestonFormula<Type> *hestonFormula, Type &a, Type &b, Type &tolerance, Type* params, int type);
private:
	int iter_;
	int nmaxiter_;
	Type wg[4]; //weights of 7 point Gauss rule
	Type wgk[8]; //weights of 15 point Kronrod rule
	Type xgk[8]; //nodes of 15 point Kronrod rule
};

template <class Type>
Integral<Type>::Integral()
{
	iter_ = 0;
	nmaxiter_ = 10000;
	wg[0] = 0.417959183673469;
	wg[1] = 0.381830050505119;
	wg[2] = 0.279705391489277;
	wg[3] = 0.12948496616887;

	wgk[0] = 0.209482141084728;
	wgk[1] = 0.204432940075298;
	wgk[2] = 0.190350578064785;
	wgk[3] = 0.169004726639267;
	wgk[4] = 0.140653259715525;
	wgk[5] = 0.10479001032225;
	wgk[6] = 0.063092092629979;
	wgk[7] = 0.022935322010529;

	xgk[0] = 0;
	xgk[1] = 0.207784955007898;
	xgk[2] = 0.405845151377397;
	xgk[3] = 0.586087235467691;
	xgk[4] = 0.741531185599394;
	xgk[5] = 0.864864423359769;
	xgk[6] = 0.949107912342758;
	xgk[7] = 0.991455371120813;
};
static int depth = 0;


template <class Type>
Type Integral<Type>::GetIntegrationResult(HestonFormula<Type> *hestonFormula, Type &a, Type &b, Type &tolerance, Type* params, int type)
{

	Type half ;//half length of the Type
	Type two;
	two = 2.0;
	half = (b - a) / two;
	Type center;
	center = (a + b) / two;
	Type gsum, k15, t, fsum, fc;
	fc = hestonFormula->HestonPIntegrand(center, params, type);
	//return res;


	depth++;
	gsum = fc * wg[0];
	k15 = fc * wgk[0];
	int j,j2;
	j2 = 2;

	for (j=1;j<=3;j++)
	{
		t = half * xgk[j2];
		fsum = hestonFormula->HestonPIntegrand(center - t,params, type) + hestonFormula->HestonPIntegrand(center + t,params, type);
		gsum = gsum + fsum * wg[j];
		k15 = k15 + fsum * wgk[j2];
		j2 = j2 + 2;
	}


	for (j2 = 1;j2<=7;j2=j2+2)
	{
		t = half * xgk[j2];
		fsum = hestonFormula->HestonPIntegrand(center - t,params, type) + hestonFormula->HestonPIntegrand(center + t,params, type);
		k15 = k15 + fsum * wgk[j2];
	}

	gsum = half * gsum;
	k15 = half * k15;
	iter_ = iter_ + 15;

	//cout << "fsum: " << currentTape[fsum.idx].value << endl;
	//cout << "currentTape[k15.idx].value " << currentTape[k15.idx].value << endl;
	//cout << "gsum " << currentTape[gsum.idx].value << endl;
	//cout << "tolerancevalue " << currentTape[tolerance.idx].value << endl;
	Type res;
	//cout << "fabs: " << fabs(k15 - gsum) << endl;
//	cout << "tolerance: " << tolerance << endl;
	if (fabs(k15 - gsum) < tolerance) 
	{
		res = k15;
		return res;
	}
	else
	{
		if ( (iter_ + 30) > nmaxiter_ ) 
		{
			std::cout << "maximum number of function evaluations exceeded";
			return 0;
		}
		
		res = GetIntegrationResult(hestonFormula, a, center, tolerance, params, type) + GetIntegrationResult(hestonFormula, center, b, tolerance, params, type);
		depth--;
	}return res;

};


#endif
