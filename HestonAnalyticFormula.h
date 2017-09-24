#ifndef HESTON_FORMULA_H
#define HESTON_FORMULA_H

#include "Integral.h"
#include "Complex.h"

#include <math.h>
// ok
#define M_PI       3.14159265358979323846

template <class Value>
class HestonFormula{
public:


	HestonFormula() {};

	~HestonFormula() {};

	ComplexType<Value> Hestf(ComplexType<Value> phi, ComplexType<Value> kappa, ComplexType<Value> theta, ComplexType<Value> sigma,
		ComplexType<Value> rho,ComplexType<Value> v0, ComplexType<Value> r,ComplexType<Value> T, ComplexType<Value> S0, int type)
	{
		//cout << "phi " << phi << endl;
	//	cout << "kappa " << kappa << endl;
	//	cout << "theta " << theta << endl;
		//cout << "sigma " << sigma << endl;
		///cout << "rho " << rho << endl;
		//cout << "v0 " << v0 << endl;

		Value one;
		Value two;
		Value zero;
		Value half;

		one = 1.0;
		two = 2.0;
		zero = 0;
		half = 0.5;

		ComplexType<Value> u;
		ComplexType<Value> b;
		ComplexType<Value> res;
		ComplexType<Value> c1 = ComplexType<Value>(one, zero);
		ComplexType<Value> c2 = ComplexType<Value>(two, zero);
		ComplexType<Value> sig2 = sigma*sigma;

		if (type == 1)
		{
			u = ComplexType<Value>(half, zero);
			b = kappa - rho*sigma;
		}
		else
		{
			u = ComplexType<Value>(-half, zero);
			b = kappa;
		}

		ComplexType<Value> a = kappa*theta;
		ComplexType<Value> x;
		x = log(S0);
		ComplexType<Value> tmp1;

		tmp1 = ComplexType<Value>(zero, (rho*sigma*phi).real());

		ComplexType<Value> tmp2;
		tmp2 = ComplexType<Value>(zero, (c2*u*phi).real());

		ComplexType<Value> d = pow((tmp1-b)*(tmp1-b)-sig2*(tmp2-phi*phi), half);

		ComplexType<Value> g = (b-tmp1+d)/(b-tmp1-d);
		tmp2 = ComplexType<Value>(zero, (r*phi*T).real());

		ComplexType<Value> C = tmp2 + (a/sig2) * ((b- tmp1 + d)*T - c2*log((c1-g*exp(d*T)) / (c1-g)));

		ComplexType<Value> D = (b-tmp1+d) / sig2 * ((c1 - exp<Value>(d*T))/ (c1 - g*exp<Value>(d*T)));

		tmp2 = ComplexType<Value>(zero, (phi*x).real());

		ComplexType<Value> f = exp<Value>(C + D*v0 + tmp2);

		return f;
	};

	Value HestonPIntegrand(Value phi, Value* params, int type)
	{
		Value res;
		res = 0;
		Value kappa = params[0];
		Value theta = params[1];
		Value sigma = params[2];
		Value rho = params[3];
		Value v0 = params[4];
		Value r = params[5];
		Value T = params[6];
		Value S0 = params[7];
		Value K = params[8];

		Value zero;
		zero = 0;
		ComplexType<Value> tmp1;
		tmp1 = ComplexType<Value>(zero, -phi*log(K));
		//cout << "tmp1: " << tmp1.real() << endl;
		ComplexType<Value> tmp2;
		tmp2 = Hestf(
			ComplexType<Value>(phi, zero),
			ComplexType<Value>(kappa, zero),
			ComplexType<Value>(theta, zero),
			ComplexType<Value>(sigma, zero),
			ComplexType<Value>(rho, zero),
			ComplexType<Value>(v0, zero),
			ComplexType<Value>(r, zero),
			ComplexType<Value>(T, zero),
			ComplexType<Value>(S0, zero),
			type);
//		cout << "tmp2: " << currentTape[tmp2.real().idx].value << ", " << currentTape[tmp2.imag().idx].value << endl;
		ComplexType<Value> tmp3 = ComplexType<Value>(zero, phi);

		res = (exp(tmp1)*tmp2/tmp3).real();
		//cout << "res: " << currentTape[res.idx].value << ", " << currentTape[res.idx].value << endl;
		//cout << counter << " res: " << res << endl;
		return res;
	}

	Value HestonP(int type, Value kappa, Value theta, Value sigma,
		Value rho, Value v0, Value r, Value T, Value S0, Value K)
	{
		Integral<Value> integral;
		Value res;
		res = 0;

		Value params[10];
		params[0] = kappa;
		params[1] = theta;
		params[2] = sigma;
		params[3] = rho;
		params[4] = v0;
		params[5] = r;
		params[6] = T;
		params[7] = S0;
		params[8] = K;
		params[9] = type;
		params[0] = kappa;

		Value error;
		error = 0.001;
		Value a;
		Value b;
		a = 0.01;
		b = 100;
		Value area = integral.GetIntegrationResult(this, a, b, error, params, type);

		Value pi = M_PI;

		res = 0.5 + (1 / pi) * area;

		//cout << "res: " << res << endl;
		//	delete params;
		return res;
	}
};

template <class T1, class T2>
class ComplexType2
{
public:
	ComplexType2(T1 a, T2 b) {};
	~ComplexType2() {};
};


template <class F1,
	class F2,
	class F3,
	class F4,
	class F5,
	class F6,
	class F7,
	class F8,
	class F9>
auto HestonCallQuad(F1 &kappa, F2 &theta, F3 &sigma,
			F4 &rho, F5 &v0, F6 &r, F7 &T, F8 &S0, F9 &K)
	{

//	auto c1 = ComplexType2<decltype(res2), decltype(T)>(res2, T);
		//auto res = K*K*kappa+theta*sigma+S0;
		return T*K;
	}

template <int N>
class ForLoop
{
public:
	ForLoop() {};
	~ForLoop() {};

	template <class Type1, class Type2>
	auto add(Type1 a, Type2 b ) { return (a + obj.add(a, b)); };

	template <class Type>
	auto sub(Type type) { return (type - obj.add(type)); };


	template <class Type>
	auto multiply(Type type) { return (type * obj.add(type)); };


	template <class Type>
	auto divide(Type type) { return (type / obj.add(type)); };

	ForLoop<N-1> obj;
};

template <>
class ForLoop<0>
{
public:
	ForLoop() {};
	~ForLoop() {};
	template <class Type1, class Type2>
	auto add(Type1 a, Type2 b) { return b; };

	template <class Type>
	auto sub(Type type) { return type; };

	template <class Type>
	auto multiply(Type type) { return type; };

	template <class Type>
	auto divide(Type type) { return type; };
};


#endif
