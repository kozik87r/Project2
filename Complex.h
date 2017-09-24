// An implementation of C++ std::ComplexType for use on CUDA devices.
// Written by John C. Travers <jtravs@gmail.com> (2012)
//
// Missing:
//  - long double support (not supported on CUDA)
//  - some integral pow functions (due to lack of C++11 support on CUDA)
//
// Heavily derived from the LLVM libcpp project (svn revision 147853).
// Based on libcxx/include/ComplexType.
// The git history contains the complete change history from the original.
// The modifications are licensed as per the original LLVM license below.
//
// -*- C++ -*-
//===--------------------------- ComplexType ----------------------------------===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is dual licensed under the MIT and the University of Illinois Open
// Source Licenses. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#ifndef COMPLEX_H
#define COMPLEX_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#define M_INFINITY 0x7f800000
//#ifndef copysign
//#define copysign _copysign
//#endif

#define NAN  0x7f80ffff
#include <cmath>
#include <limits>
#define isnan(x) ((x) != (x))
#endif

/*
    ComplexType synopsis

template<class T>
class ComplexType
{
public:
    typedef T value_type;

    ComplexType(const T& re = T(), const T& im = T());
    ComplexType(const ComplexType&);
    template<class X> ComplexType(const ComplexType<X>&);

    T real() const;
    T imag() const;

    void real(T);
    void imag(T);

    ComplexType<T>& operator= (const T&);
    ComplexType<T>& operator+=(const T&);
    ComplexType<T>& operator-=(const T&);
    ComplexType<T>& operator*=(const T&);
    ComplexType<T>& operator/=(const T&);

    ComplexType& operator=(const ComplexType&);
    template<class X> ComplexType<T>& operator= (const ComplexType<X>&);
    template<class X> ComplexType<T>& operator+=(const ComplexType<X>&);
    template<class X> ComplexType<T>& operator-=(const ComplexType<X>&);
    template<class X> ComplexType<T>& operator*=(const ComplexType<X>&);
    template<class X> ComplexType<T>& operator/=(const ComplexType<X>&);
};

template<>
class ComplexType<float>
{
public:
    typedef float value_type;

    constexpr ComplexType(float re = 0.0f, float im = 0.0f);
    explicit constexpr ComplexType(const ComplexType<double>&);

    constexpr float real() const;
    void real(float);
    constexpr float imag() const;
    void imag(float);

    ComplexType<float>& operator= (float);
    ComplexType<float>& operator+=(float);
    ComplexType<float>& operator-=(float);
    ComplexType<float>& operator*=(float);
    ComplexType<float>& operator/=(float);

    ComplexType<float>& operator=(const ComplexType<float>&);
    template<class X> ComplexType<float>& operator= (const ComplexType<X>&);
    template<class X> ComplexType<float>& operator+=(const ComplexType<X>&);
    template<class X> ComplexType<float>& operator-=(const ComplexType<X>&);
    template<class X> ComplexType<float>& operator*=(const ComplexType<X>&);
    template<class X> ComplexType<float>& operator/=(const ComplexType<X>&);
};

template<>
class ComplexType<double>
{
public:
    typedef double value_type;

    constexpr ComplexType(double re = 0.0, double im = 0.0);
    constexpr ComplexType(const ComplexType<float>&);

    constexpr double real() const;
    void real(double);
    constexpr double imag() const;
    void imag(double);

    ComplexType<double>& operator= (double);
    ComplexType<double>& operator+=(double);
    ComplexType<double>& operator-=(double);
    ComplexType<double>& operator*=(double);
    ComplexType<double>& operator/=(double);
    ComplexType<double>& operator=(const ComplexType<double>&);

    template<class X> ComplexType<double>& operator= (const ComplexType<X>&);
    template<class X> ComplexType<double>& operator+=(const ComplexType<X>&);
    template<class X> ComplexType<double>& operator-=(const ComplexType<X>&);
    template<class X> ComplexType<double>& operator*=(const ComplexType<X>&);
    template<class X> ComplexType<double>& operator/=(const ComplexType<X>&);
};

// 26.3.6 operators:
template<class T> ComplexType<T> operator+(const ComplexType<T>&, const ComplexType<T>&);
template<class T> ComplexType<T> operator+(const ComplexType<T>&, const T&);
template<class T> ComplexType<T> operator+(const T&, const ComplexType<T>&);
template<class T> ComplexType<T> operator-(const ComplexType<T>&, const ComplexType<T>&);
template<class T> ComplexType<T> operator-(const ComplexType<T>&, const T&);
template<class T> ComplexType<T> operator-(const T&, const ComplexType<T>&);
template<class T> ComplexType<T> operator*(const ComplexType<T>&, const ComplexType<T>&);
template<class T> ComplexType<T> operator*(const ComplexType<T>&, const T&);
template<class T> ComplexType<T> operator*(const T&, const ComplexType<T>&);
template<class T> ComplexType<T> operator/(const ComplexType<T>&, const ComplexType<T>&);
template<class T> ComplexType<T> operator/(const ComplexType<T>&, const T&);
template<class T> ComplexType<T> operator/(const T&, const ComplexType<T>&);
template<class T> ComplexType<T> operator+(const ComplexType<T>&);
template<class T> ComplexType<T> operator-(const ComplexType<T>&);
template<class T> bool operator==(const ComplexType<T>&, const ComplexType<T>&);
template<class T> bool operator==(const ComplexType<T>&, const T&);
template<class T> bool operator==(const T&, const ComplexType<T>&);
template<class T> bool operator!=(const ComplexType<T>&, const ComplexType<T>&);
template<class T> bool operator!=(const ComplexType<T>&, const T&);
template<class T> bool operator!=(const T&, const ComplexType<T>&);

template<class T, class charT, class traits>
  basic_istream<charT, traits>&
  operator>>(basic_istream<charT, traits>&, ComplexType<T>&);
template<class T, class charT, class traits>
  basic_ostream<charT, traits>&
  operator<<(basic_ostream<charT, traits>&, const ComplexType<T>&);

// 26.3.7 values:

template<class T>              T real(const ComplexType<T>&);
                          double real(double);
template<Integral T>      double real(T);
                          float  real(float);

template<class T>              T imag(const ComplexType<T>&);
                          double imag(double);
template<Integral T>      double imag(T);
                          float  imag(float);

template<class T> T abs(const ComplexType<T>&);

template<class T>              T arg(const ComplexType<T>&);
                          double arg(double);
template<Integral T>      double arg(T);
                          float  arg(float);

template<class T>              T norm(const ComplexType<T>&);
                          double norm(double);
template<Integral T>      double norm(T);
                          float  norm(float);

template<class T>      ComplexType<T>           conj(const ComplexType<T>&);
                       ComplexType<double>      conj(double);
template<Integral T>   ComplexType<double>      conj(T);
                       ComplexType<float>       conj(float);

template<class T>    ComplexType<T>           proj(const ComplexType<T>&);
                     ComplexType<double>      proj(double);
template<Integral T> ComplexType<double>      proj(T);
                     ComplexType<float>       proj(float);

template<class T> ComplexType<T> polar(const T&, const T& = 0);

// 26.3.8 transcendentals:
template<class T> ComplexType<T> acos(const ComplexType<T>&);
template<class T> ComplexType<T> asin(const ComplexType<T>&);
template<class T> ComplexType<T> atan(const ComplexType<T>&);
template<class T> ComplexType<T> acosh(const ComplexType<T>&);
template<class T> ComplexType<T> asinh(const ComplexType<T>&);
template<class T> ComplexType<T> atanh(const ComplexType<T>&);
template<class T> ComplexType<T> cos (const ComplexType<T>&);
template<class T> ComplexType<T> cosh (const ComplexType<T>&);
template<class T> ComplexType<T> exp (const ComplexType<T>&);
template<class T> ComplexType<T> log (const ComplexType<T>&);
template<class T> ComplexType<T> log10(const ComplexType<T>&);

template<class T> ComplexType<T> pow(const ComplexType<T>&, const T&);
template<class T> ComplexType<T> pow(const ComplexType<T>&, const ComplexType<T>&);
template<class T> ComplexType<T> pow(const T&, const ComplexType<T>&);

template<class T> ComplexType<T> sin (const ComplexType<T>&);
template<class T> ComplexType<T> sinh (const ComplexType<T>&);
template<class T> ComplexType<T> sqrt (const ComplexType<T>&);
template<class T> ComplexType<T> tan (const ComplexType<T>&);
template<class T> ComplexType<T> tanh (const ComplexType<T>&);

template<class T, class charT, class traits>
  basic_istream<charT, traits>&
  operator>>(basic_istream<charT, traits>& is, ComplexType<T>& x);

template<class T, class charT, class traits>
  basic_ostream<charT, traits>&
  operator<<(basic_ostream<charT, traits>& o, const ComplexType<T>& x);

*/

#include <math.h>
#include <sstream>

template<class Value> class  ComplexType;

template<class Value> ComplexType<Value> operator*(const ComplexType<Value>& __z, const ComplexType<Value>& __w);
template<class Value> ComplexType<Value> operator/(const ComplexType<Value>& __x, const ComplexType<Value>& __y);

template<class Value>
class  ComplexType
{
private:
    Value re;
    Value im;
public:
    ComplexType(const Value __re = 0, const Value  __im = 0);

    ComplexType(const ComplexType<Value> &  __c);

    const Value& real() const;
    const Value& imag() const;

    void real(Value __re);
    void imag(Value __im);

    ComplexType& operator= (const Value& __re);

    ComplexType& operator+=(const Value& __re);

    ComplexType& operator-=(const Value& __re);
    ComplexType& operator*=(const Value& __re);
    ComplexType& operator/=(const Value& __re);

   ComplexType& operator= (const ComplexType<Value>& __c);

   ComplexType& operator+=(const ComplexType<Value>& __c);

    ComplexType& operator-=(const ComplexType<Value>& __c);

    ComplexType& operator*=(const ComplexType<Value>& __c);
    ComplexType& operator/=(const ComplexType<Value>& __c);

};


template <class Value>
inline int  isinf(Value x);


template <class Value>
inline int  isfinite(Value x);

// 26.3.6 operators:


double maxExpression(double a, double b);

template<class Value>
inline ComplexType<Value> operator+(const ComplexType<Value>& __x, const ComplexType<Value>& __y);

template<class Value>
inline ComplexType<Value> operator+(const ComplexType<Value>& __x, const Value& __y);

template<class Value>
inline ComplexType<Value> operator+(const Value& __x, const ComplexType<Value>& __y);
template<class Value>
inline ComplexType<Value> operator- (const ComplexType<Value>& __x, const ComplexType<Value>& __y);

template<class Value>
inline ComplexType<Value> operator-(const ComplexType<Value>& __x, const Value& __y);

template<class Value>
inline ComplexType<Value> operator-(const Value& __x, const ComplexType<Value>& __y);

template<class Value>
ComplexType<Value> operator*(const ComplexType<Value>& __z, const ComplexType<Value>& __w);

template<class Value>
inline ComplexType<Value> operator*(const ComplexType<Value>& __x, const Value& __y);

template<class Value>
inline ComplexType<Value> operator*(const Value& __x, const ComplexType<Value>& __y);


template<class Value>
ComplexType<Value> operator/(const ComplexType<Value>& __z, const ComplexType<Value>& __w);

template<class Value>
inline  ComplexType<Value> operator/(const ComplexType<Value>& __x, const Value& __y);

template<class Value>
inline ComplexType<Value> operator/(const Value& __x, const ComplexType<Value>& __y);

template<class Value>
inline ComplexType<Value> operator+(const ComplexType<Value>& __x);

template<class Value>
inline ComplexType<Value>
operator-(const ComplexType<Value>& __x);

template<class Value>
inline bool operator==(const ComplexType<Value>& __x, const ComplexType<Value>& __y);

template<class Value>
inline bool operator==(const ComplexType<Value>& __x, const Value& __y);

template<class Value>
inline bool operator==(const Value& __x, const ComplexType<Value>& __y);

template<class Value>
inline bool operator!=(const ComplexType<Value>& __x, const ComplexType<Value>& __y);

template<class Value>
inline bool operator!=(const ComplexType<Value>& __x, const Value& __y);

template<class Value>
inline bool operator!=(const Value& __x, const ComplexType<Value>& __y);

// 26.3.7 values:

// real

template<class Value>
inline Value real(const ComplexType<Value>& __c);


// imag

template<class Value>
inline Value imag(const ComplexType<Value>& __c);

inline double imag(double __re);

inline float imag(float  __re);

// abs

template<class Value>
inline Value abs(const ComplexType<Value>& __c);
// arg

template<class Value>
inline Value arg(const ComplexType<Value>& __c);

inline double arg(double __re);


// norm

template<class Value>
inline Value norm(const ComplexType<Value>& __c);

inline double norm(double __re);

// conj

template<class Value>
inline ComplexType<Value> conj(const ComplexType<Value>& __c);
// proj

template<class Value>
inline ComplexType<Value> proj(const ComplexType<Value>& __c);

// polar

template<class Value> 
ComplexType<Value> polar(const Value& __rho, const Value& __theta);


template<class Value> 
ComplexType<Value> polar(const Value& __rho, const Value& __theta)
{
    if (isnan(__rho) || signbit(__rho))
        return ComplexType<Value>(Value(NAN), Value(NAN));
    if (isnan(__theta))
    {
        if (isinf(__rho))
            return ComplexType<Value>(__rho, __theta);
        return ComplexType<Value>(__theta, __theta);
    }
    if (isinf(__theta))
    {
        if (isinf(__rho))
            return ComplexType<Value>(__rho, Value(NAN));
        return ComplexType<Value>(Value(NAN), Value(NAN));
    }
    Value __x = __rho * cos(__theta);
    if (isnan(__x))
        __x = 0;
    Value __y = __rho * sin(__theta);
    if (isnan(__y))
        __y = 0;
    return ComplexType<Value>(__x, __y);
}

// log

template<class Value>
inline ComplexType<Value> log(const ComplexType<Value>& __x);
// log10

template<class Value>
inline ComplexType<Value> log10(const ComplexType<Value>& __x);
// sqrt

template<class Value>
ComplexType<Value> sqrt(const ComplexType<Value>& __x);

// exp

template<class Value>
ComplexType<Value> exp(const ComplexType<Value>& __x);
// pow

template<class Value>
inline ComplexType<Value> pow(const ComplexType<Value>& __x, const ComplexType<Value>& __y);

template<class Value>
inline ComplexType<Value> pow(const ComplexType<Value>& __x, const Value& __y);

template<class Value>
inline ComplexType<Value> pow(const Value& __x, const ComplexType<Value>& __y);

// asinh

template<class Value>
ComplexType<Value> asinh(const ComplexType<Value>& __x);

// acosh

template<class Value>
ComplexType<Value> acosh(const ComplexType<Value>& __x);

// atanh

template<class Value>
ComplexType<Value> atanh(const ComplexType<Value>& __x);

// sinh

template<class Value>
ComplexType<Value> sinh(const ComplexType<Value>& __x);
// cosh

template<class Value>
ComplexType<Value> cosh(const ComplexType<Value>& __x);

// tanh

template<class Value>
ComplexType<Value> tanh(const ComplexType<Value>& __x);
// asin

template<class Value>
ComplexType<Value> asin(const ComplexType<Value>& __x);

// acos

template<class Value>
ComplexType<Value> acos(const ComplexType<Value>& __x);

// atan

template<class Value>
ComplexType<Value> atan(const ComplexType<Value>& __x);

// sin

template<class Value>
ComplexType<Value> sin(const ComplexType<Value>& __x);
// cos

template<class Value>
inline ComplexType<Value> cos(const ComplexType<Value>& __x);

template<class Value>
ComplexType<Value> tan(const ComplexType<Value>& __x);

template<class Value, class _CharT, class _Traits>
std::basic_istream<_CharT, _Traits>& operator>>(std::basic_istream<_CharT, _Traits>& __is, ComplexType<Value>& __x);

template<class Value, class _CharT, class _Traits>
std::basic_ostream<_CharT, _Traits>& operator<<(std::basic_ostream<_CharT, _Traits>& __os, const ComplexType<Value>& __x);


template<class Value>
ComplexType<Value>::ComplexType(const Value __re, const Value __im)
	: re(__re), im(__im) {}

template<class Value>
ComplexType<Value>::ComplexType(const ComplexType<Value> & __c)
	: re(__c.real()), im(__c.imag()) {}

template<class Value>
const Value& ComplexType<Value>::real() const {return re;}
template<class Value>
const Value& ComplexType<Value>::imag() const {return im;}
template<class Value>
void ComplexType<Value>::real(Value __re) {re = __re;}
template<class Value>
void ComplexType<Value>::imag(Value __im) {im = __im;}
template<class Value>
ComplexType<Value>& ComplexType<Value>::operator= (const Value& __re)
{re = __re; im = 0; return *this;}
template<class Value>
ComplexType<Value>& ComplexType<Value>::operator+=(const Value& __re) {re += __re; return *this;}
template<class Value>
ComplexType<Value>& ComplexType<Value>::operator-=(const Value& __re) {re -= __re; return *this;}
template<class Value>
ComplexType<Value>& ComplexType<Value>::operator*=(const Value& __re) {re *= __re; im *= __re; return *this;}
template<class Value>
ComplexType<Value>& ComplexType<Value>::operator/=(const Value& __re) {re /= __re; im /= __re; return *this;}

template<class Value> ComplexType<Value>& ComplexType<Value>::operator= (const ComplexType<Value>& __c)
{
	re = __c.real();
	im = __c.imag();
	return *this;
}

template<class Value> ComplexType<Value>& ComplexType<Value>::operator+=(const ComplexType<Value>& __c)
{
	re += __c.real();
	im += __c.imag();
	return *this;
}

template<class Value> ComplexType<Value>& ComplexType<Value>::operator-=(const ComplexType<Value>& __c)
{
	re -= __c.real();
	im -= __c.imag();
	return *this;
}

template<class Value> ComplexType<Value>& ComplexType<Value>::operator*=(const ComplexType<Value>& __c)
{
	*this = *this * __c;
	return *this;
}

template<class Value> ComplexType<Value>& ComplexType<Value>::operator/=(const ComplexType<Value>& __c)
{
	*this = *this / __c;
	return *this;
}




template <class Value>
inline int  isinf(Value x) { return (x == std::numeric_limits<Value>::infinity()); };


template <class Value>
inline int  isfinite(Value x) { return (x != std::numeric_limits<Value>::infinity()); };

// 26.3.6 operators:


double maxExpression(double a, double b)
{
	return (a > b) ? a : b;
};

template<class Value>
inline ComplexType<Value> operator+(const ComplexType<Value>& __x, const ComplexType<Value>& __y)
{
	ComplexType<Value> __t(__x);
	__t += __y;
	return __t;
}

template<class Value>
inline ComplexType<Value> operator+(const ComplexType<Value>& __x, const Value& __y)
{
	ComplexType<Value> __t(__x);
	__t += __y;
	return __t;
}

template<class Value>
inline ComplexType<Value> operator+(const Value& __x, const ComplexType<Value>& __y)
{
	ComplexType<Value> __t(__y);
	__t += __x;
	return __t;
}

template<class Value>
inline ComplexType<Value> operator- (const ComplexType<Value>& __x, const ComplexType<Value>& __y)
{
	ComplexType<Value> __t(__x);
	__t -= __y;
	return __t;
}

template<class Value>
inline ComplexType<Value> operator-(const ComplexType<Value>& __x, const Value& __y)
{
	ComplexType<Value> __t(__x);
	__t -= __y;
	return __t;
}

template<class Value>
inline ComplexType<Value> operator-(const Value& __x, const ComplexType<Value>& __y)
{
	ComplexType<Value> __t(-__y);
	__t += __x;
	return __t;
}

template<class Value>
ComplexType<Value> operator*(const ComplexType<Value>& __z, const ComplexType<Value>& __w)
{
	Value __a = __z.real();
	Value __b = __z.imag();
	Value __c = __w.real();
	Value __d = __w.imag();
	Value __ac = __a * __c;
	Value __bd = __b * __d;
	Value __ad = __a * __d;
	Value __bc = __b * __c;
	Value __x = __ac - __bd;
	Value __y = __ad + __bc;
	/*    if (isnan(__x) && isnan(__y))
	{
	bool __recalc = false;
	if (isinf(__a) || isinf(__b))
	{
	__a = copysign(isinf(__a) ? Value(1, 1) : Value(0), __a);
	__b = copysign(isinf(__b) ? Value(1, 1) : Value(0), __b);
	if (isnan(__c))
	__c = copysign(Value(0), __c);
	if (isnan(__d))
	__d = copysign(Value(0), __d);
	__recalc = true;
	}
	if (isinf(__c) || isinf(__d))
	{
	__c = copysign(isinf(__c) ? Value(1, 1) : Value(0), __c);
	__d = copysign(isinf(__d) ? Value(1, 1) : Value(0), __d);
	if (isnan(__a))
	__a = copysign(Value(0), __a);
	if (isnan(__b))
	__b = copysign(Value(0), __b);
	__recalc = true;
	}
	if (!__recalc && (isinf(__ac) || isinf(__bd) ||
	isinf(__ad) || isinf(__bc)))
	{
	if (isnan(__a))
	__a = copysign(Value(0), __a);
	if (isnan(__b))
	__b = copysign(Value(0), __b);
	if (isnan(__c))
	__c = copysign(Value(0), __c);
	if (isnan(__d))
	__d = copysign(Value(0), __d);
	__recalc = true;
	}
	if (__recalc)
	{
	__x = Value(M_INFINITY, M_INFINITY) * (__a * __c - __b * __d);
	__y = Value(M_INFINITY, M_INFINITY) * (__a * __d + __b * __c);
	}
	}*/
	return ComplexType<Value>(__x, __y);
}

template<class Value>
inline ComplexType<Value> operator*(const ComplexType<Value>& __x, const Value& __y)
{
	ComplexType<Value> __t(__x);
	__t *= __y;
	return __t;
}

template<class Value>
inline ComplexType<Value> operator*(const Value& __x, const ComplexType<Value>& __y)
{
	ComplexType<Value> __t(__y);
	__t *= __x;
	return __t;
}


template<class Value>
ComplexType<Value> operator/(const ComplexType<Value>& __z, const ComplexType<Value>& __w)
{
	int __ilogbw = 0;
	Value __a = __z.real();
	Value __b = __z.imag();
	Value __c = __w.real();
	Value __d = __w.imag();
	Value __logbw = log(maxExpression(fabs(__c), fabs(__d)));
	/*   if (isfinite(__logbw))
	{
	__ilogbw = static_cast<int>(__logbw);
	__c = scalbn(__c, -__ilogbw);
	__d = scalbn(__d, -__ilogbw);
	}*/
	Value __denom = __c * __c + __d * __d;
	Value __x = (__a * __c + __b * __d) / __denom;
	Value __y = (__b * __c - __a * __d) / __denom;
	/*    if (isnan(__x) && isnan(__y))
	{
	if ((__denom == Value(0)) && (!isnan(__a) || !isnan(__b)))
	{
	__x = copysign(Value(M_INFINITY), __c) * __a;
	__y = copysign(Value(M_INFINITY), __c) * __b;
	}
	else if ((isinf(__a) || isinf(__b)) && isfinite(__c) && isfinite(__d))
	{
	__a = copysign(isinf(__a) ? Value(1) : Value(0), __a);
	__b = copysign(isinf(__b) ? Value(1) : Value(0), __b);
	__x = Value(M_INFINITY) * (__a * __c + __b * __d);
	__y = Value(M_INFINITY) * (__b * __c - __a * __d);
	}
	else if (isinf(__logbw) && __logbw > Value(0) && isfinite(__a) && isfinite(__b))
	{
	__c = copysign(isinf(__c) ? Value(1) : Value(0), __c);
	__d = copysign(isinf(__d) ? Value(1) : Value(0), __d);
	__x = Value(0) * (__a *is __c + __b * __d);
	__y = Value(0) * (__b * __c - __a * __d);
	}
	}*/
	return ComplexType<Value>(__x, __y);
}

template<class Value>
inline  ComplexType<Value> operator/(const ComplexType<Value>& __x, const Value& __y)
{
	return ComplexType<Value>(__x.real() / __y, __x.imag() / __y);
}

template<class Value>
inline ComplexType<Value> operator/(const Value& __x, const ComplexType<Value>& __y)
{
	ComplexType<Value> __t(__x);
	__t /= __y;
	return __t;
}

template<class Value>
inline ComplexType<Value> operator+(const ComplexType<Value>& __x)
{
	return __x;
}

template<class Value>
inline ComplexType<Value>
	operator-(const ComplexType<Value>& __x)
{
	return ComplexType<Value>(-__x.real(), -__x.imag());
}

template<class Value>
inline bool operator==(const ComplexType<Value>& __x, const ComplexType<Value>& __y)
{
	return __x.real() == __y.real() && __x.imag() == __y.imag();
}

template<class Value>
inline bool operator==(const ComplexType<Value>& __x, const Value& __y)
{
	return __x.real() == __y && __x.imag() == 0;
}

template<class Value>
inline bool operator==(const Value& __x, const ComplexType<Value>& __y)
{
	return __x == __y.real() && 0 == __y.imag();
}

template<class Value>
inline bool operator!=(const ComplexType<Value>& __x, const ComplexType<Value>& __y)
{
	return !(__x == __y);
}

template<class Value>
inline bool operator!=(const ComplexType<Value>& __x, const Value& __y)
{
	return !(__x == __y);
}

template<class Value>
inline bool operator!=(const Value& __x, const ComplexType<Value>& __y)
{
	return !(__x == __y);
}

// 26.3.7 values:

// real

template<class Value>
inline Value real(const ComplexType<Value>& __c)
{
	return __c.real();
}


// imag

template<class Value>
inline Value imag(const ComplexType<Value>& __c)
{
	return __c.imag();
}

inline double imag(double __re)
{
	return 0;
}

inline float imag(float  __re)
{
	return 0;
}

// abs

template<class Value>
inline Value abs(const ComplexType<Value>& __c)
{
	return sqrt(__c.real()*__c.real()+__c.imag()*__c.imag());
}

// arg

template<class Value>
inline Value arg(const ComplexType<Value>& __c)
{
	return atan2(__c.imag(), __c.real());
}

inline double arg(double __re)
{
	return atan2(0., __re);
}



// norm

template<class Value>
inline Value norm(const ComplexType<Value>& __c)
{
	if (isinf(__c.real()))
		return fabs(__c.real());
	if (isinf(__c.imag()))
		return fabs(__c.imag());
	return __c.real() * __c.real() + __c.imag() * __c.imag();
}

inline double norm(double __re)
{
	return __re * __re;
}



// conj

template<class Value>
inline ComplexType<Value> conj(const ComplexType<Value>& __c)
{
	return ComplexType<Value>(__c.real(), -__c.imag());
}

// proj

template<class Value>
inline ComplexType<Value> proj(const ComplexType<Value>& __c)
{
	ComplexType<Value> __r = __c;
	if (isinf(__c.real()) || isinf(__c.imag()))
		__r = ComplexType<Value>(M_INFINITY, copysign(Value(0), __c.imag()));
	return __r;
}

// polar


// log

template<class Value>
inline ComplexType<Value> log(const ComplexType<Value>& __x)
{
	return ComplexType<Value>(log(abs(__x)), arg(__x));
}

// log10

template<class Value>
inline ComplexType<Value> log10(const ComplexType<Value>& __x)
{
	return log(__x) / log(10);
}

// sqrt

template<class Value>
ComplexType<Value> sqrt(const ComplexType<Value>& __x)
{
	if (isinf(__x.imag()))
		return ComplexType<Value>(Value(M_INFINITY), __x.imag());
	if (isinf(__x.real()))
	{
		if (__x.real() > 0)
			return ComplexType<Value>(__x.real(), isnan(__x.imag()) ? __x.imag() : copysign(0, __x.imag()));
		return ComplexType<Value>(isnan(__x.imag()) ? __x.imag() : 0, copysign(__x.real(), __x.imag()));
	}
	return polar(sqrt(abs(__x)), arg(__x) / Value(2));
}

// exp

template<class Value>
ComplexType<Value> exp(const ComplexType<Value>& __x)
{
	Value __i = __x.imag();
	/*  if (isinf(__x.real()))
	{
	if (__x.real() < Value(0))
	{
	if (!isfinite(__i))
	__i = Value(1, 1);
	}
	else if (__i == 0 || !isfinite(__i))
	{
	if (isinf(__i))
	__i = Value(NAN);
	return ComplexType<Value>(__x.real(), __i);
	}
	}
	else if (isnan(__x.real()) && __x.imag() == 0)
	return __x;*/
	Value __e = exp(__x.real());
	ComplexType<Value> res = ComplexType<Value>(__e * cos(__i), __e * sin(__i));
	return res;
}

// pow

template<class Value>
inline ComplexType<Value> pow(const ComplexType<Value>& __x, const ComplexType<Value>& __y)
{
	return exp(__y * log(__x));
}

template<class Value>
inline ComplexType<Value> pow(const ComplexType<Value>& __x, const Value& __y)
{   
	return pow(__x, ComplexType<Value>(__y));
}

template<class Value>
inline ComplexType<Value> pow(const Value& __x, const ComplexType<Value>& __y)
{
	return pow(ComplexType<Value>(__x), __y);
}

// asinh

template<class Value>
ComplexType<Value> asinh(const ComplexType<Value>& __x)
{
	const Value __pi(atan2(+0., -0.));
	if (isinf(__x.real()))
	{
		if (isnan(__x.imag()))
			return __x;
		if (isinf(__x.imag()))
			return ComplexType<Value>(__x.real(), copysign(__pi * Value(0.25), __x.imag()));
		return ComplexType<Value>(__x.real(), copysign(Value(0), __x.imag()));
	}
	if (isnan(__x.real()))
	{
		if (isinf(__x.imag()))
			return ComplexType<Value>(__x.imag(), __x.real());
		if (__x.imag() == 0)
			return __x;
		return ComplexType<Value>(__x.real(), __x.real());
	}
	if (isinf(__x.imag()))
		return ComplexType<Value>(copysign(__x.imag(), __x.real()), copysign(__pi/Value(2), __x.imag()));
	ComplexType<Value> __z = log(__x + sqrt(pow(__x, Value(2)) + Value(1)));
	return ComplexType<Value>(copysign(__z.real(), __x.real()), copysign(__z.imag(), __x.imag()));
}

// acosh

template<class Value>
ComplexType<Value> acosh(const ComplexType<Value>& __x)
{
	const Value __pi(atan2(+0., -0.));
	if (isinf(__x.real()))
	{
		if (isnan(__x.imag()))
			return ComplexType<Value>(fabs(__x.real()), __x.imag());
		if (isinf(__x.imag()))
			if (__x.real() > 0)
				return ComplexType<Value>(__x.real(), copysign(__pi * Value(0.25), __x.imag()));
			else
				return ComplexType<Value>(-__x.real(), copysign(__pi * Value(0.75), __x.imag()));
		if (__x.real() < 0)
			return ComplexType<Value>(-__x.real(), copysign(__pi, __x.imag()));
		return ComplexType<Value>(__x.real(), copysign(Value(0), __x.imag()));
	}
	if (isnan(__x.real()))
	{
		if (isinf(__x.imag()))
			return ComplexType<Value>(fabs(__x.imag()), __x.real());
		return ComplexType<Value>(__x.real(), __x.real());
	}
	if (isinf(__x.imag()))
		return ComplexType<Value>(fabs(__x.imag()), copysign(__pi/Value(2), __x.imag()));
	ComplexType<Value> __z = log(__x + sqrt(pow(__x, Value(2)) - Value(1)));
	return ComplexType<Value>(copysign(__z.real(), Value(0)), copysign(__z.imag(), __x.imag()));
}

// atanh

template<class Value>
ComplexType<Value> atanh(const ComplexType<Value>& __x)
{
	const Value __pi(atan2(+0., -0.));
	if (isinf(__x.imag()))
	{
		return ComplexType<Value>(copysign(Value(0), __x.real()), copysign(__pi/Value(2), __x.imag()));
	}
	if (isnan(__x.imag()))
	{
		if (isinf(__x.real()) || __x.real() == 0)
			return ComplexType<Value>(copysign(Value(0), __x.real()), __x.imag());
		return ComplexType<Value>(__x.imag(), __x.imag());
	}
	if (isnan(__x.real()))
	{
		return ComplexType<Value>(__x.real(), __x.real());
	}
	if (isinf(__x.real()))
	{
		return ComplexType<Value>(copysign(Value(0), __x.real()), copysign(__pi/Value(2), __x.imag()));
	}
	if (fabs(__x.real()) == Value(1) && __x.imag() == Value(0))
	{
		return ComplexType<Value>(copysign(Value(M_INFINITY), __x.real()), copysign(Value(0), __x.imag()));
	}
	ComplexType<Value> __z = log((Value(1) + __x) / (Value(1) - __x)) / Value(2);
	return ComplexType<Value>(copysign(__z.real(), __x.real()), copysign(__z.imag(), __x.imag()));
}

// sinh

template<class Value>
ComplexType<Value> sinh(const ComplexType<Value>& __x)
{
	if (isinf(__x.real()) && !isfinite(__x.imag()))
		return ComplexType<Value>(__x.real(), Value(NAN));
	if (__x.real() == 0 && !isfinite(__x.imag()))
		return ComplexType<Value>(__x.real(), Value(NAN));
	if (__x.imag() == 0 && !isfinite(__x.real()))
		return __x;
	return ComplexType<Value>(sinh(__x.real()) * cos(__x.imag()), cosh(__x.real()) * sin(__x.imag()));
}

// cosh

template<class Value>
ComplexType<Value> cosh(const ComplexType<Value>& __x)
{
	if (isinf(__x.real()) && !isfinite(__x.imag()))
		return ComplexType<Value>(fabs(__x.real()), Value(NAN));
	if (__x.real() == 0 && !isfinite(__x.imag()))
		return ComplexType<Value>(Value(NAN), __x.real());
	if (__x.real() == 0 && __x.imag() == 0)
		return ComplexType<Value>(Value(1), __x.imag());
	if (__x.imag() == 0 && !isfinite(__x.real()))
		return ComplexType<Value>(fabs(__x.real()), __x.imag());
	return ComplexType<Value>(cosh(__x.real()) * cos(__x.imag()), sinh(__x.real()) * sin(__x.imag()));
}

// tanh

template<class Value>
ComplexType<Value> tanh(const ComplexType<Value>& __x)
{
	if (isinf(__x.real()))
	{
		if (!isfinite(__x.imag()))
			return ComplexType<Value>(Value(1), Value(0));
		return ComplexType<Value>(Value(1), copysign(Value(0), sin(Value(2) * __x.imag())));
	}
	if (isnan(__x.real()) && __x.imag() == 0)
		return __x;
	Value __2r(Value(2) * __x.real());
	Value __2i(Value(2) * __x.imag());
	Value __d(cosh(__2r) + cos(__2i));
	return  ComplexType<Value>(sinh(__2r)/__d, sin(__2i)/__d);
}

// asin

template<class Value>
ComplexType<Value> asin(const ComplexType<Value>& __x)
{
	ComplexType<Value> __z = asinh(ComplexType<Value>(-__x.imag(), __x.real()));
	return ComplexType<Value>(__z.imag(), -__z.real());
}

// acos

template<class Value>
ComplexType<Value> acos(const ComplexType<Value>& __x)
{
	const Value __pi(atan2(+0., -0.));
	if (isinf(__x.real()))
	{
		if (isnan(__x.imag()))
			return ComplexType<Value>(__x.imag(), __x.real());
		if (isinf(__x.imag()))
		{
			if (__x.real() < Value(0))
				return ComplexType<Value>(Value(0.75) * __pi, -__x.imag());
			return ComplexType<Value>(Value(0.25) * __pi, -__x.imag());
		}
		if (__x.real() < Value(0))
			return ComplexType<Value>(__pi, signbit(__x.imag()) ? -__x.real() : __x.real());
		return ComplexType<Value>(Value(0), signbit(__x.imag()) ? __x.real() : -__x.real());
	}
	if (isnan(__x.real()))
	{
		if (isinf(__x.imag()))
			return ComplexType<Value>(__x.real(), -__x.imag());
		return ComplexType<Value>(__x.real(), __x.real());
	}
	if (isinf(__x.imag()))
		return ComplexType<Value>(__pi/Value(2), -__x.imag());
	if (__x.real() == 0)
		return ComplexType<Value>(__pi/Value(2), -__x.imag());
	ComplexType<Value> __z = log(__x + sqrt(pow(__x, Value(2)) - Value(1)));
	if (signbit(__x.imag()))
		return ComplexType<Value>(fabs(__z.imag()), fabs(__z.real()));
	return ComplexType<Value>(fabs(__z.imag()), -fabs(__z.real()));
}

// atan

template<class Value>
ComplexType<Value> atan(const ComplexType<Value>& __x)
{
	ComplexType<Value> __z = atanh(ComplexType<Value>(-__x.imag(), __x.real()));
	return ComplexType<Value>(__z.imag(), -__z.real());
}

// sin

template<class Value>
ComplexType<Value> sin(const ComplexType<Value>& __x)
{
	ComplexType<Value> __z = sinh(ComplexType<Value>(-__x.imag(), __x.real()));
	return ComplexType<Value>(__z.imag(), -__z.real());
}

// cos

template<class Value>
inline ComplexType<Value> cos(const ComplexType<Value>& __x)
{
	return cosh(ComplexType<Value>(-__x.imag(), __x.real()));
}

// tan

template<class Value>
ComplexType<Value> tan(const ComplexType<Value>& __x)
{
	ComplexType<Value> __z = tanh(ComplexType<Value>(-__x.imag(), __x.real()));
	return ComplexType<Value>(__z.imag(), -__z.real());
}

template<class Value, class _CharT, class _Traits>
std::basic_istream<_CharT, _Traits>& operator>>(std::basic_istream<_CharT, _Traits>& __is, ComplexType<Value>& __x)
{
	if (__is.good())
	{
		ws(__is);
		if (__is.peek() == _CharT('('))
		{
			__is.get();
			Value __r;
			__is >> __r;
			if (!__is.fail())
			{
				ws(__is);
				_CharT __c = __is.peek();
				if (__c == _CharT(','))
				{
					__is.get();
					Value __i;
					__is >> __i;
					if (!__is.fail())
					{
						ws(__is);
						__c = __is.peek();
						if (__c == _CharT(')'))
						{
							__is.get();
							__x = ComplexType<Value>(__r, __i);
						}
						else
							__is.setstate(std::ios_base::failbit);
					}
					else
						__is.setstate(std::ios_base::failbit);
				}
				else if (__c == _CharT(')'))
				{
					__is.get();
					__x = ComplexType<Value>(__r, Value(0));
				}
				else
					__is.setstate(std::ios_base::failbit);
			}
			else
				__is.setstate(std::ios_base::failbit);
		}
		else
		{
			Value __r;
			__is >> __r;
			if (!__is.fail())
				__x = ComplexType<Value>(__r, Value(0));
			else
				__is.setstate(std::ios_base::failbit);
		}
	}
	else
		__is.setstate(std::ios_base::failbit);
	return __is;
}

template<class Value, class _CharT, class _Traits>
std::basic_ostream<_CharT, _Traits>& operator<<(std::basic_ostream<_CharT, _Traits>& __os, const ComplexType<Value>& __x)
{
	std::basic_ostringstream<_CharT, _Traits> __s;
	__s.flags(__os.flags());
	__s.imbue(__os.getloc());
	__s.precision(__os.precision());
	__s << '(' << __x.real() << ',' << __x.imag() << ')';
	return __os << __s.str();
}


#endif  // COMPLEX_H
