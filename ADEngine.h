#ifndef ADENGINE_H
#define ADENGINE_H

#ifndef NUM_VARIABLES
#define NUM_VARIABLES 13
#endif

#ifndef FUNC_INLINE
#define FUNC_INLINE inline
#endif

#define MAX_MACRO(a, b) (((a) > (b)) ? (a) : (b))
#define MIN_MACRO(a, b) (((a) < (b)) ? (a) : (b))

#ifndef DATA_TYPE
#define DATA_TYPE double
#endif

DATA_TYPE values[NUM_VARIABLES];
DATA_TYPE sensitivities[NUM_VARIABLES];


// Variable

template <int index>
class Variable
{
public:
	Variable();
	Variable(DATA_TYPE v);
	~Variable();

	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE DATA_TYPE diff();
	Variable<index> operator=(Variable<index> rhs);

	DATA_TYPE value;
	DATA_TYPE accum;
};

template <int index>
Variable<index>::Variable() 
{ 
	value = values[index];
	accum = sensitivities[index];
};

template <int index>
Variable<index>::Variable(DATA_TYPE v): value(v), accum(sensitivities[index]) 
{
	values[index] = v;
};

template <int index>
Variable<index>::~Variable() 
{
}

template <int index>
FUNC_INLINE DATA_TYPE Variable<index>::evaluate()
{
	//cout << "value[" << index << "]:=" << value << endl;
	return value;
};

template <int index>
FUNC_INLINE void Variable<index>::setDiff()
{
	accum = 1;
}

template <int index>
FUNC_INLINE DATA_TYPE Variable<index>::diff()
{
	sensitivities[index] += accum;
	//	cout << "accum[" << index << "]:=" << accum << endl;
	return accum;
}

template <int index>
Variable<index> Variable<index>::operator=(Variable<index> rhs)
{
	this->value = rhs.value;
	values[index] = rhs.value;
	return *this;
}

// Addition

template <class leftOperand, class rightOperand>
class Add
{
public:
	Add();
	Add(leftOperand arg1, rightOperand arg2);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	friend Add<leftOperand, rightOperand> operator+(leftOperand leftOp, rightOperand rightOp);

	leftOperand leftOp;
	rightOperand rightOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand, class rightOperand>
Add<leftOperand, rightOperand>::Add()
{
	value = 0;
	accum = 0;
}

template <class leftOperand, class rightOperand>
Add<leftOperand, rightOperand>::Add(leftOperand arg1, rightOperand arg2)
{
	leftOp = arg1;
	rightOp = arg2;
	accum = 0;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE DATA_TYPE Add<leftOperand, rightOperand>::evaluate()
{
	value = leftOp.evaluate() + rightOp.evaluate();
	return value;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE void Add<leftOperand, rightOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE void Add<leftOperand, rightOperand>::diff()
{
	leftOp.accum += accum;
	rightOp.accum += accum;
	leftOp.diff();
	rightOp.diff();
}

template <class leftOperand, class rightOperand>
Add<leftOperand, rightOperand> operator+(leftOperand leftOp, rightOperand rightOp)
{
	return Add<leftOperand, rightOperand>(leftOp, rightOp);
}


// Subtraction

template <class leftOperand, class rightOperand>
class Sub
{
public:
	Sub();
	Sub(leftOperand arg1, rightOperand arg2);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	friend Sub<leftOperand, rightOperand> operator-(leftOperand leftOp, rightOperand rightOp);

	leftOperand leftOp;
	rightOperand rightOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand, class rightOperand>
Sub<leftOperand, rightOperand>::Sub()
{
	value = 0;
	accum = 0;
}

template <class leftOperand, class rightOperand>
Sub<leftOperand, rightOperand>::Sub(leftOperand arg1, rightOperand arg2)
{
	leftOp = arg1;
	rightOp = arg2;
	accum = 0;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE DATA_TYPE Sub<leftOperand, rightOperand>::evaluate()
{
	value = leftOp.evaluate() - rightOp.evaluate();
	return value;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE void Sub<leftOperand, rightOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE void Sub<leftOperand, rightOperand>::diff()
{
	leftOp.accum += accum;
	rightOp.accum -= accum;
	leftOp.diff();
	rightOp.diff();
}

template <class leftOperand, class rightOperand>
Sub<leftOperand, rightOperand> operator-(leftOperand leftOp, rightOperand rightOp)
{
	return Sub<leftOperand, rightOperand>(leftOp, rightOp);
}

// Multiplication

template <class leftOperand, class rightOperand>
class Mul
{
public:
	Mul();
	Mul(leftOperand arg1, rightOperand arg2);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	friend Mul<leftOperand, rightOperand> operator*(leftOperand leftOp, rightOperand rightOp);

	leftOperand leftOp;
	rightOperand rightOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand, class rightOperand>
Mul<leftOperand, rightOperand>::Mul()
{
	value = 0;
	accum = 0;
}

template <class leftOperand, class rightOperand>
Mul<leftOperand, rightOperand>::Mul(leftOperand arg1, rightOperand arg2)
{
	leftOp = arg1;
	rightOp = arg2;
	accum = 0;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE DATA_TYPE Mul<leftOperand, rightOperand>::evaluate()
{
	value = leftOp.evaluate() * rightOp.evaluate();
	return value;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE void Mul<leftOperand, rightOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE void Mul<leftOperand, rightOperand>::diff()
{
	leftOp.accum += rightOp.value*accum;
	rightOp.accum += leftOp.value*accum;
	leftOp.diff();
	rightOp.diff();
}

template <class leftOperand, class rightOperand>
Mul<leftOperand, rightOperand> operator*(leftOperand leftOp, rightOperand rightOp)
{
	return Mul<leftOperand, rightOperand>(leftOp, rightOp);
}

// Division

template <class leftOperand, class rightOperand>
class Div
{
public:
	Div();
	Div(leftOperand arg1, rightOperand arg2);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	friend Div<leftOperand, rightOperand> operator*(leftOperand leftOp, rightOperand rightOp);

	leftOperand leftOp;
	rightOperand rightOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand, class rightOperand>
Div<leftOperand, rightOperand>::Div()
{
	value = 0;
	accum = 0;
}

template <class leftOperand, class rightOperand>
Div<leftOperand, rightOperand>::Div(leftOperand arg1, rightOperand arg2)
{
	leftOp = arg1;
	rightOp = arg2;
	accum = 0;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE DATA_TYPE Div<leftOperand, rightOperand>::evaluate()
{
	value = leftOp.evaluate() / rightOp.evaluate();
	return value;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE void Div<leftOperand, rightOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE void Div<leftOperand, rightOperand>::diff()
{
	leftOp.accum += accum/rightOp.value;
	rightOp.accum -= value/leftOp.value;
	leftOp.diff();
	rightOp.diff();
}

template <class leftOperand, class rightOperand>
Div<leftOperand, rightOperand> operator/(leftOperand leftOp, rightOperand rightOp)
{
	return Div<leftOperand, rightOperand>(leftOp, rightOp);
}


// Max

template <class leftOperand, class rightOperand>
class Max
{
public:
	Max();
	Max(leftOperand arg1, rightOperand arg2);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	leftOperand leftOp;
	rightOperand rightOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand, class rightOperand>
Max<leftOperand, rightOperand>::Max()
{
	value = 0;
	accum = 0;
}

template <class leftOperand, class rightOperand>
Max<leftOperand, rightOperand>::Max(leftOperand arg1, rightOperand arg2)
{
	leftOp = arg1;
	rightOp = arg2;
	accum = 0;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE DATA_TYPE Max<leftOperand, rightOperand>::evaluate()
{
	value = MAX_MACRO(leftOp.evaluate(), rightOp.evaluate());
	return value;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE void Max<leftOperand, rightOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE void Max<leftOperand, rightOperand>::diff()
{
	if(leftOp.value > rightOp.value)
	{
		leftOp.accum += accum;
		rightOp.accum += 0;
	}
	else
	{
		leftOp.accum += 0;
		rightOp.accum += accum;
	}
	leftOp.diff();
	rightOp.diff();
};

template <class leftOperand, class rightOperand>
Max<leftOperand, rightOperand> max2(leftOperand leftOp, rightOperand rightOp)
{
	return Max<leftOperand, rightOperand>(leftOp, rightOp);
}

// Min

template <class leftOperand, class rightOperand>
class Min
{
public:
	Min();
	Min(leftOperand arg1, rightOperand arg2);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	leftOperand leftOp;
	rightOperand rightOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand, class rightOperand>
Min<leftOperand, rightOperand>::Min()
{
	value = 0;
	accum = 0;
}

template <class leftOperand, class rightOperand>
Min<leftOperand, rightOperand>::Min(leftOperand arg1, rightOperand arg2)
{
	leftOp = arg1;
	rightOp = arg2;
	accum = 0;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE DATA_TYPE Min<leftOperand, rightOperand>::evaluate()
{
	value = MIN(leftOp.evaluate(), rightOp.evaluate());
	return value;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE void Min<leftOperand, rightOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand, class rightOperand>
FUNC_INLINE void Min<leftOperand, rightOperand>::diff()
{
	if(leftOp.value > rightOp.value)
	{
		leftOp.accum += 0;
		rightOp.accum += accum;
	}
	else
	{
		leftOp.accum += accum;
		rightOp.accum += 0;
	}
	leftOp.diff();
	rightOp.diff();
};

template <class leftOperand, class rightOperand>
Min<leftOperand, rightOperand> min2(leftOperand leftOp, rightOperand rightOp)
{
	return Min<leftOperand, rightOperand>(leftOp, rightOp);
}

// Sin

template <class leftOperand>
class Sin
{
public:
	Sin();
	Sin(leftOperand arg1);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	leftOperand leftOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand>
Sin<leftOperand>::Sin()
{
	value = 0;
	accum = 0;
}

template <class leftOperand>
Sin<leftOperand>::Sin(leftOperand arg1)
{
	leftOp = arg1;
	accum = 0;
}

template <class leftOperand>
FUNC_INLINE DATA_TYPE Sin<leftOperand>::evaluate()
{
	value = sin(leftOp.evaluate());
	return value;
}

template <class leftOperand>
FUNC_INLINE void Sin<leftOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand>
FUNC_INLINE void Sin<leftOperand>::diff()
{
	leftOp.accum += cos(leftOp.value)*accum;
	leftOp.diff();
};

template <class leftOperand>
Sin<leftOperand> sin(leftOperand leftOp)
{
	return Sin<leftOperand>(leftOp);
}

// Cos

template <class leftOperand>
class Cos
{
public:
	Cos();
	Cos(leftOperand arg1);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	leftOperand leftOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand>
Cos<leftOperand>::Cos()
{
	value = 0;
	accum = 0;
}

template <class leftOperand>
Cos<leftOperand>::Cos(leftOperand arg1)
{
	leftOp = arg1;
	accum = 0;
}

template <class leftOperand>
FUNC_INLINE DATA_TYPE Cos<leftOperand>::evaluate()
{
	value = cos(leftOp.evaluate());
	return value;
}

template <class leftOperand>
FUNC_INLINE void Cos<leftOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand>
FUNC_INLINE void Cos<leftOperand>::diff()
{
	leftOp.accum -= sin(leftOp.value)*accum;
	leftOp.diff();
};

template <class leftOperand>
Cos<leftOperand> cos(leftOperand leftOp)
{
	return Cos<leftOperand>(leftOp);
}

// Tan

template <class leftOperand>
class Tan
{
public:
	Tan();
	Tan(leftOperand arg1);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	leftOperand leftOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand>
Tan<leftOperand>::Tan()
{
	value = 0;
	accum = 0;
}

template <class leftOperand>
Tan<leftOperand>::Tan(leftOperand arg1)
{
	leftOp = arg1;
	accum = 0;
}

template <class leftOperand>
FUNC_INLINE DATA_TYPE Tan<leftOperand>::evaluate()
{
	value = tan(leftOp.evaluate());
	return value;
}

template <class leftOperand>
FUNC_INLINE void Tan<leftOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand>
FUNC_INLINE void Tan<leftOperand>::diff()
{
	leftOp.accum += (1/(cos(leftOp.value)*cos(leftOp.value)))*accum;
	leftOp.diff();
};

template <class leftOperand>
Tan<leftOperand> tan(leftOperand leftOp)
{
	return Tan<leftOperand>(leftOp);
}

// Ctg

template <class leftOperand>
class Ctg
{
public:
	Ctg();
	Ctg(leftOperand arg1);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	leftOperand leftOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand>
Ctg<leftOperand>::Ctg()
{
	value = 0;
	accum = 0;
}

template <class leftOperand>
Ctg<leftOperand>::Ctg(leftOperand arg1)
{
	leftOp = arg1;
	accum = 0;
}

template <class leftOperand>
FUNC_INLINE DATA_TYPE Ctg<leftOperand>::evaluate()
{
	value = 1/tan(leftOp.evaluate());
	return value;
}

template <class leftOperand>
FUNC_INLINE void Ctg<leftOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand>
FUNC_INLINE void Ctg<leftOperand>::diff()
{
	leftOp.accum += (1/(sin(leftOp.value)*sin(leftOp.value)))*accum;
	leftOp.diff();
};

template <class leftOperand>
Ctg<leftOperand> ctg(leftOperand leftOp)
{
	return Ctg<leftOperand>(leftOp);
}

// Exp

template <class leftOperand>
class Exp
{
public:
	Exp();
	Exp(leftOperand arg1);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	leftOperand leftOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand>
Exp<leftOperand>::Exp()
{
	value = 0;
	accum = 0;
}

template <class leftOperand>
Exp<leftOperand>::Exp(leftOperand arg1)
{
	leftOp = arg1;
	accum = 0;
}

template <class leftOperand>
FUNC_INLINE DATA_TYPE Exp<leftOperand>::evaluate()
{
	value = exp(leftOp.value);
	return value;
}

template <class leftOperand>
FUNC_INLINE void Exp<leftOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand>
FUNC_INLINE void Exp<leftOperand>::diff()
{
	leftOp.accum += value*accum;
	leftOp.diff();
};

template <class leftOperand>
Exp<leftOperand> exp(leftOperand leftOp)
{
	return Exp<leftOperand>(leftOp);
}

// Log

template <class leftOperand>
class Log
{
public:
	Log();
	Log(leftOperand arg1);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	leftOperand leftOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand>
Log<leftOperand>::Log()
{
	value = 0;
	accum = 0;
}

template <class leftOperand>
Log<leftOperand>::Log(leftOperand arg1)
{
	leftOp = arg1;
	accum = 0;
}

template <class leftOperand>
FUNC_INLINE DATA_TYPE Log<leftOperand>::evaluate()
{
	value = log(leftOp.value);
	return value;
}

template <class leftOperand>
FUNC_INLINE void Log<leftOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand>
FUNC_INLINE void Log<leftOperand>::diff()
{
	leftOp.accum += accum/leftOp.value;
	leftOp.diff();
};

template <class leftOperand>
Log<leftOperand> log(leftOperand leftOp)
{
	return Log<leftOperand>(leftOp);
}

// Sqrt

template <class leftOperand>
class Sqrt
{
public:
	Sqrt();
	Sqrt(leftOperand arg1);
	FUNC_INLINE DATA_TYPE evaluate();
	FUNC_INLINE void setDiff();
	FUNC_INLINE void diff();

	leftOperand leftOp;
	DATA_TYPE value;
	DATA_TYPE accum;
};

template <class leftOperand>
Sqrt<leftOperand>::Sqrt()
{
	value = 0;
	accum = 0;
}

template <class leftOperand>
Sqrt<leftOperand>::Sqrt(leftOperand arg1)
{
	leftOp = arg1;
	accum = 0;
}

template <class leftOperand>
FUNC_INLINE DATA_TYPE Sqrt<leftOperand>::evaluate()
{
	if(leftOp.value!=0)
	{
		value = sqrt(leftOp.value);
	}
	else
		value = 0;

	return value;
}

template <class leftOperand>
FUNC_INLINE void Sqrt<leftOperand>::setDiff()
{
	accum = 1;
}

template <class leftOperand>
FUNC_INLINE void Sqrt<leftOperand>::diff()
{
	leftOp.accum += (0.5/value)*accum;
	leftOp.diff();
};

template <class leftOperand>
Sqrt<leftOperand> sqrt(leftOperand leftOp)
{
	return Sqrt<leftOperand>(leftOp);
}
#endif