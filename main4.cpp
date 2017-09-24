#include <iostream>

using std::integral_constant;

template<size_t N>
struct fibonacci : integral_constant<size_t, fibonacci<N - 1>::value + fibonacci<N - 2>::value > {};

template<1> struct fibonacci<1> : integral_constant<size_t, 1> {};
template<0> struct fibonacci<0> : integral_constant<size_t, 0> {};

int main()
{
	std::cout << fibonacci<10>::value << std::endl;
	return 0;
}
