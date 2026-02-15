#include <iostream>
#include <cmath>
using namespace std;

static double quick_pow(double a, int b) {

	if (b == 0 || a == 1)
		return 1;
	if (b == 1)
		return a;
	if (b < 0)
		return 1 / quick_pow(a, -b);
	double tmp = quick_pow(a, b / 2);
	if (b % 2 == 0)
		return tmp * tmp;
	return  a * tmp * tmp;
}

static double quick_pow(double a, double b) {
	int p = round(b);
	double q = b - p;
	return pow(a, q) * quick_pow(a, p);
}