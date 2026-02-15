#include "complex.h"

static const complex i(0, 1);

complex::complex(const double& re, const double& im) {
	real = re;
	imag = im;
}

double complex::re() const { return real; }

double complex::im() const { return imag; }

complex complex::conj() const { return complex(real, -imag); }

double complex::theta() const {
	if (imag == 0) {
		if (real >= 0)
			return 0;
		return PI;
	}
	if (real == 0) {
		if (imag > 0)
			return PI_2;
		return -PI_2;
	}
	if (real == imag) {
		if (real > 0)
			return PI_4;
		return -PI_2 - PI_4;
	}
	if (real == -imag) {
		if (real > 0)
			return -PI_4;
		return PI_2 + PI_4;
	}
	const double theta = acos(real / sqrt(real * real + imag * imag));
	if (imag < 0)
		return -theta;
	return theta;
}

double abs(const complex& z) {
	if (z.real == 0 || z.imag == 0)
		return abs(z.real + z.imag);
	return sqrt(z.real * z.real + z.imag * z.imag);
}

complex operator+(const complex& z1, const complex& z2) { return complex(z1.real + z2.real, z1.imag + z2.imag); }

complex operator+(const double& n, const complex& z) { return complex(n + z.real, z.imag); }

complex operator-(const complex& z1, const complex& z2) { return complex(z1.real - z2.real, z1.imag - z2.imag); }

complex operator-(const double& n, const complex& z) { return complex(n - z.real, -z.imag); }

complex operator-(const complex& z) { return complex(-z.real, -z.imag); }

complex operator*(const double& n, const complex& z) { return complex(n * z.real, n * z.imag); }

complex operator*(const complex& z1, const complex& z2) { return complex(z1.real * z2.real - z1.imag * z2.imag, z1.imag * z2.real + z1.real * z2.imag); }

complex operator/(const double& n, const complex& z) {
	const double g = z.real * z.real + z.imag * z.imag;
	return complex(n * z.real / g, -n * z.imag / g);
}

complex operator/(const complex& z1, const complex& z2) {
	const double g = z2.real * z2.real + z2.imag * z2.imag;
	return complex((z1.real * z2.real + z1.imag * z2.imag) / g, (z1.imag * z2.real - z1.real * z2.imag) / g);
}

complex exp(const complex& z) { return exp(z.real) * complex(cos(z.imag), sin(z.imag)); }

static complex pow(const complex& z, const int& n) {
	if (n == 0 || z == 1)
		return 1;
	if (n == 1)
		return z;
	if (n < 0)
		return 1 / pow(z, -n);
	const complex tmp = pow(z, n / 2);
	if (n % 2 == 0)
		return tmp * tmp;
	return  z * tmp * tmp;
}

complex pow(const complex& z, const double& n) {
	const double t = n * z.theta();
	return pow(z.real * z.real + z.imag * z.imag, n / 2) * complex(cos(t), sin(t));
}

complex pow(const complex& z, const complex& n) {
	const double j = exp(-n.imag * z.theta());
	const double k = n.imag * log(z.real * z.real + z.imag * z.imag) / 2;
	return pow(z, n.real) * complex(j * cos(k), j * sin(k));
}

complex sqrt(const complex& z) {
	const double m = abs(z);
	if (z.imag < 0)
		return complex(sqrt((m + z.real) / 2), -sqrt((m - z.real) / 2));
	return complex(sqrt((m + z.real) / 2), sqrt((m - z.real) / 2));
}

static complex cbrt(const complex& z) {
	const double t = z.theta() / 3;
	return cbrt(abs(z)) * complex(cos(t), sin(t));
}

complex root(const complex& z, const double& n) {
	const double t = z.theta() / n;
	return pow(z.real * z.real + z.imag * z.imag, 0.5 / n) * complex(cos(t), sin(t));
}

static complex root(const complex& z, const complex& n) { return pow(z, 1 / n); }

complex log(const complex& z) { return complex(log(z.real * z.real + z.imag + z.imag) / 2, z.theta()); }

complex sin(const complex& z) { return complex(sin(z.real) * cosh(z.imag), cos(z.real) * sinh(z.imag)); }

complex cos(const complex& z) { return complex(cos(z.real) * cosh(z.imag), -sin(z.real) * sinh(z.imag)); }

static complex tan(const complex& z) { return sin(z) / cos(z); }

static complex asin(const complex& z) { return complex(0, -1) * log(i * z + sqrt(1 - z * z)); }

static complex acos(const complex& z) { return PI_2 + i * log(i * z + sqrt(1 - z * z)); }

static complex atan(const complex& z) { return i * log((1 - i * z) / (1 + i * z)) / 2; }

static complex sinh(const complex& z) { return (exp(z) - exp(-z)) / 2; }

static complex cosh(const complex& z) { return (exp(z) + exp(-z)) / 2; }

static complex tanh(const complex& z) {
	const complex A = exp(z);
	const complex B = exp(-z);
	return (A - B) / (A + B);
}

static complex asinh(const complex& z) { return log(z + sqrt(z * z + 1)); }

static complex acosh(const complex& z) { return log(z + sqrt(z * z - 1)); }

static complex atanh(const complex& z) { return log((1 + z) / (1 - z)) / 2; }

bool operator==(const complex& z1, const complex& z2) { return z1.real == z2.real && z1.imag == z2.imag; }

bool operator!=(const complex& z1, const complex& z2) { return z1.real != z2.real || z1.imag != z2.imag; }

static void operator+=(complex& z1, const complex& z2) { z1 = z1 + z2; }

static void operator-=(complex& z1, const complex& z2) { z1 = z1 - z2; }

static void operator*=(complex& z1, const complex& z2) { z1 = z1 * z2; }

static void operator/=(complex& z1, const complex& z2) { z1 = z1 / z2; }

static complex stoc(const std::string& str) {
	double re = 0;
	double im = 0;
	std::string tmp = str;
	bool neg = false;
	if (str.at(0) == '-') {
		neg = true;
		tmp = tmp.substr(1, std::string::npos);
	}
	if (tmp.back() == 'i') {
		if (str.size() == 1)
			im = 1;
		else {
			size_t plusPos = tmp.find('+');
			size_t minPos = tmp.find('-');
			if (plusPos != std::string::npos) {
				re = stod(tmp.substr(0, plusPos));
				std::string h = tmp.substr(plusPos + 1, std::string::npos);
				if (h.size() == 1)
					im = 1;
				else
					im = stod(h);
			}
			else if (minPos != std::string::npos) {
				re = stod(tmp.substr(0, minPos));
				std::string h = tmp.substr(minPos + 1, std::string::npos);
				if (h.size() == 1)
					im = -1;
				else
					im = -stod(h);
			}
			else
				im = stod(tmp);
		}
	}
	else
		re = stod(tmp);
	if (neg) {
		if (re == 0)
			im = -im;
		else
			re = -re;
	}
	return complex(re, im);
}

std::string to_string(const complex& z) {
	if (z.imag < 0)
		return to_string(z.real) + to_string(z.imag) + "i";
	return to_string(z.real) + "+" + to_string(z.imag) + "i";
}

static std::istream& operator>>(std::istream& in, complex& z) {
	std::string str;
	in >> str;
	z = stoc(str);
	return in;
}

std::ostream& operator<<(std::ostream& out, const complex& z) {
	std::stringstream ss_real;
	std::stringstream ss_imag;
	ss_real << z.real;
	ss_imag << z.imag;
	std::string realStr = ss_real.str();
	std::string imagStr = ss_imag.str();
	if (imagStr != "0" && imagStr != "-0") {
		if (realStr != "0" && realStr != "-0") {
			out << realStr;
			if (imagStr.at(0) == '-') {
				if (imagStr == "-1")
					out << "-";
				else
					out << imagStr;
			}
			else {
				out << "+";
				if (imagStr != "1")
					out << imagStr;
			}
		}
		else {
			if (imagStr == "-1")
				out << "-";
			else if (imagStr != "1")
				out << imagStr;
		}
		out << "i";
	}
	else{
		if (realStr == "-0")
			realStr = "0";
		out << realStr;
	}
	return out;
}