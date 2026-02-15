#include <iostream>
#include <vector>
#include "complex.cpp"
using namespace std;

static const complex w1(-0.5, 0.8660254037844386467637231707529361834714026269051903140279034897259665084544);
static const complex w2(-0.5, -0.8660254037844386467637231707529361834714026269051903140279034897259665084544);

static size_t find_max(double arr[], size_t size) {
	size_t result = 0;
	double val = arr[0];
	for (size_t j = 0; j < size; j++) {
		if (arr[j] > val) {
			val = arr[j];
			result = j;
		}
	}
	return result;
}

static size_t find_min(double arr[], size_t size) {
	size_t result = 0;
	double val = arr[0];
	for (size_t j = 0; j < size; j++) {
		if (arr[j] < val) {
			val = arr[j];
			result = j;
		}
	}
	return result;
}

static complex *solve_quadratic(const complex &a, const complex &b, const complex &c) {
	complex *x = new complex[2];
	const complex Delta = b * b - 4 * a * c;
	x[0] = (-b + sqrt(Delta)) / (2 * a);
	x[1] = (-b - sqrt(Delta)) / (2 * a);
	return x;
}

static complex *solve_cubic(const complex &a, const complex &b, const complex &c, const complex &d) {
	complex *x = new complex[3];
	const complex tmp = b * (2 * b * b - 9 * a * c);
	const complex p = tmp + 27 * a * a * d;
	const complex q = b * b - 3 * a * c;
	const complex Delta = c * c * (b * b - 4 * a * c) - 2 * tmp * d - 27 * pow(a * d, 2);
	const complex A1 = cbrt((-p + 3 * a * sqrt(-3 * Delta)) / 2);
	const complex A2 = (A1 == 0) ? cbrt((-p - 3 * a * sqrt(-3 * Delta)) / 2) : q / A1;
	x[0] = (-b + A1 + A2) / (3 * a);
	x[1] = (-b + w1 * A1 + w2 * A2) / (3 * a);
	x[2] = (-b + w2 * A1 + w1 * A2) / (3 * a);
	if (a.im() == 0 && b.im() == 0 && c.im() == 0 && d.im() == 0) {
		if (Delta.re() >= 0) {
			x[0] = x[0].re();
			x[1] = x[1].re();
			x[2] = x[2].re();
		} else {
			double xs[] = {abs(x[0].im()), abs(x[1].im()), abs(x[2].im())};
			size_t min = find_min(xs, 3);
			x[min] = x[min].re();
			if (a / b == c / d) {
				x[(min + 1) % 3] = x[(min + 1) % 3].im() * i;
				x[(min + 2) % 3] = x[(min + 1) % 3].conj();
			}
		}
	}
	if (d == 0) {
		double xt[] = { abs(x[0].re() + x[0].im()), abs(x[1].re() + x[1].im()), abs(x[2].re() + x[2].im())};
		size_t min = find_min(xt, 3);
		size_t max = find_max(xt, 3);
		if (c == 0) {
			x[(max + 1) % 3] = 0;
			x[(max + 2) % 3] = 0;
		}
		else
			x[min] = 0;
	}
	return x;
}

static complex *solve_quartic(const complex &a, const complex &b, const complex &c, const complex &d, const complex &e) {
	complex *x = new complex[4];
	const complex *ks = solve_cubic(1, -c, b * d - 4 * a * e, -b * b * e - a * d * d + 4 * a * c * e);
	complex k = ks[0];
	if (ks[0].im() != 0) {
		if (ks[1].im() == 0)
			k = ks[1];
		else if (ks[2].im() == 0)
			k = ks[2];
	}
	delete[] ks;
	const complex T = sqrt(b * b - 4 * a * (c - k));
	const complex U = (T == 0) ? sqrt(pow(b, 4) - 4 * a * b * b * (c + k) + 16 * a * a * (b * d - 4 * a * e + k * k)) : (b * b * b - 4 * a * b * c + 8 * a * a * d) / T;
	const complex tmp = 2 * b * b - 4 * a * (c + k);
	x[0] = (-b + T + sqrt(tmp - 2 * U)) / (4 * a);
	x[1] = (-b + T - sqrt(tmp - 2 * U)) / (4 * a);
	x[2] = (-b - T + sqrt(tmp + 2 * U)) / (4 * a);
	x[3] = (-b - T - sqrt(tmp + 2 * U)) / (4 * a);
	return x;
}

static int solve(const int& degree, const vector<complex>& coeff) {
	size_t n = coeff.size();
	complex* x = nullptr;
	switch (degree) {
	case 2:
		if (n != 3) {
			cerr << "二次方程必须有三个系数。";
			return 1;
		}
		if (coeff[0] == 0) {
			cerr << "最高次系数不能为0。";
			return 1;
		}
		x = solve_quadratic(coeff[0], coeff[1], coeff[2]);
		break;
	case 3:
		if (n != 4) {
			cerr << "三次方程必须有四个系数。";
			return 1;
		}
		if (coeff[0] == 0) {
			cerr << "最高次系数不能为0。";
			return 1;
		}
		x = solve_cubic(coeff[0], coeff[1], coeff[2], coeff[3]);
		break;
	case 4:
		if (n != 5) {
			cerr << "四次方程必须有五个系数。";
			return 1;
		}
		if (coeff[0] == 0) {
			cerr << "最高次系数不能为0。";
			return 1;
		}
		x = solve_quartic(coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
		break;
	default:
		cerr << "无效方程次数";
		return 1;
	}
	cout << endl;
	for (int j = 0; j < degree - 1; j++) {
		cout << "x" << j + 1 << "=" << x[j] << endl;
	}
	cout << "x" << degree << "=" << x[degree - 1] << endl;
	delete[] x;
	return 0;
}

static void help(const string& prog) {
	cout << "使用方法: " << endl;
	cout << "  二次方程: " << prog << " -qd <a> <b> <c>            //a*x^2+b*x+c=0" << endl;
	cout << "  三次方程: " << prog << " -cb <a> <b> <c> <d>        //a*x^3+b*x^2+c*x+d=0" << endl;
	cout << "  四次方程: " << prog << " -qt <a> <b> <c> <d> <e>    //a*x^4+b*x^3+c*x^2+d*x+e=0" << endl;
}

int main(size_t argc, char *argv[]) {
	int degree;
	vector<complex> coeff;
	string mode;
	if (argc > 1) {
		mode = argv[1];
		if (mode == "-qd")
			degree = 2;
		else if (mode == "-cb")
			degree = 3;
		else if (mode == "-qt")
			degree = 4;
		else if (mode == "-h") {
			help(argv[0]);
			return 0;
		} else {
			cerr << "无效方程类型";
			return 1;
		}
		for (size_t j = 2; j < argc; j++)
			coeff.push_back(stoc(argv[j]));
		return solve(degree, coeff);
	}
	while (true) {
		cout << "请输入方程次数(2-4)，或输入'q'退出: ";
		cin >> mode;
		if (mode == "q")
			break;
		try { degree = stoi(mode); }
		catch (const invalid_argument& e) { degree = 0; }
		catch (const out_of_range& e) { degree = 0; }
		switch (degree) {
			case 2:
				coeff.resize(3);
				cout << "a*x^2+b*x+c=0" << endl;
				cout << "a=";
				cin >> coeff[0];
				cout << "b=";
				cin >> coeff[1];
				cout << "c=";
				cin >> coeff[2];
				break;
			case 3:
				coeff.resize(4);
				cout << "a*x^3+b*x^2+c*x+d=0" << endl;
				cout << "a=";
				cin >> coeff[0];
				cout << "b=";
				cin >> coeff[1];
				cout << "c=";
				cin >> coeff[2];
				cout << "d=";
				cin >> coeff[3];
				break;
			case 4:
				coeff.resize(5);
				cout << "a*x^4+b*x^3+c*x^2+d*x+e=0" << endl;
				cout << "a=";
				cin >> coeff[0];
				cout << "b=";
				cin >> coeff[1];
				cout << "c=";
				cin >> coeff[2];
				cout << "d=";
				cin >> coeff[3];
				cout << "e=";
				cin >> coeff[4];
				break;
		}
		solve(degree, coeff);
		cout << endl;
	}
	return 0;
}