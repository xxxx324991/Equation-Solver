#include <cmath>
#include <sstream>
constexpr auto PI = 3.1415926535897932384626433832795;
constexpr auto PI_2 = 1.57079632679489661923132169163975;
constexpr auto PI_4 = 0.78539816339744830961566084581988;

class complex {
	private:
		double real, imag;
	public:
		complex(const double &re = 0, const double &im = 0);
		double re() const;
		double im() const;
		complex conj() const;
		double theta() const;

		friend double abs(const complex &z);
		friend complex operator+(const complex &z1, const complex &z2);
		friend complex operator+(const double &n, const complex &z);
		friend complex operator-(const complex &z1, const complex &z2);
		friend complex operator-(const double &n, const complex &z);
		friend complex operator-(const complex &z);
		friend complex operator*(const double &n, const complex &z);
		friend complex operator*(const complex &z1, const complex &z2);
		friend complex operator/(const double &n, const complex &z);
		friend complex operator/(const complex &z1, const complex &z2);
		friend complex exp(const complex &z);
		friend complex pow(const complex &z, const double &n);
		friend complex pow(const complex &z, const complex &n);
		friend complex sqrt(const complex &z);
		friend complex root(const complex &z, const double &n);
		friend complex log(const complex &z);
		friend complex sin(const complex &z);
		friend complex cos(const complex &z);
		friend bool operator==(const complex &z1, const complex &z2);
		friend bool operator!=(const complex &z1, const complex &z2);
		friend std::string to_string(const complex &z);
		friend std::ostream &operator<<(std::ostream &out, const complex &z);
};