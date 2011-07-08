#ifndef _vec_h_
#define _vec_h_
#include <cmath>

/**
 * This is a three dimensional vector that supports
 * some common vector operations.
 */
class vec
{
	public:
		vec() { u[0] = u[1] = u[2] = 0.0; }
		vec(double x, double y = 0.0, double z = 0.0)
		{
			u[0] = x;
			u[1] = y;
			u[2] = z;
		}
		vec(const vec& src)
		{
			u[0] = src[0];
			u[1] = src[1];
			u[2] = src[2];
		}
		vec operator=(const vec& src)
		{
			u[0] = src[0];
			u[1] = src[1];
			u[2] = src[2];
			return *this;
		}
		const double& operator[](int i) const { return  u[i]; }
		double& operator[](int i) { return  u[i]; }
		vec operator*(double a) const
		{
			vec result;
			result[0] = a*u[0];
			result[1] = a*u[1];
			result[2] = a*u[2];
			return result;
		}
		vec operator/(double a) const
		{
			vec result;
			result[0] = u[0]/a;
			result[1] = u[1]/a;
			result[2] = u[2]/a;
			return result;
		}
		vec operator+(const vec& other) const
		{
			vec result;
			result[0] = other[0]+u[0];
			result[1] = other[1]+u[1];
			result[2] = other[2]+u[2];
			return result;
		}
		vec operator-(const vec& other) const
		{
			vec result;
			result[0] = u[0]-other[0];
			result[1] = u[1]-other[1];
			result[2] = u[2]-other[2];
			return result;
		}
		vec operator+=(const vec& other)
		{
			*this = (*this)+other;
			return *this;
		}
		vec operator-=(const vec& other)
		{
			*this = (*this)-other;
			return *this;
		}
		vec operator/=(double other)
		{
			*this = (*this)/other;
			return *this;
		}
		vec operator*=(double other)
		{
			*this = (*this)*other;
			return *this;
		}
		double norm2() const
		{
			return sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
		}
		friend vec operator*(double a, const vec& b)
		{
			return b*a;
		}
		~vec(){}
	private:
		double u[3];
};

#endif
