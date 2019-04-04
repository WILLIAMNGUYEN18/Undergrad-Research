#ifndef __QUAT_NORM_HEADER__
#define __QUAT_NORM_HEADER__

#include "vec.h"
#include "mat.h"

//==========[ class QuatNorm ]=================================================

class QuatNorm {
public:

	//---[ public variables ]------------------------------

	double		x, y, z, w;

	//---[ constructors ]----------------------------------

	QuatNorm();
	QuatNorm( double x, double y, double z, double w );
	QuatNorm( Vec3d v, double w );
	QuatNorm( double theta, Vec3d v );
	QuatNorm( double thetaX, double thetaY, double thetaZ );

	//---[ conversion operator ]---------------------------

	void getMatrices( Mat4d& mat, Mat4d* deriv ) const;
	Mat4f toMatrix();
	Mat4d toMatrixD();

	//---[ normalize/conjugate ]---------------------------

	void normalize();
	void conjugate();

	//---[ overloaded operators ]--------------------------

	QuatNorm& operator+=( const QuatNorm& a );
	double &operator[]( int i );

	friend QuatNorm operator+( const QuatNorm& a, const QuatNorm& b );
	friend QuatNorm operator*( const QuatNorm& a, const QuatNorm& b );
	friend QuatNorm operator*( double s, const QuatNorm& a );

	//---[ I/O Operators ]-----------------------

	friend ostream& operator <<(ostream& s, const QuatNorm& q);
	friend istream& operator >>(istream& s, QuatNorm& v);
	void saveRIB(ostream &out);
};

inline QuatNorm& QuatNorm::operator+=( const QuatNorm& a ) {
	x += a.x;
	y += a.y;
	z += a.z;
	w += a.w;

	return *this;
}

inline double &QuatNorm::operator[]( int i ) {
	switch (i) {
	case 0:
		return x;
	case 1:
		return y;
	case 2:
		return z;
	default:
		return w;
	}
}


inline QuatNorm operator+( const QuatNorm& a, const QuatNorm& b ) {
	QuatNorm		q;

	q.x = a.x + b.x;
	q.y = a.y + b.y;
	q.z = a.z + b.z;
	q.w = a.w + b.w;

	return q;
}

inline QuatNorm operator*( double s, const QuatNorm& a ) {
	QuatNorm	qq;

	qq.w = a.w * s;
	qq.x = a.x * s;
	qq.y = a.y * s;
	qq.z = a.z * s;

	return qq;
}

inline QuatNorm operator*( const QuatNorm& a, const QuatNorm& b ) {
	QuatNorm		qq;

	qq.w = (a.w*b.w) - (a.x*b.x) - (a.y*b.y) - (a.z*b.z);
	qq.x = (a.w*b.x) + (a.x*b.w) + (a.y*b.z) - (a.z*b.y);
	qq.y = (a.w*b.y) + (a.y*b.w) + (a.z*b.x) - (a.x*b.z);
	qq.z = (a.w*b.z) + (a.z*b.w) + (a.x*b.y) - (a.y*b.x);
	
	return qq;
}

ostream& operator <<(ostream& s, const QuatNorm& q);
istream& operator >>(istream& s, QuatNorm& v);

double quatDist(QuatNorm a, QuatNorm b);
double quatDist(const double *a, const double *b);
inline double quatDist(const Vec4d &a, const Vec4d &b) {
	return quatDist(a.n, b.n);
}

QuatNorm slerp(QuatNorm q0, QuatNorm q1, double u);

#endif
