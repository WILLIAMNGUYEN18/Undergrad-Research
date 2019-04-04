#include "quatnorm.h"

const double PI = 3.14159265358979;

QuatNorm::QuatNorm() {
	x = y = z = 0.0;
	w = 1.0;
}

QuatNorm::QuatNorm( double x, double y, double z, double w ) {
	this->x = x;
	this->y = y;
	this->z = z;
	this->w = w;
}

QuatNorm::QuatNorm( Vec3d v, double w ) {
	x = v[0];
	y = v[1];
	z = v[2];
	this->w = w;
}

QuatNorm::QuatNorm( double theta, Vec3d v ) {
	double s = sin(theta/2.0);
	x = v[0] * s;
	y = v[1] * s;
	z = v[2] * s;
	this->w = cos(theta/2.0);
}

QuatNorm::QuatNorm( double thetaX, double thetaY, double thetaZ ) {
	thetaX *= 0.5;
	thetaY *= 0.5;
	thetaZ *= 0.5;

	Vec3d		xRot = Vec3d(1,0,0) * sin( thetaX );
	Vec3d		yRot = Vec3d(0,1,0) * sin( thetaY );
	Vec3d		zRot = Vec3d(0,0,1) * sin( thetaZ );
	QuatNorm	xQuat( xRot, cos(thetaX) );
	QuatNorm	yQuat( yRot, cos(thetaY) );
	QuatNorm	zQuat( zRot, cos(thetaZ) );
	QuatNorm	composite = zQuat * yQuat * xQuat;

	composite.normalize();

	*this = composite;
}

void QuatNorm::getMatrices( Mat4d& m, Mat4d* deriv ) const {
	double t1 = x*x;
	double t2 = y*y;
	double t3 = z*z;
	double t4 = w*w;
	double t5 = t1-t2-t3+t4;
	double t6 = t1+t2+t3+t4;
	double t7 = 1/t6;
	double t9 = x*y;
	double t10 = w*z;
	double t11 = t9+t10;
	double t14 = x*z;
	double t15 = w*y;
	double t16 = t14-t15;
	double t19 = t9-t10;
	double t22 = t1-t2+t3-t4;
	double t24 = w*x;
	double t25 = y*z;
	double t26 = t24+t25;
	double t29 = t15+t14;
	double t32 = -t25+t24;
	double t35 = t1+t2-t3-t4;
	double t37 = x*t7;
	double t38 = t6*t6;
	double t39 = 1/t38;
	double t40 = t5*t39;
	double t43 = y*t7;
	double t44 = 2.0*t43;
	double t45 = t11*t39;
	double t49 = z*t7;
	double t50 = 2.0*t49;
	double t51 = t16*t39;
	double t55 = t19*t39;
	double t59 = t22*t39;
	double t62 = w*t7;
	double t63 = 2.0*t62;
	double t64 = t26*t39;
	double t68 = t29*t39;
	double t72 = t32*t39;
	double t76 = t35*t39;
	double t81 = 2.0*t37;

	m[0][0] = t5*t7;
	m[0][1] = 2.0*t11*t7;
	m[0][2] = 2.0*t16*t7;
	m[0][3] = 0.0;
	m[1][0] = 2.0*t19*t7;
	m[1][1] = -t22*t7;
	m[1][2] = 2.0*t26*t7;
	m[1][3] = 0.0;
	m[2][0] = 2.0*t29*t7;
	m[2][1] = -2.0*t32*t7;
	m[2][2] = -t35*t7;
	m[2][3] = 0.0;
	m[3][0] = 0.0;
	m[3][1] = 0.0;
	m[3][2] = 0.0;
	m[3][3] = 1.0;
	deriv[0][0][0] = 2.0*t37-2.0*t40*x;
	deriv[0][0][1] = t44-4.0*t45*x;
	deriv[0][0][2] = t50-4.0*t51*x;
	deriv[0][0][3] = 0.0;
	deriv[0][1][0] = t44-4.0*t55*x;
	deriv[0][1][1] = -2.0*t37+2.0*t59*x;
	deriv[0][1][2] = t63-4.0*t64*x;
	deriv[0][1][3] = 0.0;
	deriv[0][2][0] = t50-4.0*t68*x;
	deriv[0][2][1] = -t63+4.0*t72*x;
	deriv[0][2][2] = -2.0*t37+2.0*t76*x;
	deriv[0][2][3] = 0.0;
	deriv[0][3][0] = 0.0;
	deriv[0][3][1] = 0.0;
	deriv[0][3][2] = 0.0;
	deriv[0][3][3] = 0.0;
	deriv[1][0][0] = -2.0*t43-2.0*t40*y;
	deriv[1][0][1] = t81-4.0*t45*y;
	deriv[1][0][2] = -t63-4.0*t51*y;
	deriv[1][0][3] = 0.0;
	deriv[1][1][0] = t81-4.0*t55*y;
	deriv[1][1][1] = 2.0*t43+2.0*t59*y;
	deriv[1][1][2] = t50-4.0*t64*y;
	deriv[1][1][3] = 0.0;
	deriv[1][2][0] = t63-4.0*t68*y;
	deriv[1][2][1] = t50+4.0*t72*y;
	deriv[1][2][2] = -2.0*t43+2.0*t76*y;
	deriv[1][2][3] = 0.0;
	deriv[1][3][0] = 0.0;
	deriv[1][3][1] = 0.0;
	deriv[1][3][2] = 0.0;
	deriv[1][3][3] = 0.0;
	deriv[2][0][0] = -2.0*t49-2.0*t40*z;
	deriv[2][0][1] = t63-4.0*t45*z;
	deriv[2][0][2] = t81-4.0*t51*z;
	deriv[2][0][3] = 0.0;
	deriv[2][1][0] = -t63-4.0*t55*z;
	deriv[2][1][1] = -2.0*t49+2.0*t59*z;
	deriv[2][1][2] = t44-4.0*t64*z;
	deriv[2][1][3] = 0.0;
	deriv[2][2][0] = t81-4.0*t68*z;
	deriv[2][2][1] = t44+4.0*t72*z;
	deriv[2][2][2] = 2.0*t49+2.0*t76*z;
	deriv[2][2][3] = 0.0;
	deriv[2][3][0] = 0.0;
	deriv[2][3][1] = 0.0;
	deriv[2][3][2] = 0.0;
	deriv[2][3][3] = 0.0;
	deriv[3][0][0] = 2.0*t62-2.0*t40*w;
	deriv[3][0][1] = t50-4.0*t45*w;
	deriv[3][0][2] = -t44-4.0*t51*w;
	deriv[3][0][3] = 0.0;
	deriv[3][1][0] = -t50-4.0*t55*w;
	deriv[3][1][1] = 2.0*t62+2.0*t59*w;
	deriv[3][1][2] = t81-4.0*t64*w;
	deriv[3][1][3] = 0.0;
	deriv[3][2][0] = t44-4.0*t68*w;
	deriv[3][2][1] = -t81+4.0*t72*w;
	deriv[3][2][2] = 2.0*t62+2.0*t76*w;
	deriv[3][2][3] = 0.0;
	deriv[3][3][0] = 0.0;
	deriv[3][3][1] = 0.0;
	deriv[3][3][2] = 0.0;
	deriv[3][3][3] = 0.0;
}

Mat4f QuatNorm::toMatrix() {
	Mat4f m;
	double t1 = x*x;
	double t2 = y*y;
	double t3 = z*z;
	double t4 = w*w;
	double t5 = t1-t2-t3+t4;
	double t6 = t1+t2+t3+t4;
	double t7 = 1/t6;
	double t9 = x*y;
	double t10 = w*z;
	double t11 = t9+t10;
	double t14 = x*z;
	double t15 = w*y;
	double t16 = t14-t15;
	double t19 = t9-t10;
	double t22 = t1-t2+t3-t4;
	double t24 = w*x;
	double t25 = y*z;
	double t26 = t24+t25;
	double t29 = t15+t14;
	double t32 = -t25+t24;
	double t35 = t1+t2-t3-t4;
	double t37 = x*t7;
	double t38 = t6*t6;
	double t39 = 1/t38;
	double t40 = t5*t39;
	double t43 = y*t7;
	double t44 = 2.0*t43;
	double t45 = t11*t39;
	double t49 = z*t7;
	double t50 = 2.0*t49;
	double t51 = t16*t39;
	double t55 = t19*t39;
	double t59 = t22*t39;
	double t62 = w*t7;
	double t63 = 2.0*t62;
	double t64 = t26*t39;
	double t68 = t29*t39;
	double t72 = t32*t39;
	double t76 = t35*t39;
	double t81 = 2.0*t37;

	m[0][0] = t5*t7;
	m[0][1] = 2.0*t11*t7;
	m[0][2] = 2.0*t16*t7;
	m[0][3] = 0.0;
	m[1][0] = 2.0*t19*t7;
	m[1][1] = -t22*t7;
	m[1][2] = 2.0*t26*t7;
	m[1][3] = 0.0;
	m[2][0] = 2.0*t29*t7;
	m[2][1] = -2.0*t32*t7;
	m[2][2] = -t35*t7;
	m[2][3] = 0.0;
	m[3][0] = 0.0;
	m[3][1] = 0.0;
	m[3][2] = 0.0;
	m[3][3] = 1.0;
	return m;
}

Mat4d QuatNorm::toMatrixD() {
	Mat4d m;
	double t1 = x*x;
	double t2 = y*y;
	double t3 = z*z;
	double t4 = w*w;
	double t5 = t1-t2-t3+t4;
	double t6 = t1+t2+t3+t4;
	double t7 = 1/t6;
	double t9 = x*y;
	double t10 = w*z;
	double t11 = t9+t10;
	double t14 = x*z;
	double t15 = w*y;
	double t16 = t14-t15;
	double t19 = t9-t10;
	double t22 = t1-t2+t3-t4;
	double t24 = w*x;
	double t25 = y*z;
	double t26 = t24+t25;
	double t29 = t15+t14;
	double t32 = -t25+t24;
	double t35 = t1+t2-t3-t4;
	double t37 = x*t7;
	double t38 = t6*t6;
	double t39 = 1/t38;
	double t40 = t5*t39;
	double t43 = y*t7;
	double t44 = 2.0*t43;
	double t45 = t11*t39;
	double t49 = z*t7;
	double t50 = 2.0*t49;
	double t51 = t16*t39;
	double t55 = t19*t39;
	double t59 = t22*t39;
	double t62 = w*t7;
	double t63 = 2.0*t62;
	double t64 = t26*t39;
	double t68 = t29*t39;
	double t72 = t32*t39;
	double t76 = t35*t39;
	double t81 = 2.0*t37;

	m[0][0] = t5*t7;
	m[0][1] = 2.0*t11*t7;
	m[0][2] = 2.0*t16*t7;
	m[0][3] = 0.0;
	m[1][0] = 2.0*t19*t7;
	m[1][1] = -t22*t7;
	m[1][2] = 2.0*t26*t7;
	m[1][3] = 0.0;
	m[2][0] = 2.0*t29*t7;
	m[2][1] = -2.0*t32*t7;
	m[2][2] = -t35*t7;
	m[2][3] = 0.0;
	m[3][0] = 0.0;
	m[3][1] = 0.0;
	m[3][2] = 0.0;
	m[3][3] = 1.0;
	return m;
}

void QuatNorm::normalize() {
	double		magnitude = sqrt( x*x + y*y + z*z + w*w );

	if( magnitude == 0.0 )
		return;

	x /= magnitude;
	y /= magnitude;
	z /= magnitude;
	w /= magnitude;
}

void QuatNorm::conjugate() {
	x = -x;
	y = -y;
	z = -z;
}

ostream& operator <<( ostream& os, const QuatNorm & q )
{
	return os << q.x << " " << q.y << " " << q.z << " " << q.w;
}

istream& operator >>( istream& is, QuatNorm & q )
{
	return is >> q.x >> q.y >> q.z >> q.w;
}

double quatDist(QuatNorm a, QuatNorm b) {
	a.normalize();
	b.normalize();
	
	double d = fabs(a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w);
	if (d > 1)
		d = 1;
	d = acos(d);
	return d;
}

double quatDist(const double *a, const double *b) {
	// assume normalized
	double d = fabs(a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]);
	if (d > 1)
		d = 1;
	d = acos(d);
	return d;
}

QuatNorm slerp(QuatNorm from, QuatNorm to, double t) {
	from.normalize();
	to.normalize();

	// from GamaSutra:
	//  http://www.gamasutra.com/features/19980703/quaternions_01.htm
	QuatNorm res;
	float to1[4];
	double omega, cosom, sinom, scale0, scale1;

	// calc cosine
	cosom = from.x * to.x + from.y * to.y + from.z * to.z
		   + from.w * to.w;

	// adjust signs (if necessary)
	if (cosom < 0.0){ 
		cosom = -cosom;
		to1[0] = -to.x;
		to1[1] = - to.y;
		to1[2] = - to.z;
		to1[3] = - to.w;
	}
	else  {
		to1[0] = to.x;
		to1[1] = to.y;
		to1[2] = to.z;
		to1[3] = to.w;
	}

	// calculate coefficients
	if ( (1.0 - cosom) > 0.0001 ) {
		// standard case (slerp)
		omega = acos(cosom);
		sinom = sin(omega);
		scale0 = sin((1.0 - t) * omega) / sinom;
		scale1 = sin(t * omega) / sinom;
	}
	else {        
		// "from" and "to" quaternions are very close 
		//  ... so we can do a linear interpolation
		scale0 = 1.0 - t;
		scale1 = t;
	}

	// calculate final values
	res.x = scale0 * from.x + scale1 * to1[0];
	res.y = scale0 * from.y + scale1 * to1[1];
	res.z = scale0 * from.z + scale1 * to1[2];
	res.w = scale0 * from.w + scale1 * to1[3];
	return res;
}

void QuatNorm::saveRIB(ostream &out) {
	Vec3d v(x, y, z);
	if (v.iszero())
		v = Vec3d(1, 0, 0);
	else
		v.normalize();
	double angle = -acos(w) * 2.0 * 180.0 / PI;
	out << "Rotate " << angle << " " << v[0] << " " << v[1] << " " << v[2] << endl;
}