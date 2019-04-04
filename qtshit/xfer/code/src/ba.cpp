/*
	ba.cpp

	Various convenient constants and functions.

	Brett Allen
	2002
*/

#include <iostream>
#include <fstream>
using namespace std;
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/glu.h>
#include <time.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/timeb.h>

#include "ba.h"

double zero = 0.0;
double NaN = 0.0 / zero;

//double baPow(double x, double y) {
//	return pow(x, y);
//}

Mat4d invertTransform(Mat4d &m) {
	Mat4d rot = m.transpose();
	Mat4d trans;

	trans[0][3] = -rot[3][0];
	trans[1][3] = -rot[3][1];
	trans[2][3] = -rot[3][2];
	rot[3][0] = 0;
	rot[3][1] = 0;
	rot[3][2] = 0;

	return rot * trans;
}
/*
QuatNorm matToQuat(Mat3d &R) {
	double tmp = fabs(R[0][0] + R[1][1] + R[3][3] + 1);
    double s_ = 0.5*sqrt(tmp);
	Vec3d v_;

      if(s_ > 1e-10)
      {
         v_[0] = (R[2,1]-R[1,2])/(4*s_);
         v_[1] = (R[0,2]-R[2,0])/(4*s_);
         v_[2] = (R[1,0]-R[0,1])/(4*s_);
      }
      else
      {
         // |w| <= 1/2
         static int s_iNext[3] = { 1, 2, 0 };
         int i = 1;
         if ( R[1,1] > R[0,0] )
            i = 2;
         if ( R[2,2] > R[1,1] )
            i = 3;
         int j = s_iNext[i-1];
         int k = s_iNext[j-1];

         double fRoot = sqrt(R[i][i]-R[j][j]-R[k][k] + 1.0);

         double *tmp[3] = { &v_[0], &v_[1], &v_[2] };
         *tmp[i-1] = 0.5*fRoot;
         fRoot = 0.5/fRoot;
         s_ = (R[k][j]-R[j][k])*fRoot;
         *tmp[j-1] = (R[j][i]+R[i][j])*fRoot;
         *tmp[k-1] = (R[k][i]+R[i][k])*fRoot;
      }
	  QuatNorm q(v_[0], v_[1], v_[2], s_);
	  q.normalize();
	  return q;
}
*/
/*
QuatNorm matToQuat(Mat3d &matrix) {
	double x, y, z, w;

	// Find the diagonal of the matrix by adding up its diagonal indices.
	// This is also known as the "trace", but I will call the variable diagonal.
	float diagonal = matrix[0][0] + matrix[1][1] + matrix[2][2] + 1;
	float scale = 0.0f;

	// Below we check if the diagonal is greater than zero.  To avoid accidents with
	// floating point numbers, we substitute 0 with 0.00000001.  If the diagonal is
	// great than zero, we can perform an "instant" calculation, otherwise we will need
	// to identify which diagonal element has the greatest value.  Note, that it appears
	// that %99 of the time, the diagonal IS greater than 0 so the rest is rarely used.

	// If the diagonal is greater than zero
	if(diagonal > 0.00000001)
	{
		// Calculate the scale of the diagonal
		scale = float(sqrt(diagonal ) * 2);

		// Calculate the x, y, x and w of the quaternion through the respective equation
		x = ( matrix[2][1] - matrix[1][2] ) / scale;
		y = ( matrix[0][2] - matrix[2][0] ) / scale;
		z = ( matrix[1][0] - matrix[0][1] ) / scale;
		w = 0.25f * scale;
	}
	else 
	{
		// If the first element of the diagonal is the greatest value
		if ( matrix[0][0] > matrix[1][1] && matrix[0][0] > matrix[2][2] )  
		{	
			// Find the scale according to the first element, and double that value
			scale  = (float)sqrt( 1.0f + matrix[0][0] - matrix[1][1] - matrix[2][2] ) * 2.0f;

			// Calculate the x, y, x and w of the quaternion through the respective equation
			x = 0.25f * scale;
			y = (matrix[1][0] + matrix[0][1] ) / scale;
			z = (matrix[0][2] + matrix[2][0] ) / scale;
			w = (matrix[2][1] - matrix[1][2] ) / scale;	
		} 
		// Else if the second element of the diagonal is the greatest value
		else if ( matrix[1][1] > matrix[2][2] ) 
		{
			// Find the scale according to the second element, and double that value
			scale  = (float)sqrt( 1.0f + matrix[1][1] - matrix[0][0] - matrix[2][2] ) * 2.0f;
			
			// Calculate the x, y, x and w of the quaternion through the respective equation
			x = (matrix[1][0] + matrix[0][1] ) / scale;
			y = 0.25f * scale;
			z = (matrix[2][1] + matrix[1][2] ) / scale;
			w = (matrix[0][2] - matrix[2][0] ) / scale;
		} 
		// Else the third element of the diagonal is the greatest value
		else 
		{	
			// Find the scale according to the third element, and double that value
			scale  = (float)sqrt( 1.0f + matrix[2][2] - matrix[0][0] - matrix[1][1] ) * 2.0f;

			// Calculate the x, y, x and w of the quaternion through the respective equation
			x = (matrix[0][2] + matrix[2][0] ) / scale;
			y = (matrix[2][1] + matrix[1][2] ) / scale;
			z = 0.25f * scale;
			w = (matrix[1][0] - matrix[0][1] ) / scale;
		}
	}
	QuatNorm q(x, y, z, -w);
	q.normalize();
	return q;
}
*/
	
QuatNorm matToQuat(Mat3d &mati) {
	QuatNorm q;
	Mat3d mat = mati.transpose();

	double d0 = mat[0][0], d1 = mat[1][1], d2 = mat[2][2];
	double xx = 1.0 + d0 - d1 - d2;               // from the diagonal of rotation
	double yy = 1.0 - d0 + d1 - d2;               // matrix, find the terms in
	double zz = 1.0 - d0 - d1 + d2;               // each Quaternion compoment
	double rr = 1.0 + d0 + d1 + d2;

	double max = rr;                              // find the maximum of all
	if (xx > max) max = xx;                               // diagonal terms.
	if (yy > max) max = yy;
	if (zz > max) max = zz;

	if (rr == max) {
		double r4 = sqrt(rr * 4.0);
		q.x = (mat[1][2] - mat[2][1]) / r4;     // find other components from
		q.y = (mat[2][0] - mat[0][2]) / r4;     // off diagonal terms of
		q.z = (mat[0][1] - mat[1][0]) / r4;     // rotation matrix.
		q.w = r4 / 4.0;
	}
	else if (xx == max) {
		double x4 = sqrt(xx * 4.0);
		q.x = x4 / 4.0;
		q.y = (mat[0][1] + mat[1][0]) / x4;
		q.z = (mat[0][2] + mat[2][0]) / x4;
		q.w = (mat[1][2] - mat[2][1]) / x4;
	}
	else if (yy == max) {
		double y4 = sqrt(yy * 4.0);
		q.x = (mat[0][1] + mat[1][0]) / y4;
		q.y =  y4 / 4.0;
		q.z = (mat[1][2] + mat[2][1]) / y4;
		q.w = (mat[2][0] - mat[0][2]) / y4;
	}
	else {
		double z4 = sqrt(zz * 4.0);
		q.x = (mat[0][2] + mat[2][0]) / z4;
		q.y = (mat[1][2] + mat[2][1]) / z4;
		q.z =  z4 / 4.0;
		q.w = (mat[0][1] - mat[1][0]) / z4;
	}
	q.w = -q.w;
	q.normalize();
	return q;
}

Vec3d matToEuler(Mat4d m) {
	// adapted from http://www.darwin3d.com/gamedev/quat2eul.cpp
	Vec3d euler;
	double sy, cy, cx, sx, cz, sz;

	sy = -m[2][0];
	cy = sqrt(1.0 - (sy * sy));
	euler[1] = atan2(sy,cy);

	if (fabs(fabs(sy)-1.0) > 1e-8) {
		cx = m[2][2] / cy;
		sx = m[2][1] / cy;
		euler[0] = atan2(sx,cx);

		cz = m[0][0] / cy;
		sz = m[1][0] / cy;
		euler[2] = atan2(sz,cz);
	}
	else {
		// if cos(y)=0, do things differently
		cx = m[1][1];
		sx = -m[1][2];
		euler[0] = atan2(sx,cx);

		cz = 1.0f;
		sz = 0.0f;
		euler[2] = atan2(sz,cz);
	}
	return euler;
}

Vec3d matToEulerYZX(Mat4d m) {
	Vec3d euler;
	double sy, cy, cx, sx, cz, sz;

	sz = m[1][0];
	cz = sqrt(1.0 - (sz * sz));
	euler[1] = atan2(sz, cz);

	if (fabs(fabs(sz)-1.0) > 1e-8) {
		cx = m[1][1] / cz;
		sx = -m[1][2] / cz;
		euler[2] = atan2(sx,cx);

		cy = m[0][0] / cz;
		sy = -m[2][0] / cz;
		euler[0] = atan2(sy,cy);
	}
	else {
		// if cos(z)=0, do things differently
		euler[2] = asin(sqrt(0.5 * (m[2][2] + 1)));

		cy = 1.0f;
		sy = 0.0f;
		euler[0] = atan2(sy,cy);
	}
	return euler;
}

Vec3d matToEulerXZY(Mat4d m) {
	Vec3d euler;
	double sy, cy, cx, sx, cz, sz;

	sz = -m[0][1];
	cz = sqrt(1.0 - (sz * sz));
	euler[1] = atan2(sz, cz);

	if (fabs(fabs(sz)-1.0) > 1e-8) {
		cx = m[1][1] / cz;
		sx = m[2][1] / cz;
		euler[0] = atan2(sx,cx);

		cy = m[0][0] / cz;
		sy = m[0][2] / cz;
		euler[2] = atan2(sy,cy);
	}
	else {
		// if cos(z)=0, do things differently
		cx = m[1][0];
		sx = -m[1][2];
		euler[0] = atan2(sx,cx);

		cy = 1.0f;
		sy = 0.0f;
		euler[2] = atan2(sy,cy);
	}
	return euler;
}

void glbDirectedCyl(Vec3d dir, double len, double botWidth, double topWidth) {
	static GLUquadricObj *cyl = NULL;

	if (cyl == NULL)
		cyl = gluNewQuadric();

	dir.normalize();
	Vec3d normal = dir ^ Vec3d(0, 0, 1);

	glPushMatrix();
		if (normal.length() < 0.001)
			normal = Vec3d(1, 0, 0);
		glRotated(-acos(dir[2]) * RAD_TO_DEG, normal[0], normal[1], normal[2]);
		gluCylinder(cyl, botWidth, topWidth, len, 16, 2);
	glPopMatrix();
}

void glbSphere(Vec3d pos, double rad, int quality) {
	static GLUquadricObj *sphere = NULL;

	if (sphere == NULL)
		sphere = gluNewQuadric();

	glPushMatrix();
		pos.glTranslate();
		gluSphere(sphere, rad, quality, quality);
	glPopMatrix();
}

void glbBox(Vec3d pos, Vec3d rad) {
	// top and bottom
	glNormal3d(0,0,-1);
	glBegin(GL_POLYGON);
		glVertex3d(pos[0]-rad[0], pos[1]-rad[1], pos[2]-rad[2]);
		glVertex3d(pos[0]-rad[0], pos[1]+rad[1], pos[2]-rad[2]);
		glVertex3d(pos[0]+rad[0], pos[1]+rad[1], pos[2]-rad[2]);
		glVertex3d(pos[0]+rad[0], pos[1]-rad[1], pos[2]-rad[2]);
	glEnd();
	glNormal3d(0,0,1);
	glBegin(GL_POLYGON);
		glVertex3d(pos[0]-rad[0], pos[1]-rad[1], pos[2]+rad[2]);
		glVertex3d(pos[0]+rad[0], pos[1]-rad[1], pos[2]+rad[2]);
		glVertex3d(pos[0]+rad[0], pos[1]+rad[1], pos[2]+rad[2]);
		glVertex3d(pos[0]-rad[0], pos[1]+rad[1], pos[2]+rad[2]);
	glEnd();
	// left and right
	glNormal3d(0,-1,0);
	glBegin(GL_POLYGON);
		glVertex3d(pos[0]-rad[0], pos[1]-rad[1], pos[2]-rad[2]);
		glVertex3d(pos[0]+rad[0], pos[1]-rad[1], pos[2]-rad[2]);
		glVertex3d(pos[0]+rad[0], pos[1]-rad[1], pos[2]+rad[2]);
		glVertex3d(pos[0]-rad[0], pos[1]-rad[1], pos[2]+rad[2]);
	glEnd();
	glNormal3d(0,1,0);
	glBegin(GL_POLYGON);
		glVertex3d(pos[0]-rad[0], pos[1]+rad[1], pos[2]-rad[2]);
		glVertex3d(pos[0]-rad[0], pos[1]+rad[1], pos[2]+rad[2]);
		glVertex3d(pos[0]+rad[0], pos[1]+rad[1], pos[2]+rad[2]);
		glVertex3d(pos[0]+rad[0], pos[1]+rad[1], pos[2]-rad[2]);
	glEnd();
	// front and back
	glNormal3d(-1,0,0);
	glBegin(GL_POLYGON);
		glVertex3d(pos[0]-rad[0], pos[1]-rad[1], pos[2]-rad[2]);
		glVertex3d(pos[0]-rad[0], pos[1]-rad[1], pos[2]+rad[2]);
		glVertex3d(pos[0]-rad[0], pos[1]+rad[1], pos[2]+rad[2]);
		glVertex3d(pos[0]-rad[0], pos[1]+rad[1], pos[2]-rad[2]);
	glEnd();
	glNormal3d(1,0,0);
	glBegin(GL_POLYGON);
		glVertex3d(pos[0]+rad[0], pos[1]-rad[1], pos[2]-rad[2]);
		glVertex3d(pos[0]+rad[0], pos[1]+rad[1], pos[2]-rad[2]);
		glVertex3d(pos[0]+rad[0], pos[1]+rad[1], pos[2]+rad[2]);
		glVertex3d(pos[0]+rad[0], pos[1]-rad[1], pos[2]+rad[2]);
	glEnd();
}

void glbAxes() {
	Vec3d v;

	v = Vec3d(1, 0, 0);
	v.glColor();
	glbDirectedCyl(v, 1, 0.05, 0.05);
	v.glTranslate();
	glbDirectedCyl(v, 0.1, 0.1, 0.0);
	(-v).glTranslate();

	v = Vec3d(0, 1, 0);
	v.glColor();
	glbDirectedCyl(v, 1, 0.05, 0.05);
	v.glTranslate();
	glbDirectedCyl(v, 0.1, 0.1, 0.0);
	(-v).glTranslate();

	v = Vec3d(0, 0, 1);
	v.glColor();
	glbDirectedCyl(v, 1, 0.05, 0.05);
	v.glTranslate();
	glbDirectedCyl(v, 0.1, 0.1, 0.0);
	(-v).glTranslate();
}

long baTimer() {
#ifdef WIN32
	static long curTime = -1;

	struct _timeb timebuffer;
	long newTime, retTime;

	_ftime(&timebuffer);
	newTime = timebuffer.time * 1000 + timebuffer.millitm;
	if (curTime == -1)
		retTime = 0;
	else
		retTime = newTime - curTime;

	curTime = newTime;
	return retTime;
#else
	return 0;
#endif
}

void initRand() {
	srand( (unsigned)time(NULL) );
}

double boundedRand( double lo, double hi ) {
	double	r = (double)rand();

	// get a number between 0.0 and 1.0
	r /= double(RAND_MAX);

	r *= hi - lo;
	r += lo;

	return r;
}

bool openIFStream(ifstream *is, const char *fname, const char *msg) {
	is->clear();
	is->open(fname);

	if (*is)
		return true;

	if (msg)
		cerr << "WARNING: can't open " << msg << " '" << fname << "' for reading" << endl;
	else
		cerr << "WARNING: can't open '" << fname << "' for reading" << endl;

	return false;
}

bool openFile(FILE **f, const char *fname, const char *mode, const char *msg) {
	*f = fopen(fname, mode);

	if (*f)
		return true;

	if (msg)
		cerr << "WARNING: can't open " << msg << " '" << fname << "'" << endl;
	else
		cerr << "WARNING: can't open '" << fname << "'" << endl;

	return false;
}

bool openOFStream(ofstream *os, const char *fname, const char *msg) {
	os->clear();
	os->open(fname);

	if (*os)
		return true;

	if (msg)
		cerr << "WARNING: can't open " << msg << " '" << fname << "' for writing" << endl;
	else
		cerr << "WARNING: can't open '" << fname << "' for writing" << endl;

	return false;
}

bool baAssert(bool test, char *message, bool fatal) {
	if (!test) {
		if (message) {
			if (fatal)
				cerr << "FATAL ASSERTION ERROR: ";
			else
				cerr << "ASSERTION ERROR: ";
			cerr << message << endl;
		}

		if (fatal)
			exit(0);
	}

	return test;
}

Vec3d hotCold(double val, double min, double max) {
	val = (val - ((max + min) / 2.0)) / (0.25 * (max - min));
	if (val <= -2.0)
		return Vec3d(0, 1, 1);
	else if (val < -1.0)
		return Vec3d(0, 1, -(val + 1));
	else if (val < 0)
		return Vec3d((val + 1), 1, (val + 1));
	else if (val < 1)
		return Vec3d(1, (1 - val), (1 - val));
	else if (val < 2)
		return Vec3d(1, (val - 1), 0);
	else
		return Vec3d(1, 1, 0);
}


void rstrip(char *s) {
	int i = (int)strlen(s) - 1;
	while (i >= 0 && s[i] <= ' ') {
		s[i] = 0;
		i--;
	}
}

double lerp(double d1, double d2, double v) {
	return (1.0 - v) * d1 + v * d2;
}

bool sphereIntersection(Vec3d pnt, Vec3d dir, Vec3d center, double radius) {
	double t;
	//Calculate the coefficients for the quadratic that defines where the line intersects the sphere
	double a = dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2];
	double b = 2*pnt[0]*dir[0] - 2*dir[0]*center[0] + 2*pnt[1]*dir[1] - 2*dir[1]*center[1] +
			   2*pnt[2]*dir[2] - 2*dir[2]*center[2];
	double c = pnt[0]*pnt[0] - 2*pnt[0]*center[0] + center[0]*center[0] +
			   pnt[1]*pnt[1] - 2*pnt[1]*center[1] + center[1]*center[1] +
			   pnt[2]*pnt[2] - 2*pnt[2]*center[2] + center[2]*center[2] - 
			   radius*radius;
	//if the determinant of the quadratic is greater than or equal to 0, there is an intersection, else not
	return b*b - 4*a*c >= 0;
}
