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
		// if cos(y)=0, we're in trouble
		cx = m[1][1];
		sx = -m[1][2];
		euler[0] = atan2(sx,cx);

		cz = 1.0f;
		sz = 0.0f;
		euler[2] = atan2(sz,cz);
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
	int i = strlen(s) - 1;
	while (i >= 0 && s[i] <= ' ') {
		s[i] = 0;
		i--;
	}
}

double lerp(double d1, double d2, double v) {
	return (1.0 - v) * d1 + v * d2;
}
