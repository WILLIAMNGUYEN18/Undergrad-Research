/*
	ba.h

	Various convenient constants and functions.

	Brett Allen
	2002
*/

#ifndef BA_H
#define BA_H

#include "mat.h"
#include "vec.h"
#include "quatnorm.h"
#include <fstream>
using namespace std;
#include <stdio.h>

const double PI = 3.14159625358979;
const double RAD_TO_DEG = 180.0 / PI;
const double DEG_TO_RAD = PI / 180.0;
extern double NaN;

typedef char char80[80];

inline double sqr(double x) {
	return x*x;
}

inline int sqr(int x) {
	return x*x;
}

Mat4d invertTransform(Mat4d &m);
QuatNorm matToQuat(Mat3d &mati);
Vec3d matToEuler(Mat4d m);
Vec3d matToEulerXZY(Mat4d m);

void glbDirectedCyl(Vec3d dir, double len, double botWidth, double topWidth);
void glbSphere(Vec3d pos, double rad, int quality = 10);
void glbBox(Vec3d pos, Vec3d rad);
void glbAxes();

long baTimer();

void initRand();
double boundedRand(double lo, double hi);

bool openIFStream(ifstream *is, const char *fname, const char *msg = NULL);
bool openOFStream(ofstream *os, const char *fname, const char *msg = NULL);
bool openFile(FILE **f, const char *fname, const char *mode = "rb", const char *msg = NULL);

bool baAssert(bool test, char *message = NULL, bool fatal = true);

Vec3d hotCold(double val, double min = -1, double max = 1);

double lerp(double d1, double d2, double v);

bool sphereIntersection(Vec3d pnt, Vec3d dir, Vec3d center, double radius);

template <typename T>
inline T baMin(const T a, const T b) {
	if (a < b)
		return a;
	return b;
}

template <typename T>
inline T baMax(const T a, const T b) {
	if (a > b)
		return a;
	return b;
}

void rstrip(char *s);

#ifndef WIN32
#define _finite finite
#endif

#endif
