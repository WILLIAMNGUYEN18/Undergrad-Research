#ifndef MARKERS_H
#define MARKERS_H

#include <iostream>
#include <vector>
using namespace std;

#include "vec.h"
#include "mat.h"
#include "quatnorm.h"
#include "saveload.h"

class KCoord;
class SkelTransform;
class TriMesh;

const int mkrPOS = 0;
const int mkrPT = 1;
const int mkrBARY = 2;

// ------------- MarkerData ------------

class Marker : public SLInterface {
public:
	char name[40];
	Vec3d pos;
	Vec3d color;
	int kind;
	
	int skelFrame;
	double boneOfs;
//	KCoord *coord;
	SkelTransform *trans;

	int baryVerts[3];
	Vec3d baryPos;
	TriMesh *mesh;

	Marker();
	Vec3d curPos();
	Vec3d offset();

	void drawGL();

	static void registerClass(SaveLoad *sl);
	static SLInterface *newInstance(int);
};

class MarkerSet : public SLInterface {
public:
	Marker *markers;
	int numMarkers;

	MarkerSet();

	void init(int numMarkers);
	Vec3d &v(int marker);
	Vec3d &c(int marker);
	Vec3d offset(int marker);
	Vec3d curPos(int marker);

	void drawGL();

	void load(const char *fname);
	void save(const char *fname);

	void loadFromLandmarks(const char *fname);
	void loadText(const char *fname);
	void loadFromMkr(const char *fname);
	void saveToMkr(const char *fname);

	static void registerClass(SaveLoad *sl);
	static SLInterface *newInstance(int);
};

#endif