#ifndef CHAR_VEC
#define CHAR_VEC

#include <vector>
#include "mat.h"
#include "vec.h"
using namespace std;

class Skeleton;
class TriMesh;
class Skin;

class CharVec;

class CharVecPt {
public:
	CharVec *cv;

	// current local position (and template's local position)
	Vec3d localPos, templateLocalPos;
	// current vector given joint DOFs
	Vec3d *curComponents;

	Vecd data, backup;
	int dataParts;
	int numInfluences;
	int *influences;
	double *infWeights;

	CharVecPt();
	~CharVecPt();
	void init(CharVec *iCV, int nInf);
	void initData();

	double *calcMultipliers(double *jointDofs);
	void updateCurComponents(double *jointDofs = NULL);
	void updateLocalPos(double *n, int nSize = -1);

	int getInfIndex(int inf);
};

class CharVec {
public:
	bool vis, update;
	int numPts;
	int numComponents;
	CharVecPt *cvPts;
	Skeleton *skel;

	int numPreMirror;
	int *mirrorMap;

	TriMesh *tm;
	vector<int> *neighbors;

	double *dofMin, *dofMax;

	CharVec() {
		vis = false;
		update = false;
		dofMin = NULL;
		dofMax = NULL;
	}

	void createCylArm(int lSegs, int rSegs, double rad);
	void createFromSkin(Skin &skin, int iNumComponents);
	void createFromMesh(TriMesh *mesh);
	void setFromLocalTM(Mat4d *localIFrames);

	void backupData();

	void updateCurComponents(double *jointDofs);
	void updateLocalPos(double *n, int nSize = -1);
	void updateTM(Mat4d *localFrames = NULL);
	void render(int viewMode, Vec3d bkg);

	static void dofsFromSkel(Skeleton *skel, double *dofs);
	static void CharVec::mirrorDofs(double *dofs);
};

#endif