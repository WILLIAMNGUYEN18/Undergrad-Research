#ifndef SKIN_CALC_H
#define SKIN_CALC_H

#include "mat.h"
#include "skeleton.h"

class TriMesh;
extern TriMesh *scMesh;
extern Skeleton *scSkel;

//#define POINT_TYPE MWESkinCalcPt
//#define POINT_TYPE LBSSkinCalcPt
#define POINT_TYPE RTSkinCalcPt
#define USE_POSITIONS

class SkinCalcPt {
public:
	int numVars;

	int numTInf;
	int *tTransforms;
	double *tWeights;

	SkinCalcPt();
	virtual ~SkinCalcPt();

	virtual void free(bool delFirst = true) { }
	virtual void initTrans(int iNumTInf) = 0;
	virtual void initLocal(Skeleton *skel, Vec3d &v) = 0;
	virtual void copyFromVars(double *vars) = 0;
	virtual void copyToVars(double *vars) = 0;
	virtual double calcGrad(double *grad, Vec3d &target, Skeleton *skel, double weight=1) = 0;

	virtual void updateGlobal(Skeleton *skel) = 0;
	virtual void getPos(Vec3d &v) = 0;

	double getTransWeight(int t);

	virtual void write(FILE *f) = 0;
	virtual void read(FILE *f) = 0;
};

class LBSSkinCalcPt : public SkinCalcPt {
public:
	Vec3d globalPos;
	Vec3d *tLocalPos;

	virtual void free(bool delFirst = true);
	virtual void initTrans(int iNumTInf);
	virtual void initLocal(Skeleton *skel, Vec3d &v);
	virtual void copyFromVars(double *vars);
	virtual void copyToVars(double *vars);
	virtual double calcGrad(double *grad, Vec3d &target, Skeleton *skel, double weight=1);

	virtual void updateGlobal(Skeleton *skel);
	virtual void getPos(Vec3d &v);

	virtual void write(FILE *f);
	virtual void read(FILE *f);
};

class MWESkinCalcPt : public SkinCalcPt {
public:
	Vec3d localPos;
	Mat4d localMat;

	typedef double double12[12];
	double12 *tMat;

	virtual void free(bool delFirst = true);
	virtual void initTrans(int iNumTInf);
	virtual void initLocal(Skeleton *skel, Vec3d &v);
	virtual void copyFromVars(double *vars);
	virtual void copyToVars(double *vars);
	virtual double calcGrad(double *grad, Vec3d &target, Skeleton *skel, double weight=1);

	virtual void updateGlobal(Skeleton *skel);
	virtual void getPos(Vec3d &v);

	virtual void write(FILE *f);
	virtual void read(FILE *f);
};

class RTSkinCalcPt : public SkinCalcPt {
public:
	Mat4d localMat;
	Vec3d localPos;
	double *tPositions;

	virtual void free(bool delFirst = true);
	virtual void initTrans(int iNumTInf);
	virtual void initLocal(Skeleton *skel, Vec3d &v);
	virtual void copyFromVars(double *vars);
	virtual void copyToVars(double *vars);
	virtual double calcGrad(double *grad, Vec3d &target, Skeleton *skel, double weight=1);

	virtual void updateGlobal(Skeleton *skel);
	virtual void getPos(Vec3d &v);

	virtual void write(FILE *f);
	virtual void read(FILE *f);
};

class Skin {
public:
	int numPts;
	POINT_TYPE *points;
	Skeleton *skel;
	Vec3d *baseVerts;

	int numPreMirror;
	int *mirrorMap;

	void init(int iNumPts);
	void updatePoints();
	void updateMesh(TriMesh *tm);

	void renderPoints();

	void save(char *fname);
	void load(char *fname);
};

void loadMarkerAnalysis(char *path, char *maFName);
void setMatchSkelFrame(int f);
void renderMarkers(int frame);
void saveSkelMatrices();

void initSCMesh(char *fname, char *initPose, char *initMesh);
void calcTransInfluences();

#endif