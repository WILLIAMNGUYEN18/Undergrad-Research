#ifndef SKELETON_H
#define SKELETON_H

#include "saveload.h"
#include "kinematics.h"

class SkelTransform;

const int MAX_NAME_LEN = 40;
const int COPY_POSE = 1;
const int COPY_INT = 2;
const int COPY_ALL = 3;
const int MAX_COMBINED = 3;

class Skeleton : public SLInterface {
public:
	NameTableT<SkelTransform>transforms;
	int numIntrinsicDofs, numPoseDofs, numDofs;
	bool *dofIntrins;

	Skeleton();

	void init();
	static Skeleton* load(char *fname);

	void copyVals(Skeleton *sk, int mode = COPY_ALL);
	void interpVals(Skeleton *sk0, Skeleton *sk1, double interp, bool poseOnly = false);
	void updateCoords();
	void updateDerivs();
	void updateGlobalDerivs();

	void drawGL(SkelTransform *tr = NULL, double alpha = 1.0);
	void drawRIB(ostream &rib);

	void zero();
	void loadPose(istream &in);
	void savePose(ostream &out);

	static void registerClass(SaveLoad *sl);
	static SLInterface *newInstance(int);

	void copyFrom(Skeleton *skel);
	void mirrorLR();

	void allocDerivs(int numDofs);
};

class SkelTransform : public SLInterface {
public:
	int index;
	int dofInd;
	char name[MAX_NAME_LEN];
	Vec3d color;
	bool isIntrinsic;

	Mat4d *globalDerivs;

	KCoord curCoord, globalCoord;

	char parent[MAX_NAME_LEN];

	SkelTransform *parentPtr;
	vector<SkelTransform *> children;

	SkelTransform();
	virtual ~SkelTransform();

	void addChild(SkelTransform *st);

	virtual void initRefs(Skeleton *skel) { }
	virtual int numDofs();
	virtual void loadDofs(double *v, double scale=1.0);
	virtual void unloadDofs(double *v, double scale=1.0);
	virtual double &getDofAddr(int dofI);
	virtual void updateCoord() = 0;
	virtual void updateDerivs() { };
	virtual void updateGlobalDerivs(int maxDof);

	virtual void copyVal(SkelTransform *k) { }
	virtual void interpVal(SkelTransform *k0, SkelTransform *k1, double interp) { }
	virtual void mirrorVal(SkelTransform *k) { }
	virtual void zero() { }
	virtual void loadPose(istream &in);
	virtual void savePose(ostream &out);
	virtual double distance(SkelTransform *k) { return 0; }
	virtual void normalize() { }

	virtual void drawGL(double alpha = 1.0) { }
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);

	static void registerProps(SLClassInfo *ci);

	virtual void copyFrom(SkelTransform *s);
	virtual SkelTransform *clone()=0;
};

class SkelVec3d : public SkelTransform {
public:
	Vec3d curVal;
};

class SkelTranslation : public SkelVec3d {
public:
	static const int CONSTANT_X;
	static const int CONSTANT_Y;
	static const int CONSTANT_Z;

	int xyz;
	Vec3d box;

	SkelTranslation();

	virtual int numDofs();
	virtual void loadDofs(double *v, double scale=1.0);
	virtual void unloadDofs(double *v, double scale=1.0);
	virtual double &getDofAddr(int dofI);
	virtual void updateCoord();
	virtual void updateDerivs();

	virtual void copyVal(SkelTransform *k);
	virtual void interpVal(SkelTransform *k0, SkelTransform *k1, double interp);
	virtual void zero();
	virtual void loadPose(istream &in);
	virtual double distance(SkelTransform *k);

	virtual void drawGL(double alpha = 1.0);
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);

	static void registerClass(SaveLoad *sl);
	static SLInterface *newInstance(int);

	virtual void copyFrom(SkelTransform *s);
	virtual SkelTransform *clone();
};

class SkelQuatRotation : public SkelTransform {
public:
	QuatNorm curQuat;

	SkelQuatRotation();

	virtual int numDofs();
	virtual void loadDofs(double *v, double scale=1.0);
	virtual void unloadDofs(double *v, double scale=1.0);
	virtual double &getDofAddr(int dofI);
	virtual void updateCoord();
	virtual void updateDerivs();

	virtual void copyVal(SkelTransform *k);
	virtual void interpVal(SkelTransform *k0, SkelTransform *k1, double interp);
	virtual void mirrorVal(SkelTransform *k);
	virtual void zero();
	virtual void loadPose(istream &in);
	virtual void normalize();
	virtual double distance(SkelTransform *k);

	virtual void drawGL(double alpha = 1.0);
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);

	static void registerClass(SaveLoad *sl);
	static SLInterface *newInstance(int);

	virtual void copyFrom(SkelTransform *s);
	virtual SkelTransform *clone();
};

class SkelEulerRotation : public SkelTransform {
public:
	int axis;
	double curAngle;
	double min, max;
	QuatNorm curQuat;

	SkelEulerRotation();

	virtual int numDofs();
	virtual void loadDofs(double *v, double scale=1.0);
	virtual void unloadDofs(double *v, double scale=1.0);
	virtual double &getDofAddr(int dofI);
	virtual void updateCoord();
	virtual void updateDerivs();

	virtual void copyVal(SkelTransform *k);
	virtual void interpVal(SkelTransform *k0, SkelTransform *k1, double interp);
	virtual void mirrorVal(SkelTransform *k);
	virtual void zero();
	virtual void loadPose(istream &in);
	virtual void normalize();
	virtual double distance(SkelTransform *k);

	virtual void drawGL(double alpha = 1.0);
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);

	static void registerClass(SaveLoad *sl);
	static SLInterface *newInstance(int);

	virtual void copyFrom(SkelTransform *s);
	virtual SkelTransform *clone();
};

class SkelPolarAxisRotation : public SkelTransform {
public:
	double curPhi, curPsi, curTheta;
	Vec3d curAxis;
	QuatNorm curQuat;
	double min, max;

	SkelPolarAxisRotation();

	virtual int numDofs();
	virtual void loadDofs(double *v, double scale=1.0);
	virtual void unloadDofs(double *v, double scale=1.0);
	virtual double &getDofAddr(int dofI);
	virtual void updateCoord();
	virtual void updateDerivs();

	virtual void copyVal(SkelTransform *k);
	virtual void interpVal(SkelTransform *k0, SkelTransform *k1, double interp);
	virtual void mirrorVal(SkelTransform *k);
	virtual void zero();
	virtual void loadPose(istream &in);
	virtual void normalize();
	virtual double distance(SkelTransform *k);

	virtual void drawGL(double alpha = 1.0);
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);

	static void registerClass(SaveLoad *sl);
	static SLInterface *newInstance(int);

	virtual void copyFrom(SkelTransform *s);
	virtual SkelTransform *clone();
};

class SkelSymmetricTranslation : public SkelVec3d {
public:
	SkelTranslation *orig;
	char origTranslation[MAX_NAME_LEN];
	Vec3d axisMult;
	Vec3d box;

	SkelSymmetricTranslation();

	virtual void initRefs(Skeleton *skel);
	virtual void updateCoord();
	virtual void updateGlobalDerivs(int maxDof);

	void loadPose(istream &in);
	virtual double distance(SkelTransform *k);

	virtual void drawGL(double alpha = 1.0);
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);

	static void registerClass(SaveLoad *sl);
	static SLInterface *newInstance(int);

	virtual void copyFrom(SkelTransform *s);
	virtual SkelTransform *clone();
};

class SkelPartialTransform : public SkelTransform {
public:
	SkelTransform *orig;
	char origTransform[MAX_NAME_LEN];
	double factor;

	SkelPartialTransform();

	virtual void initRefs(Skeleton *skel);
	virtual void updateCoord();
	virtual void updateGlobalDerivs(int maxDof);

	void loadPose(istream &in);
	virtual double distance(SkelTransform *k);

	virtual void drawGL(double alpha = 1.0);
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);

	static void registerClass(SaveLoad *sl);
	static SLInterface *newInstance(int);

	virtual void copyFrom(SkelTransform *s);
	virtual SkelTransform *clone();
};

class SkelCombinedTransform : public SkelTransform {
public:
	SkelTransform *orig[MAX_COMBINED];
	char origTransform[MAX_COMBINED][256];

	SkelCombinedTransform();

	virtual void initRefs(Skeleton *skel);
	virtual void updateCoord();
	virtual void updateGlobalDerivs(int maxDof);

	static void registerClass(SaveLoad *sl);
	static SLInterface *newInstance(int);

	virtual void copyFrom(SkelTransform *s);
	virtual SkelTransform *clone();
};

#endif
