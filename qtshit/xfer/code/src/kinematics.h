#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "mat.h"
#include "vec.h"
#include "quatnorm.h"

// ------------- Coord -----------------

class KCoord {
public:
	Mat4d mat;
	QuatNorm q;
	Vec3d v;

	typedef Mat4d *Mat4dPtr;
	int numDofs;
	Mat4dPtr *deriv;

	KCoord();
	~KCoord();
	void initDeriv(int iNumDofs);
	void setDeriv(int dof, Mat4d &val);

	//friend KCoord operator*(const KCoord& a, const KCoord& b);
	KCoord operator*(const KCoord& a);
	KCoord& operator=(const KCoord& a);
};


/*
#include <iostream>
#include <vector>
using namespace std;

#include <FL/Fl_Output.h>
#include <FL/FL_Value_Input.h>
#include "vec.h"
#include "mat.h"
#include "quatnorm.h"

class KinematicModel;

// ------------- Data ------------------

class KDofSet {
public:
	Vecd data;

	int globalSize;
	int perFrameSize;
	int numFrames;
	int numDofs;
	int nonMarkerSize;

	KDofSet();
	~KDofSet();

	void flush();
	KDofSet *clone();
	void copyVals(KDofSet *d);

	int calcIndex(int frame, int ofs);
	double &v(int ofs);
	double &v(int frame, int ofs);

	friend ostream& operator <<(ostream& os, KDofSet &k);
	friend istream& operator >>(istream& is, KDofSet &k);
};

class KMarkerData {
public:
	Vec3d *data;
	Vec3d *colors;
	int numMarkers;
	int numFrames;
	bool variableSize;
	bool colorData;

	KMarkerData();
	~KMarkerData();

	void flush();
	void init(int nMarkers, int nFrames);
	Vec3d &v(int frame, int marker);
	Vec3d &c(int frame, int marker);

	friend ostream& operator <<(ostream& os, KMarkerData &k);
	friend istream& operator >>(istream& is, KMarkerData &k);
};

// ------------- Marker ----------------

class KMarker {
public:
	int dofIndex;
	int markerIndex;

	Vec3d curVal, curPos;
	int transform;
	KinematicModel *model;

	void loadDofs(KDofSet *dofs);
	void unloadDofs(KDofSet *dofs);
	void updatePos();

	friend ostream& operator <<(ostream& os, KMarker &k);
	friend istream& operator >>(istream& is, KMarker &k);
};


// ------------- Transform -------------

class KTransform {
public:
	bool isCopy;
	bool perFrame;

	int index, dofGIndex, dofFIndex;
	char name[40];
	Vec3d color;

	KCoord curCoord;

	int parent;
	vector<int> children;
	KinematicModel *model;

	// scaling
	double lo, hi, importance;

	bool isConst;

	KTransform();
	virtual ~KTransform();

	void addChild(int index);

//	Mat4d &getTransform(int frame);
//	Mat4d &getDeriv(int frame, int dof);

	virtual int type() = 0;
	virtual int numTotalDofs();
	virtual int numGDofs() = 0;
	virtual int numFDofs() = 0;
	virtual void loadRaw(double *raw) = 0;
	virtual void unloadRaw(double *raw) = 0;
	virtual void loadDofs(KDofSet *dofs, int frame);
	virtual void unloadDofs(KDofSet *dofs, int frame) ;
	virtual void updateTransform() = 0;
	virtual void updateDerivs() = 0;

	void loadAllDofs(KDofSet *dofs, int frame);
	void unloadAllDofs(KDofSet *dofs, int frame);
	void calcRecursiveTransforms(KCoord *coords);
	void calcRecursiveDerivs(KCoord *coords, int frame);

	virtual void copyVal(KTransform *k) { }
	virtual void normalize() { }

	virtual void put(ostream &os);
	virtual void get(istream &is);
	virtual void putVal(ostream &os);
	virtual void getVal(istream &is);
	friend ostream& operator <<(ostream& os, KTransform &k);
	friend istream& operator >>(istream& is, KTransform &k);

	virtual void renderSkel(bool showMarkers, bool fancy = false);
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);
	void renderAllSkel(bool showMarkers, bool fancy = false);
	void renderAllSkelRIB(bool showMarkers, ostream &rib);

	virtual void updateUI(Fl_Value_Input *val0, Fl_Value_Input *val1, 
		Fl_Value_Input *val2, Fl_Value_Input *val3, Fl_Output *out, int frame);

	double scaleToDof(double x);
	double scaleFromDof(double x);
	double scaleDeriv();
};

typedef KTransform *KTransformPtr;

const int KT_TRANSLATION = 0;
const int KT_QUATROTATION = 1;
const int KT_EULERROTATION = 2;
const int KT_POLARAXISROTATION = 3;
const int KT_SYMMETRICTRANSLATION = 4;

class KTranslation : public KTransform {
public:
	static const int CONSTANT_X;
	static const int CONSTANT_Y;
	static const int CONSTANT_Z;
	int constAxes;
	Vec3d constVal;
	Vec3d curVal;

	KTranslation();

	virtual int type();
	virtual int numGDofs();
	virtual int numFDofs();
	virtual void loadRaw(double *raw);
	virtual void unloadRaw(double *raw);
	virtual void updateTransform();
	virtual void updateDerivs();

	virtual void copyVal(KTransform *k);

	virtual void putVal(ostream &os);
	virtual void getVal(istream &is);

	virtual void renderSkel(bool showMarkers, bool fancy = false);
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);

	virtual void updateUI(Fl_Value_Input *val0, Fl_Value_Input *val1, 
		Fl_Value_Input *val2, Fl_Value_Input *val3, Fl_Output *out, int frame);
};

class KQuatRotation : public KTransform {
public:
	QuatNorm constQuat;
	QuatNorm curQuat;

	virtual int type();
	virtual int numGDofs();
	virtual int numFDofs();
	virtual void loadRaw(double *raw);
	virtual void unloadRaw(double *raw);
	virtual void updateTransform();
	virtual void updateDerivs();

	virtual void copyVal(KTransform *k);
	virtual void normalize();

	virtual void putVal(ostream &os);
	virtual void getVal(istream &is);

	virtual void renderSkel(bool showMarkers, bool fancy = false);
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);

	virtual void updateUI(Fl_Value_Input *val0, Fl_Value_Input *val1, 
		Fl_Value_Input *val2, Fl_Value_Input *val3, Fl_Output *out, int frame);
};

class KEulerRotation : public KTransform {
public:
	int axis;
	double curAngle;
	double constAngle;
	double min, max;
	QuatNorm curQuat;

	KEulerRotation();

	virtual int type();
	virtual int numGDofs();
	virtual int numFDofs();
	virtual void loadRaw(double *raw);
	virtual void unloadRaw(double *raw);
	virtual void updateTransform();
	virtual void updateDerivs();

	virtual void copyVal(KTransform *k);
	virtual void normalize();

	virtual void putVal(ostream &os);
	virtual void getVal(istream &is);

	virtual void renderSkel(bool showMarkers, bool fancy = false);
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);

	virtual void updateUI(Fl_Value_Input *val0, Fl_Value_Input *val1, 
		Fl_Value_Input *val2, Fl_Value_Input *val3, Fl_Output *out, int frame);
};

class KPolarAxisRotation : public KTransform {
public:
	double curPhi, curPsi, curTheta;
	double constPhi, constPsi, constTheta;
	Vec3d curAxis;
	QuatNorm curQuat;
	double min, max;
	bool constAxis;

	KPolarAxisRotation();

	virtual int type();
	virtual int numGDofs();
	virtual int numFDofs();
	virtual void loadDofs(KDofSet *dofs, int frame);
	virtual void unloadDofs(KDofSet *dofs, int frame);
	virtual void loadRaw(double *raw);
	virtual void unloadRaw(double *raw);
	virtual void updateTransform();
	virtual void updateDerivs();

	virtual void copyVal(KTransform *k);
	virtual void normalize();

	virtual void putVal(ostream &os);
	virtual void getVal(istream &is);

	virtual void renderSkel(bool showMarkers, bool fancy = false);
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);

	virtual void updateUI(Fl_Value_Input *val0, Fl_Value_Input *val1, 
		Fl_Value_Input *val2, Fl_Value_Input *val3, Fl_Output *out, int frame);
};

class KSymmetricTranslation : public KTransform {
protected:
	KTranslation *orig();

public:
	int origTranslation;
	Vec3d axisMult;

	KSymmetricTranslation();

	void updateFromOrig();

	virtual int type();
	virtual int numGDofs();
	virtual int numFDofs();
	virtual void loadRaw(double *raw);
	virtual void unloadRaw(double *raw);
	virtual void updateTransform();
	virtual void updateDerivs();

	virtual void copyVal(KTransform *k);

	virtual void putVal(ostream &os);
	virtual void getVal(istream &is);

	virtual void renderSkel(bool showMarkers, bool fancy = false);
	virtual void renderSkelRIB(bool showMarkers, ostream &rib);

	virtual void updateUI(Fl_Value_Input *val0, Fl_Value_Input *val1, 
		Fl_Value_Input *val2, Fl_Value_Input *val3, Fl_Output *out, int frame);
};

// ------------- Model -----------------

class KinematicModel {
public:
	int numMeshes, numMarkers, numTransforms;

	KMarker    *markers;
	KTransformPtr *transforms;

	KDofSet     *dofs;
	KMarkerData *markerData;
	KMarkerData *rawMarkerData;

	vector<KDofSet*> capDofs;

	KMarkerData curMarkerPos;

	typedef char Str80[80];
	Str80 *meshNames;
	Str80 prefix;

	KinematicModel();

	void flush();

	void initMT(int nMeshes, int nMarkers, int nTransforms);

	void putKin(ostream &os);
	void getKin(istream &is);

	void assignDofIndices();

	void loadDofs(int frame);
	void unloadDofs(int frame);

	void calcTransforms();
	void calcDerivs();

	// coord stuff:
	KCoord *coordData;
	void initCoords();
	void flushCoords();
	Mat4d &trans(int transform, int frame);
	Mat4d* &deriv(int transform, int frame, int dof);
	void setDeriv(int transform, int frame, int dof, Mat4d &m);

	void renderMarkers(bool showDeltas, int frame);

	void loadAll(char *lPrefix);
	bool loadKin(char *fname);
	bool saveKin(char *fname);
	bool loadDofs(char *fname);
	bool saveDofs(char *fname);
	bool loadDofs(istream &in);
	bool saveDofs(ostream &out);
	bool loadPose(char *fname);
	bool savePose(char *fname);

	void fitMarkers(int frame);

	void copyVals(KinematicModel *k);
	void normalizeTransforms();

//	KinematicModel *clone();
};
*/
#endif
