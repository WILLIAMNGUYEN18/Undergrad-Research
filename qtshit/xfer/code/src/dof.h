#ifndef DOF_H
#define DOF_H

#include "skeleton.h"
#include "markers.h"
#include "solver.h"
#include "saveload.h"

class SkelTransformDof : public SLInterface {
protected:
	SkelTransform *trans;
	int index;
public:
	SkelTransformDof();

	virtual Mat4d deriv(SkelTransform *t);
	virtual Vec4d transDeriv(SkelTransform *t);

	void setTransform(SkelTransform *t, int ind);
	virtual void initSkel();
	bool hasDependency(SkelTransform *t);

	virtual void skelToDof(double &value);
	virtual void dofToSkel(double value);
};

class SkelTelescopingDof : public SkelTransformDof {
protected:
	SkelTranslation *tr;
	Vec3d direction;

public:
	virtual Mat4d deriv(SkelTransform *t);
	virtual Vec4d transDeriv(SkelTransform *t);

	virtual void initSkel();

	virtual void skelToDof(double &value);
	virtual void dofToSkel(double value);
};

class DofSet : public IDifferentiableFunction, SLInterface {
protected:
public:
	Skeleton *skel, *origSkel;
	static const int LOCAL_TRANS;
	static const int GLOBAL_TRANS;
	static const int MIXED_TRANS;

	vector<SkelTransformDof*> globalTransDofs;
	vector<SkelTransformDof*> localTransDofs;
	
	MarkerSet *calculatedMarkers;
	vector<int> calcToTrueMap;
	vector<MarkerSet*> trueMarkers;

	int capRate, stepCount;
	double lastErr;

	Vecd variables;

	int numFrames, numDofs, numGlobalDofs, numLocalDofs, numTransDofs;

	DofSet();

	void addTrans(char *name, int kind);
	void setSkel(Skeleton *sk);
	void init();
	void initSkel();
	void loadLandmarks(char *fname);
	void skelValsFromOrig(int frame = -1);
	void skelValsToOrig(int frame);

	void interp(NameTableX<double> *weights);

	void skelToDofs(Vecd &values, int frame);
	void dofsToSkel(Vecd &values, int frame);
	void updateTransDerivs();

	void evaluateGradient(Vecd &values, Vecd &gradient);
	double evaluateFunction(Vecd &values);
	void solverStep();

	void drawGL(int frame = -1);

	static void registerClass(SaveLoad *sl);
	static SLInterface *newInstance(int);
};

#endif
