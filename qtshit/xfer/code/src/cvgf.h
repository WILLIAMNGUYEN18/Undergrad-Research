#ifndef _CVGF_
#define _CVGF_

#include "charVec.h"
#include "solver.h"

class TriMesh;
class LBFGSSolver;
class Skin;
class Keeper;
class VLMatd;
class VLVecd;


class CVExample {
public:
	Vec3d *points;
	double *conf;
	int numPoints;
	Mat4d *trans, *iTrans;
	int numTPoints;
	double *dofs;
	int numDofs;

	int charID;
	char fname[80];

	void init(int nPoints);
	void buildExample(TriMesh *tm, Skin *skin);
	void save(char *fname);
	void load(char *fname);
};

class CVGoalFunction : public IDifferentiableFunction {
public:
	CharVec *cv;
	CVExample *examples;
	int numChars;

	double distDeviation;
	double matchDeviation;

	VLVecd *curMu;
	VLMatd *curPhi;
	double curDErr, curMErr, curErr;

	bool lockShape;
	bool lockPose;
	bool useCovariance, useRegularization;

	int numExamples;

	LBFGSSolver *solver;
	int solvingPt, solverXYZ;
	int solvingPtIndex;
	Vecd solveGrad;

	CVGoalFunction();

	// set examples before calling:
	void init(int iNumExamples, int iNumChars);

	int numObsPts();

	void runEM(int numIterations);
	void eStep();

	void optimizeWb(bool test = false);
	void regularize();
	void normalizeN();

	virtual double evaluateFunction(Vecd& variables);
	virtual void evaluateGradient(Vecd& variables, Vecd& gradient);

	double calcFreeEnergy(bool verbose);
};

class CVSkinningGF : public IDifferentiableFunction {
public:
	CVGoalFunction *cvGF;
	Skin *skin;
	Skeleton *matchPoses;
	Vecd vars, grad;
	int numVars;
	double lastErr;


	void init(CVGoalFunction *gf, Skin *sk, Skeleton *mp);
	void varsToSkin(Vecd &variables);

	virtual double evaluateFunction(Vecd& variables);
	virtual void evaluateGradient(Vecd& variables, Vecd& gradient);
	virtual void solverStep();
};

#endif
