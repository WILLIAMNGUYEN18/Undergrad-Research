#ifndef U_SOLVER_H
#define U_SOLVER_H

#include <vector>
using namespace std;
#include "vec.h"
#include "solver.h"
#include "uSkin.h"

class UDataSet;
class TriMesh;
class EdgeList;
class VLVecd;
class VLMatd;
class MarkerSet;


class UExample {
public:
	char fname[80];
	int character;
	char charName[5], poseCh;
	bool minPose;

	int numPts;
	double *ptsConf;
	Vec3d *pts;
	Vec3d *normals;

	TriMesh *mesh;
	char *vertList;
	EdgeList *edgeList;

	int *lookup;
	int luSize;
	Vec3d *luOffsets;

	UExample();
	void init(int iNumPts);
	void init(const char *meshName);
	void initLookup(int *lu, int lus);

	bool save(const char *fname);
	bool load(const char *fname);

	bool getPt(int ind, Vec3d *v, double *conf = NULL);
};

class UDataSet {
public:
	int numExamples;
	int numCharacters;
	UExample *examples;

	int *charIndex;
	VLVecd *charMu;
//	VLMatd *charPhi;
	double phiLogEntropy;

	UDataSet();

	void init(int iNumExamples, int iNumCharacters);
	void initLookup(MarkerSet &mrefs);
};

class NeighborRelation {
public:
	int v0, v1, ind0, ind1;
	double dist;
	void set(double sDist, int sV0, int sInd0, int sV1=-1, int sInd1 = 0) {
		dist = sDist; v0 = sV0; ind0 = sInd0; v1 = sV1; ind1 = sInd1;
	}
};

class USolver : public IDifferentiableFunction {
protected:
	int dressDofPos, weightDofPos;

public:
	UDataSet *dataSet;
	USkin *skin;
	TriMesh *uMesh;
	int numComponents, numPDComponents;
	double lastErr;
	double *markerAssignments;
	int numMarkers;
	double *pddMask;
	Vec3d *dfdv;

	int numOrigPts, numPts;
	int *mirrorMap, *mirrorTrans;
	int maxSolveEx;
	int singleComp, minComp, maxComp, maxPComp;

	vector<int> *weightPriorVerts;
	vector<double> *weightPriorCoeffs;

	// regularization sigmas
	double numWeightNeigh, numPddNeigh, numPdnNeigh, numDressMuNeigh, numDressWNeigh, numVert, numNorm, numNM, numMarkerSpring, numNMNeigh;
	double weightNeighSigma, pddNeighSigma, pdnNeighSigma, dressMuNeighSigma, dressWNeighSigma, vertSigma, normSigma, nmSigma, markerSpringSigma, nmNeighSigma;
	double weightNeighErr, pddNeighErr, pdnNeighErr, dressMuNeighErr, dressWNeighErr, vertErr, normErr, nmErr, markerSpringErr, nmNeighErr;
	// other sigmas
	double minPoseSigma, pddSigma;

	// neighbor table
	vector<int> *neighbors;
	vector<NeighborRelation> neighTable, pddNeighTable, pdnNeighTable;

	// marker springs
	vector<int> mSpringVerts;
	vector<int> mSpringTrans;
	vector<double> mSpringWeights;
	vector<Vec3d> mSpringVerts2;
	int mSpringNumEx;

	// variables
	bool optDress, optInt, optPose, optWeight, optPDD, optX, optNM, optNMP;
	Vecd vDress, vInt, vPose, vWeight, vPDD, vNMP;
	Vec3d *vNM;
	int dpDress, dpInt, dpPose, dpWeight, dpPDD, dpX, dpNM, dpNMP;
	Vecd curVars, grad;
	bool globalPoseOnly, lockPoseZero, ignorePoints;

	USolver();
	void init(UDataSet *iDataSet, USkin *iSkin, TriMesh *iMesh, const char *mirrorFN, int iNumComponents, int iNumPDComponents);
	void dumpDofPos();

	void save(FILE *f);
	void load(FILE *f);
	bool saveDress(char *fname);
	bool loadDress(char *fname, char *fname2 = NULL);
	bool savePoses(char *fname, int mode = 0);
	bool loadPoses(char *fname, bool intOnly = false);

	void weightsToVars();
	void dressToVars(int comp = -1);

	void updateSkel(int ex);
	void updateSkel(double *comps, int numComps, int poseInd = -2);
	void updateWeights();
	void updatePoints(int ch);
	void updatePoints(double *comps, int numComps, bool pddOnly = false);

	void rayPoint(Vec3d pt, Vec3d norm, UExample *ex, double &curDistance, Vec3d &curClosestPt, double &curSurfaceWeight);

	void adjustWeights();

//	double evaluateFunctionExp();

	void buildNeighborTable(TriMesh *tm);
	void buildPDNeighTable();
	void initWeightStencil();
	void initVars(bool verbose = false);
	void toVars();
	void fromVars(Vecd &vars);
//	void calcPtErrAndGrad(int pt, Vec3d &target, int ex, double mult);
	void applyDfdv(int ex);
	void calcNMErr(Vec3d &curVal, Vec3d &targetVal, int ex, int pt, int gradOfs, double yMult);
	virtual double evaluateFunction(Vecd& variables);
	virtual void evaluateGradient(Vecd& variables, Vecd& gradient);
	virtual void solverStep();
	void evaluateWeightPrior();
	void evaluateNMPrior();
	void evaluatePDDPrior();
	void evaluatePDNPrior();

	void renormalizeX();
	void autoWeight();
	void eStep();

	// accessors
	inline Vec3d varDressPt(int pt, int comp) {
		int ofs;
		if (pt >= numOrigPts) {
			ofs = (comp * numOrigPts + mirrorMap[pt]) * 3;
			return Vec3d(vDress[ofs+0], -vDress[ofs+1], vDress[ofs+2]);
		}
		ofs = (comp * numOrigPts + pt) * 3;
		return Vec3d(vDress[ofs+0], vDress[ofs+1], vDress[ofs+2]);
	}

	inline double &varWeight(int pt, int inf) {
		return vWeight[pt*(skin->maxInf-1) + inf];
	}
};


#endif
