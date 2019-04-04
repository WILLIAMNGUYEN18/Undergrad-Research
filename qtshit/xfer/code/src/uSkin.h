#ifndef U_SKIN_H
#define U_SKIN_H

#include "ba.h"
#include "rbf.h"
#include <vector>

class Skeleton;
class SkelTransform;
class UExample;
class USolver;
class VLMatd;
class USkin;
class TriMesh;

class UPoseDepDef {
public:
	int transInd;
	USkin *skin;
	SkelTransform *transform;
	bool isMirror;

	int numSamples;
	RBF *rbf;

	UPoseDepDef() { isMirror = false; }

	Vec3d *load(USkin *iSkin, TriMesh *mesh, ifstream &f);
	void mirrorFrom(UPoseDepDef *origPDD, int iTransInd);
};


class USkin {
public:
	Skeleton *skel;
	int *dofIndex;

	int numPts;
	Vec3d *dressPts, *pddDressPts;
	Vec3d *dressJoints;
	char80 *ptNames;

	int numTransInit;
	char80 *tiFrames;
	int *tiMarkers;

	Mat4d *curMats;
	QuatNorm *curJointQuat;
	Vec3d *curJointPos;
	Vec3d *curPts;
	Vec4d *curJointDerivs;

	int maxInf;
	int *infJoints;
	double *infWeights;

	// number of transforms:
	int numTrans;
	// for each point, we'll record the offset into the pddKeys array:
	vector<int> *pddNMIndex, *pddPtIndex;
	// offset values for each sample point:
	int numPddNMKeys, numPddPtKeys;
	Vec3d *pddNMKeys, *pddPtKeys;

	vector<UPoseDepDef*> pdds;
	Vec3d *baseNormalMap;

	USkin();
	void init(Skeleton *iSkel, int iNumPts, int iMaxInf);

	void transInit(UExample *ex);

	int getPddJoint(int pt, int inf);
	void getSkinMat(int pt, Mat3d &m, Vec3d &base);

	void updateMats();
	void updateJoints();
	void updateRBFs();
	void updatePts(int minPt = 0, int maxPt = -1, bool updateNM = true);
//	void calcGrad(int pt, Vec3d v, 
//		double *dressGrad, double *intGrad, 
//		double *poseGrad,
//		double *weightGrad, double *pddGrad);
	void calcVecGrad(int pt, 
		Vec3d *dressGrad, Vec3d *intGrad, 
		Vec3d *poseGrad,
		Vec3d *weightGrad, Vec3d *pddGrad);
//	double calcGradExp(int pt, Vec3d v, 
//					 double *dressGrad, double *intGrad,
//					 double *poseGrad, double *weightGrad, double *pddGrad,
//					 VLMatd &phi, USolver *uSolver);
};

#endif