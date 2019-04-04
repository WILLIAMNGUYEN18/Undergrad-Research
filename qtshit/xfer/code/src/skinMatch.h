#ifndef SKIN_MATCH_H
#define SKIN_MATCH_H

#include "surfaceDef.h"

class USkin;

class SkinMatchGF : public LADeformationGoalFunction {
public:
	TriMesh *edgeMatchTM;
	int numEdges;
	MeshEdge *templateEdges;
	Vec3d *undefVerts, *origDressPts, *targets;
	double *conf;
	bool *restrict;
	double *baryRecon;
	int *baryPts;
	bool findTargets;

	USkin *skin;

	MarkerSet *srcMarkers;

	QuatNorm *ptRots;
	QuatNorm **ptRotDerivs;
	vector<int> *tmNeigh;

	SkinMatchGF(TriMesh *mesh, vector <int>* neigh, USkin *iSkin);
	void newMatch();

	void unfold(int reps);

	virtual void applyDef(Vecd &variables);
	virtual void zeroDeformation();

	virtual void addGradientVec(int index, Vec3d v);

	virtual double evaluateFunction(Vecd& variables);
};

#endif