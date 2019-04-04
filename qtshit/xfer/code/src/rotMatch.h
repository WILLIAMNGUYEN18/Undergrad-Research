#ifndef ROT_MATCH_H
#define ROT_MATCH_H

#include "surfaceDef.h"

Mat3d ptsRot(Vec3d *v0, Vec3d *v1, int numPts);

class RotMatchGF : public LADeformationGoalFunction {
public:
	TriMesh *edgeMatchTM;
	int numEdges;
	MeshEdge *templateEdges;
	vector<TMNeigh> *geodesics;
	Vec3d *undefVerts, *matchVerts, *targets;
	double *conf;
	bool *restrict;
	double *baryRecon;
	int *baryPts;
	bool findTargets;

	MarkerSet *srcMarkers;

	QuatNorm *ptRots;
	QuatNorm **ptRotDerivs;
	vector<int> *tmNeigh;

	RotMatchGF(TriMesh *mesh, vector <int>* neigh, TriMesh *matchTM);
	void newMatch(TriMesh *mesh, TriMesh *matchTM);

	void unfold(int reps);

	virtual void applyDef(Vecd &variables);
	virtual void zeroDeformation();

	virtual void addGradientVec(int index, Vec3d v);

	virtual double evaluateFunction(Vecd& variables);
};

#endif