#ifndef SURFACE_DEF_H
#define SURFACE_DEF_H

#include "vec.h"
#include "solver.h"
#include "edgelist.h"
#include "trimesh_util.h"

class TriMesh;
class MarkerSet;
class MeshEdge;

class LADeformationGoalFunction : public IDifferentiableFunction {
public:
	TriMesh *targetTM, *curMesh;
	vector <int> *curMeshNeigh;
	Vec3d *origVerts;
	double lastErr;
	double *restAreas;

	vector<int> markerRefs;
	MarkerSet *markers;

	int capRate, stepCount;

	bool enableSurfaceMatch, enableSmoothness, enableMarkerMatch;
	bool smatchShowError;
	double surfaceMatchWeight, smoothnessWeight, markerMatchWeight;
	double *neighWeights, *surfWeights;

	double curDistance;
	Vec3d curGradient;
	Vec3d curClosestPt;
	double curSurfaceWeight;

	char *vertList;
	EdgeList *edgeList;

	Vecd vars, grad;
	int varsPerVert;

	bool lockShape;

	double bendWeight;
	
	LADeformationGoalFunction(TriMesh *mesh, vector <int>* neigh = NULL);
	~LADeformationGoalFunction();

	void prepareTriMesh(TriMesh *dtm);
	virtual void applyDef(Vecd &variables);
	virtual void zeroDeformation();

	void updateCurrent(int pt);
	virtual void addGradientVec(int index, Vec3d v);

	virtual double evaluateFunction(Vecd& variables);
	virtual void evaluateGradient(Vecd& variables, Vecd& gradient);
	virtual void solverStep();
};

class EdgeMatchGF : public LADeformationGoalFunction {
public:
	TriMesh *edgeMatchTM;
	int numEdges;
	MeshEdge *templateEdges;
	vector<TMNeigh> *geodesics;
	Vec3d *undefVerts;
	double *conf;
	bool *restrict;
	double *baryRecon;
	int *baryPts;

	EdgeMatchGF(TriMesh *mesh, vector <int>* neigh, TriMesh *matchTM);

	void unfold(int reps);

	virtual void applyDef(Vecd &variables);
	virtual void zeroDeformation();

	virtual void addGradientVec(int index, Vec3d v);

	virtual double evaluateFunction(Vecd& variables);
};

class Skin;
class Skeleton;

/*
class SkinMatchGF : public LADeformationGoalFunction {
public:
	int numExamples;
	int *varIndex;
	Skin *skin;
	Skeleton *poses;
	TriMesh *meshArray;

	char **vertListArray;
	EdgeList **edgeListArray;
	EdgeList templateEdgeList;
	TriMesh *edgeMatchTM;

	int numVars;

	SkinMatchGF(int numEx, TriMesh *mesh, Skin *sk, Skeleton *mp);

	void prepareTriMesh(TriMesh *dtm);

	virtual void applyDef(Vecd &variables);
	virtual void zeroDeformation();

	virtual void addGradientVec(int index, Vec3d v);

	virtual double evaluateFunction(Vecd& variables);
};
*/
class MeshEdge {
public:
	int v0, v1;
	int f0, f1;
	int opp0, opp1;
	Vec3d n0, n1, e;
	double n0Len, n1Len;
	double eCurLen, eCurAngle;
	double eRestLen, eRestAngle;
	double he0, he1;
	double restBendFactor;

	MeshEdge();
	void update(TriMesh *tm);
	void update(Vec3d &x0, Vec3d &x1, Vec3d &x2, Vec3d &x3);
	void init(TriMesh *tm);
};

MeshEdge *buildEdges(TriMesh *tm, int &numEdges);

#endif