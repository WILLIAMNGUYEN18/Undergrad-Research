#ifndef PCA_FIT_H
#define PCA_FIT_H

#include "vec.h"
#include "solver.h"
#include "fast_trimesh.h"
#include <vector>
using namespace std;

// PCAData: a class to store the raw PCA components
class PCAData {
public:
	int numComponents, activeComponents;
	float *average;
	float *components;
	float *sigma2;

	void reconstructMesh(FastTriMesh *mesh, Vecd &weights, vector<int> *markerIndices = NULL);
};


// PCAFitGF: a "goal function" class to find the most 
// likely PCA weights to meet the given criteria
class PCAFitGF : public IDifferentiableFunction {
public:
	double lastErr;

	double markerVariance;

	Vecd vars, grad;

	FastTriMesh *curMesh;
	PCAData *pcaData;
	
	vector<Vec3d> markerPositions;
	vector<int> markerIndices;

	PCAFitGF(FastTriMesh *iMesh, PCAData *iData);
	~PCAFitGF();

	void deleteMarker(int meshPos);

	void solve(int maxIter = 1000000);

	void applyDef(Vecd &variables, bool updateAll = true);
	void zeroDeformation(bool recalc = true);

	void updateCurrent(int pt);
	void addGradientVec(int index, Vec3d v);

	virtual double evaluateFunction(Vecd& variables);
	virtual void evaluateGradient(Vecd& variables, Vecd& gradient);
//	virtual void solverStep();
};



#endif