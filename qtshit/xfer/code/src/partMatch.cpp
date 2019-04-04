#include <fstream>
#include "partMatch.h"
#include "markers.h"
#include "trimesh.h"
using namespace std;


static double sqr(double x) { return x*x; }


// PartMatchGF functions ============================================

PartMatchGF::PartMatchGF(PCAPart *iPart, TriMesh *iTargetMesh, MarkerSet *iTargetMkr) {
	part = iPart;
	targetMesh = iTargetMesh;
	targetMkr = iTargetMkr;

	markerVariance = sqr(0.01);

//	vars.resize(iData->activeComponents, true);
//	grad.resize(vars.size());
}

PartMatchGF::~PartMatchGF() {
}

void PartMatchGF::applyDef(Vecd &variables, bool updateAll) {
/*	if (&variables != &vars)
		vars = variables;

	if (updateAll)
		pcaData->reconstructMesh(curMesh, vars);
	else
		pcaData->reconstructMesh(curMesh, vars, &markerIndices);*/
}

void PartMatchGF::zeroDeformation(bool recalc) {
/*	if (vars.size() != pcaData->activeComponents) {
		vars.resize(pcaData->activeComponents);
		grad.resize(vars.size());
	}

	vars.zeroElements();
	if (recalc) {
		pcaData->reconstructMesh(curMesh, vars);
		curMesh->calcNormals();
	}*/
}

void PartMatchGF::addGradientVec(int index, Vec3d v) {
/*	int i;
	for (i=0; i < pcaData->activeComponents; i++) {
		int ofs = (i*curMesh->numPts + index) * 3;
		grad[i] += v * Vec3d(pcaData->components[ofs+0], pcaData->components[ofs+1], pcaData->components[ofs+2]);
	}*/
}

double PartMatchGF::evaluateFunction(Vecd& variables) {
	double ret = 0;
/*	double pcaErr = 0;
	double markerErr = 0;
	int i;

	double markerDev = 0.01;

	// reconstruct the current shape
	applyDef(variables, false);

	// zero out gradient
	grad.zeroElements();

	// pca log likelihood term
	pcaErr = 0;
	for (i=0; i < pcaData->activeComponents; i++) {	
		pcaErr += sqr(vars[i]) / pcaData->sigma2[i];
		grad[i] += 2.0 * vars[i] / pcaData->sigma2[i];
	}
	ret += pcaErr;

	// marker-matching term
	for (i=0; i < markerPositions.size(); i++) {
		int n = markerIndices[i];

		Vec3d delta = Vec3d(curMesh->verts[n*3+0], curMesh->verts[n*3+1], curMesh->verts[n*3+2]) - markerPositions[i];
		markerErr += delta.length2() / markerVariance;
		addGradientVec(n, 2.0 * delta / markerVariance);
	}
	ret += markerErr;

//	cout << ret << " (pca: " << pcaErr << "; markers: " << markerErr << ")" << endl;

	lastErr = ret;*/
	return ret;
}

void PartMatchGF::evaluateGradient(Vecd& variables, Vecd& gradient) {
	gradient = grad;
}
