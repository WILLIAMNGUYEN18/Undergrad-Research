#include <fstream>
#include "pcaFitGF.h"
using namespace std;


static double sqr(double x) { return x*x; }


// PCAData functions ================================================

void PCAData::reconstructMesh(FastTriMesh *mesh, Vecd &weights, vector<int> *markerIndices) {
	int i, j;

	if (markerIndices == NULL) {
		memcpy(mesh->verts, average, sizeof(float)*mesh->numPts*3);
		for (j=0; j < activeComponents; j++) {
			int ofs = j * mesh->numPts * 3;
			for (i=0; i < mesh->numPts*3; i++) {
				mesh->verts[i] += weights[j] * components[ofs + i];
			}
		}
	}
	else {
		for (i=0; i < markerIndices->size(); i++) {
			int index = (*markerIndices)[i] * 3;
			mesh->verts[index+0] = average[index+0];
			mesh->verts[index+1] = average[index+1];
			mesh->verts[index+2] = average[index+2];

			for (j=0; j < activeComponents; j++) {
				int ofs = j * mesh->numPts * 3;
				mesh->verts[index+0] += weights[j] * components[ofs + index+0];
				mesh->verts[index+1] += weights[j] * components[ofs + index+1];
				mesh->verts[index+2] += weights[j] * components[ofs + index+2];
			}
		}
	}
	mesh->dirtyVerts = true;
//	mesh->calcNormals();
}


// PCAFitGF functions ===============================================

PCAFitGF::PCAFitGF(FastTriMesh *iMesh, PCAData *iData) {
	curMesh = iMesh;
	pcaData = iData;

	markerVariance = sqr(0.01);

	vars.resize(iData->activeComponents, true);
	grad.resize(vars.size());
}

PCAFitGF::~PCAFitGF() {
}

void PCAFitGF::solve(int maxIter) {
	LBFGSSolver *solver;

	solver = new LBFGSSolver(this);
	solver->solve(1e+3, 1e-5, vars, maxIter);
	delete solver;
	solver = NULL;

	applyDef(vars);
}

void PCAFitGF::applyDef(Vecd &variables, bool updateAll) {
	if (&variables != &vars)
		vars = variables;

	if (updateAll)
		pcaData->reconstructMesh(curMesh, vars);
	else
		pcaData->reconstructMesh(curMesh, vars, &markerIndices);
}

void PCAFitGF::zeroDeformation(bool recalc) {
	if (vars.size() != pcaData->activeComponents) {
		vars.resize(pcaData->activeComponents);
		grad.resize(vars.size());
	}

	vars.zeroElements();
	if (recalc) {
		pcaData->reconstructMesh(curMesh, vars);
		curMesh->calcNormals();
	}
}

void PCAFitGF::addGradientVec(int index, Vec3d v) {
	int i;
	for (i=0; i < pcaData->activeComponents; i++) {
		int ofs = (i*curMesh->numPts + index) * 3;
		grad[i] += v * Vec3d(pcaData->components[ofs+0], pcaData->components[ofs+1], pcaData->components[ofs+2]);
	}
}

double PCAFitGF::evaluateFunction(Vecd& variables) {
	double ret = 0;
	double pcaErr = 0;
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

	lastErr = ret;
	return ret;
}

void PCAFitGF::evaluateGradient(Vecd& variables, Vecd& gradient) {
	gradient = grad;
}

void PCAFitGF::deleteMarker(int meshPos){
	//Strategy: swap then pop.  Swap the position to be deleted with the last position in the vector,
	//then pop the new last item off the list to delete it.
	for(int i = 0; i < markerIndices.size(); i++){
		if(markerIndices[i] == meshPos){
			meshPos = i;
			break;
		}
	}
	Vec3d tempPos = markerPositions[meshPos];
	int tempIndex = markerIndices[meshPos];
	int lastIndex = (int)markerPositions.size() - 1;
	markerPositions[meshPos] = markerPositions[lastIndex];
	markerIndices[meshPos] = markerIndices[lastIndex];
	markerPositions.pop_back();
	markerIndices.pop_back();
}
/*void PCAFitGF::solverStep() {
//		if (stepCount % 4 == 0)
//			smoothnessStrength *= 0.5;

	if (capRate > 0 && (stepCount < 8 || stepCount % capRate == 0)) {
		cout << "frame " << stepCount << "; error " << lastErr << endl;
//			defMesh->calcNormals();

		if (recordOptimization) {
			char fname[80];
			sprintf(fname, "anim/%04d.tga", stepCount);
			redrawV();
			uiScreenshot(fname);
		}
		else {
			redrawV();
			uiWait();
		}
	}

	stepCount++;
}*/
