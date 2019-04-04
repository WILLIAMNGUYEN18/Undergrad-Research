#include "doppel2.h"
#include "skinMatch.h"
#include "quatnorm.h"
#include "markers.h"
#include "uSkin.h"
#include <float.h>
#include "vl/VLd.h"

static const double CONF_MULT = 5.0;



// SkinMatchGF =============================================

SkinMatchGF::SkinMatchGF(TriMesh *mesh, vector <int>* neigh, USkin *iSkin) : LADeformationGoalFunction(mesh, neigh) {
	int i, j;

	skin = iSkin;
	origDressPts = new Vec3d[skin->numPts];

	restrict = NULL;
	baryRecon = NULL;
	ptRots = NULL;
	tmNeigh = NULL;
	targets = NULL;
	findTargets = true;

	varsPerVert = 3;

	vars.resize(skin->numPts * varsPerVert, true);
	grad.resize(skin->numPts * varsPerVert);

	templateEdges = buildEdges(mesh, numEdges);
	restAreas = new double[skin->numPts];
	enableSurfaceMatch = true;
	enableSmoothness = true;
	enableMarkerMatch = true;
	conf = new double[curMesh->numPts()];

	if (tmNeigh)
		delete tmNeigh;
	tmNeigh = findTMNeighbors(curMesh);

	newMatch();
}


void SkinMatchGF::newMatch() {
	int i, j;

	zeroDeformation();

	for (i=0; i < skin->numPts; i++)
		origDressPts[i] = skin->dressPts[i];

/*	for (i=0; i < numEdges; i++) {
		templateEdges[i].init(mesh); //matchTM);
	}

	for (i=0; i < matchTM->numTris(); i++) {
		int v0 = matchTM->getTri(i, 0);
		int v1 = matchTM->getTri(i, 1);
		int v2 = matchTM->getTri(i, 2);

		Vec3d v = matchTM->getPt(v2)-matchTM->getPt(v0);
		Vec3d w = matchTM->getPt(v1)-matchTM->getPt(v0);
		restAreas[i] = 0.5 * (v ^ w).length();
	}
*/
	applyDef(vars);
}

void SkinMatchGF::unfold(int reps) {
}

void SkinMatchGF::applyDef(Vecd &variables) {
	int i;

	if (&variables != &vars)
		vars = variables;

	for (i=0; i < skin->numPts; i++) {
		skin->dressPts[i] = origDressPts[i] + Vec3d(vars[i*3+0],vars[i*3+1],vars[i*3+2]);
	}
	skin->updatePts();

	for (i=0; i < curMesh->numPts(); i++) {
		curMesh->getPt(i) = skin->curPts[i];
	}
	curMesh->calcNormals();

}

void SkinMatchGF::zeroDeformation() {
	memset(vars.n, 0, sizeof(double)*skin->numPts);
}

void SkinMatchGF::addGradientVec(int index, Vec3d v) {
//	skin->calcGrad(index, v, &grad[index*3], NULL, NULL, NULL, NULL);
}

double SkinMatchGF::evaluateFunction(Vecd& variables) {
	double ret = 0;
	double surfaceErr = 0;
	double smoothErr = 0;
	double markerErr = 0;
	int i, j, k;

	static double *vColors = NULL;
	static int *vCounts;
	if (vColors == NULL) {
		vColors = new double[curMesh->numPts()];
		vCounts = new int[curMesh->numPts()];
	}
	memset(vColors, 0, sizeof(double)*curMesh->numPts());
	memset(vCounts, 0, sizeof(int)*curMesh->numPts());

	uiWait();

	// update curMesh
	applyDef(variables);

	// zero out gradient
	grad.zeroElements();

	if (findTargets) {
		if (!targets)
			targets = new Vec3d[curMesh->numPts()];

		for (i = 0; i < curMesh->numPts(); i++) {
			targets[i] = Vec3d();
			conf[i] = 0;

			if (surfWeights[i] > 0) {
				updateCurrent(i);
				curSurfaceWeight *= surfWeights[i];
				conf[i] = curSurfaceWeight;
				targets[i] = curClosestPt;
			}
		}
		findTargets = false;
	}

	for (i = 0; i < curMesh->numPts(); i++) {
		if (restrict && restrict[i])
			continue;

		// surface term
		if (surfaceMatchWeight > 0 && conf[i] > 0) {
			curSurfaceWeight = baMin(1.0,CONF_MULT*conf[i]) * surfaceMatchWeight;

			Vec3d delta = curMesh->getPt(i) - targets[i];
			double dist = delta.length();
			surfaceErr += curSurfaceWeight * sqr(dist);

			addGradientVec(i, curSurfaceWeight * 2.0 * dist * delta);
		}
		curMesh->getPtColor(i)[1] = baMin(1.0,CONF_MULT*conf[i]);
	}

	// smoothness term
	if (smoothnessWeight > 0 || bendWeight > 0) {
		for (i=0; i < numEdges; i++) {
			MeshEdge &ei = templateEdges[i];
			ei.update(curMesh);

			if (restrict && (restrict[ei.v0] || restrict[ei.v1]))
				continue;

			for (j=0; j < 3; j++) {
				double delta = vars[ei.v0*3 + j] - vars[ei.v1*3 + j];
				double weight = smoothnessWeight * neighWeights[ei.v0] * neighWeights[ei.v1];

				smoothErr += weight * sqr(delta);
				grad[ei.v0*3 + j] += 2.0 * weight * delta;
				grad[ei.v1*3 + j] -= 2.0 * weight * delta;
			}
		}
	}

	if (enableSurfaceMatch) {
		ret += surfaceErr;
	}

	if (enableSmoothness) {
		ret += smoothErr;
	}

	if (enableMarkerMatch && markerMatchWeight > 0) {
		for (i=0; i < baMin(srcMarkers->numMarkers, markers->numMarkers); i++) {
			Vec3d v = srcMarkers->curPos(i);
			if (v.iszero() || markers->v(i).iszero())
				continue;

			Vec3d delta = v - markers->v(i);
			markerErr += delta.length2();

			addGradientVec(srcMarkers->markers[i].baryVerts[0], 
				srcMarkers->markers[i].baryPos[0] * markerMatchWeight * 2.0 * delta);
			if (srcMarkers->markers[i].baryVerts[1] >= 0)
				addGradientVec(srcMarkers->markers[i].baryVerts[1], 
					srcMarkers->markers[i].baryPos[1] * markerMatchWeight * 2.0 * delta);
			if (srcMarkers->markers[i].baryVerts[2] >= 0)
				addGradientVec(srcMarkers->markers[i].baryVerts[2], 
					srcMarkers->markers[i].baryPos[2] * markerMatchWeight * 2.0 * delta);
//			addGradientVec(n, markerMatchWeight * 2.0 * delta);
		}

		markerErr *= markerMatchWeight;
		ret += markerErr;
	}

	cout << ret << " (surface: " << surfaceErr << "; smoothness: " << smoothErr << "; markers: " << markerErr << ")" << endl;


	double max = 0;
	for (i=0; i < grad.size(); i++) {
		if (!_finite(grad[i]))
			grad[i] = 0;
		if (fabs(grad[i]) > max)
			max = fabs(grad[i]);
	}


	lastErr = ret;
	return ret;
}
