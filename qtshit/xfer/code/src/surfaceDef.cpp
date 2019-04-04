#include "doppel2.h"
#include <float.h>
#include "surfaceDef.h"
#include "trimesh.h"
#include "trimesh_util.h"
#include "markers.h"
#include "ba.h"
#include "quatnorm.h"
//#include "skinCalc.h"

double springCoeff = 1;
double spandex = 0.7;
//bool recordOptimization = true;
bool recordOptimization = false;

//#define SPLIT_TRANSFORMATION
#define USE_FINITE_DIFF_DERIVS

double angleComp(double a1, double &a2, bool update = false) {
	double diff = a1 - a2;
	if (fabs(a1 - (a2 + 2.0 * PI)) < diff) {
		diff = a1 - (a2 + 2.0 * PI);
		if (update)
			a2 += 2.0 * PI;
	}
	else if (fabs(a1 - (a2 - 2.0 * PI)) < diff) {
		diff = a1 - (a2 - 2.0 * PI);
		if (update)
			a2 -= 2.0 * PI;
	}
	return diff;
}

// LADeformationGoalFunction function ===============================

LADeformationGoalFunction::LADeformationGoalFunction(TriMesh *mesh, vector<int> *neigh) {
	curMesh = mesh;
	if (neigh)
		curMeshNeigh = neigh;
	else
		curMeshNeigh = findTMNeighbors(curMesh);

	origVerts = new Vec3d[curMesh->numPts()];
	int i;
	for (i=0; i < curMesh->numPts(); i++)
		origVerts[i] = curMesh->getPt(i);

	neighWeights = new double[curMesh->numPts()];
	surfWeights = new double[curMesh->numPts()];
	for (i=0; i < curMesh->numPts(); i++) {
		neighWeights[i] = 1;
		surfWeights[i] = 1;
	}

	lockShape = true;
	enableSurfaceMatch = true;
	enableSmoothness = true;
	enableMarkerMatch = true;
	smatchShowError = false;

	surfaceMatchWeight = 1.0;
	smoothnessWeight = 1.0;
	markerMatchWeight = 1.0;

	stepCount = 0;
	capRate = 1;

	vertList = NULL;
	markers = NULL;

	edgeList = NULL;

#ifdef SPLIT_TRANSFORMATION
	varsPerVert = 10;
#else
	varsPerVert = 12;
#endif

	vars.resize(curMesh->numPts() * varsPerVert, true);
	grad.resize(curMesh->numPts() * varsPerVert);
}

LADeformationGoalFunction::~LADeformationGoalFunction() {
	if (vertList)
		delete []vertList;
}

void LADeformationGoalFunction::prepareTriMesh(TriMesh *dtm) {
	targetTM = dtm;

	if (edgeList)
		delete edgeList;
	edgeList = new EdgeList();
	edgeList->buildFromTriMesh(*targetTM);
	if (vertList)
		delete []vertList;
	vertList = new char[edgeList->numVerts];
	edgeList->markVerts(vertList);
}

void LADeformationGoalFunction::applyDef(Vecd &variables) {
	int i;

	if (&variables != &vars)
		vars = variables;

	for (i=0; i < variables.size(); i += varsPerVert) {
		Vec3d cp = origVerts[i/varsPerVert];

#ifdef SPLIT_TRANSFORMATION
		// scale
		Vec3d result = prod(Vec3d(variables[i+4], variables[i+5], variables[i+6]), cp);
		// rotation
		result = QuatNorm(variables[i+0], variables[i+1], variables[i+2], variables[i+3]).toMatrixD() * result;
		// translation
		result[0] += variables[i+7];
		result[1] += variables[i+8];
		result[2] += variables[i+9];

		curMesh->getPt(i/varsPerVert) = result;
#else
		curMesh->getPt(i/varsPerVert) = Vec3d( 
			variables[i+0] * cp[0] + variables[i+1] * cp[1] + variables[i+ 2] * cp[2] + variables[i+3],
			variables[i+4] * cp[0] + variables[i+5] * cp[1] + variables[i+ 6] * cp[2] + variables[i+7],
			variables[i+8] * cp[0] + variables[i+9] * cp[1] + variables[i+10] * cp[2] + variables[i+11]);
#endif
	}
	curMesh->calcNormals();
}

void LADeformationGoalFunction::zeroDeformation() {
	if (vars.size() != curMesh->numPts() * varsPerVert) {
		vars.resize(curMesh->numPts() * varsPerVert);
		grad.resize(curMesh->numPts() * varsPerVert);
	}

	vars.zeroElements();

	int i;
	for (i=0; i < curMesh->numPts(); i++) {
#ifdef SPLIT_TRANSFORMATION
		vars[i*varsPerVert + 3] = 1;
		vars[i*varsPerVert + 4] = 1;
		vars[i*varsPerVert + 5] = 1;
		vars[i*varsPerVert + 6] = 1;
#else
		vars[i*varsPerVert + 0] = 1;
		vars[i*varsPerVert + 5] = 1;
		vars[i*varsPerVert + 10] = 1;
#endif
	}
}

const double CONF_MULT = 5.0;
const double NORMAL_TOL = cos(70 * DEG_TO_RAD);

void LADeformationGoalFunction::updateCurrent(int pt) {
	curClosestPt = Vec3d();

	curDistance = 0;
	curSurfaceWeight = 0;
	curGradient = Vec3d();

	targetTM->closestRestrictNormal = true;
//		targetTM->closestRestrictNormal = false;
	targetTM->closestNormalRestriction = curMesh->getPtNormal(pt);
//		targetTM->closestNormalRestriction.normalize();

//#define CAST_RAYS
	/*
if (pt >= 32*32*76 && pt < 32*32*77) {
//#ifdef CAST_RAYS
	if (targetTM->calcRayIntersection(curMesh->evalPts[pt].dispPos, -curMesh->evalPts[pt].dispNorm)) {
		if (targetTM->hitPos && (!targetTM->hitNeg || targetTM->tPos < fabs(targetTM->tNeg)) && targetTM->tPos < 0.1) {
//				curClosestPt = curMesh->evalPts[pt].dispPos + hitPos*curMesh->evalPts[pt].dispNorm;
			curDistance = targetTM->tPos;
			curGradient = targetTM->tPos * curMesh->evalPts[pt].dispNorm;
		}
		else if (targetTM->hitNeg && -targetTM->tNeg < 0.1) {
//				curClosestPt = curMesh->evalPts[pt].dispPos + hitNeg*curMesh->evalPts[pt].dispNorm;
			curDistance = -targetTM->tNeg;
			curGradient = targetTM->tNeg*curMesh->evalPts[pt].dispNorm;
		}
	}
}
//#else
else {*/
	if (targetTM->calcClosestPoint(curMesh->getPt(pt), 1.10)) {
		curClosestPt = targetTM->closestPt;

		Vec3d delta = targetTM->closestPt - curMesh->getPt(pt);
		delta.normalize();

		// if closest point is a vertex, check if we need to flip the sign
		if (targetTM->closestTri[1] == -1) {
//				curSurfaceWeight = 1.0 - vertList[targetTM->closestTri[0]];
			curSurfaceWeight = targetTM->getPtConf(targetTM->closestTri[0]);
//				if (curSurfaceWeight < 0 || curSurfaceWeight > 1) {
//					cout << curSurfaceWeight << ": " << targetTM->closestBary[0] << ", " << vertList[targetTM->closestTri[0]] << endl;
//				}

			double dotProd = delta * targetTM->getPtNormal(targetTM->closestTri[0]);
			if (dotProd < 0)
				targetTM->closestDist *= -1;

			// check if this vertex adjoins a hole
			if (vertList[targetTM->closestTri[0]]) {
				curSurfaceWeight = 0;
				return;
			}

			// check normal match
			if (targetTM->closestNormalRestriction * targetTM->getPtNormal(targetTM->closestTri[0]) < NORMAL_TOL) {
				curSurfaceWeight = 0;
				return;
			}
		}
		// if closest point is an edge, check if we need to flip the sign
		else if (targetTM->closestTri[2] == - 1) {
			curSurfaceWeight = targetTM->closestBary[0] * targetTM->getPtConf(targetTM->closestTri[0]) + targetTM->closestBary[1] * targetTM->getPtConf(targetTM->closestTri[1]);
//				curSurfaceWeight = 1.0 - (targetTM->closestBary[0] * vertList[targetTM->closestTri[0]] + targetTM->closestBary[1] * vertList[targetTM->closestTri[1]]);
//				if (curSurfaceWeight < 0 || curSurfaceWeight > 1) {
//					cout << curSurfaceWeight << ": " << targetTM->closestBary[0] << ", " << vertList[targetTM->closestTri[0]] << "; " << targetTM->closestBary[1] << ", " << vertList[targetTM->closestTri[1]] << endl;
//				}

			EdgeInfo *ei = edgeList->findEdge(targetTM->closestTri[0], targetTM->closestTri[1]);
			if (!ei) {
				cout << "error: unknown edge!!" << endl;
				return;
			}
			double dotProd = delta * ei->normal;
			if (dotProd < 0)
				targetTM->closestDist *= -1;

			// check if this edge adjoins a hole
			if (ei->count == 1) {
				curSurfaceWeight = 0;
				return;
			}

			// check normal match
			if (targetTM->closestNormalRestriction * ei->normal < NORMAL_TOL) {
				curSurfaceWeight = 0;
				return;
			}
		}
		else {
			curSurfaceWeight = 
				targetTM->closestBary[0] * targetTM->getPtConf(targetTM->closestTri[0]) + 
				targetTM->closestBary[1] * targetTM->getPtConf(targetTM->closestTri[1]) + 
				targetTM->closestBary[2] * targetTM->getPtConf(targetTM->closestTri[2]);
//				curSurfaceWeight = 1.0 - (targetTM->closestBary[0] * vertList[targetTM->closestTri[0]] + targetTM->closestBary[1] * vertList[targetTM->closestTri[1]] + targetTM->closestBary[2] * vertList[targetTM->closestTri[2]]);
//				if (curSurfaceWeight < 0 || curSurfaceWeight > 1) {
//					cout << curSurfaceWeight << ": " << targetTM->closestBary[0] << ", " << vertList[targetTM->closestTri[0]] << "; " << targetTM->closestBary[1] << ", " << vertList[targetTM->closestTri[1]] << "; " << targetTM->closestBary[2] << ", " << vertList[targetTM->closestTri[2]] << endl;
//				}

			Vec3d verts[3];
			verts[0] = targetTM->getPt(targetTM->closestTri[0]);
			verts[1] = targetTM->getPt(targetTM->closestTri[1]);
			verts[2] = targetTM->getPt(targetTM->closestTri[2]);
			Vec3d norm = -(verts[1] - verts[0]) ^ (verts[2] - verts[0]);
			norm.normalize();

			// check normal match
			if (targetTM->closestNormalRestriction * norm < NORMAL_TOL) {
				curSurfaceWeight = 0;
				return;
			}
		}

		curDistance = targetTM->closestDist;
		curGradient = curMesh->getPt(pt) - targetTM->closestPt;
	}
//}

//		double x = max(0, baMin(1.0, curDistances[pt] / 0.1 + 0.5));
//		defMesh->m_colors[pt] = Vec3d(x, x, x);
}

void LADeformationGoalFunction::addGradientVec(int index, Vec3d v) {
	Vec3d cp = origVerts[index];
	int ofs = index*varsPerVert;

#ifdef SPLIT_TRANSFORMATION
	int i;
	Mat4d qVal;
	Mat4d qCurDerivs[4];
	QuatNorm(vars[ofs+0], vars[ofs+1], vars[ofs+2], vars[ofs+3]).getMatrices(qVal, qCurDerivs);
	Vec3d scaleVec = Vec3d(vars[ofs+4], vars[ofs+5], vars[ofs+6]);
	Vec3d scalePt = prod(scaleVec, cp);

	grad[ofs+4 + 0] = (qVal * prod(Vec3d(1, 0, 0), cp)) * v;
	grad[ofs+4 + 1] = (qVal * prod(Vec3d(0, 1, 0), cp)) * v;
	grad[ofs+4 + 2] = (qVal * prod(Vec3d(0, 0, 1), cp)) * v;

//	for (i=0; i < 4; i++)
//		grad[ofs + i] = (qCurDerivs[i] * scalePt) * v;

	for (i=0; i < 3; i++)
		grad[ofs+7 + i] = v[i];
#else
		grad[ofs + 0] += cp[0] * v[0];
		grad[ofs + 1] += cp[1] * v[0];
		grad[ofs + 2] += cp[2] * v[0];
	grad[ofs + 3] += v[0];

		grad[ofs + 4] += cp[0] * v[1];
		grad[ofs + 5] += cp[1] * v[1];
		grad[ofs + 6] += cp[2] * v[1];
	grad[ofs + 7] += v[1];

		grad[ofs + 8] += cp[0] * v[2];
		grad[ofs + 9] += cp[1] * v[2];
		grad[ofs + 10] += cp[2] * v[2];
	grad[ofs + 11] += v[2];
#endif
}

double LADeformationGoalFunction::evaluateFunction(Vecd& variables) {
	double ret = 0;
	double surfaceErr = 0;
	double smoothErr = 0;
	double markerErr = 0;
	int i, j;

	uiWait();

	// update curMesh
	applyDef(variables);

	// zero out gradient
	grad.zeroElements();

	for (i = 0; i < curMesh->numPts(); i++) {
		// surface term
		if (enableSurfaceMatch) {
			updateCurrent(i);

			// hack
			if (lockShape)
				curSurfaceWeight *= surfWeights[i];
			curMesh->getPtColor(i)[1] = curSurfaceWeight;

			curSurfaceWeight = baMin(1.0,CONF_MULT*curSurfaceWeight);

			surfaceErr += curSurfaceWeight * sqr(curDistance);

			addGradientVec(i, surfaceMatchWeight * curSurfaceWeight * 2.0 * curGradient);

			if (smatchShowError) {
				curMesh->getPtColor(i)[0] = 0.2 + baMin(0.8, sqr(curDistance) * 4000.0);
			}
		}

		// smoothness term
		if (enableSmoothness) {
			double curSmoothnessWeight = 1 + (1.0 - curSurfaceWeight);
			for (j=0; j < curMeshNeigh[i].size(); j++) {
				int n = curMeshNeigh[i][j];

				double nW = neighWeights[i] + neighWeights[n]; //curMesh->evalPts[i].neighWeights[j];

				int k;
				double sum = 0;
				for (k=0; k < varsPerVert; k++) {
					double dist = variables[i*varsPerVert + k] - variables[n*varsPerVert + k];
					sum += curSmoothnessWeight * nW * sqr(dist);
					grad[i*varsPerVert + k] += curSmoothnessWeight * nW * smoothnessWeight * 2.0 * dist;
					grad[n*varsPerVert + k] -= curSmoothnessWeight * nW * smoothnessWeight * 2.0 * dist;
				}
				smoothErr += sum;

				if (smatchShowError) {
					curMesh->getPtColor(i)[1] = 0.2 + baMin(0.8, sum * 20000.0);
				}
			}
		}

//		if (smatchShowError) {
//			curMesh->getPtColor(i)[2] = 0.2;
//		}
/*
#ifdef SPLIT_TRANSFORMATION
		// penalize translation
		ret += sqr(variables[i*varsPerVert + 7]) + sqr(variables[i*varsPerVert + 8]) + sqr(variables[i*varsPerVert + 9]);
		grad[i*varsPerVert + 7] += 2.0 * variables[i*varsPerVert + 7];
		grad[i*varsPerVert + 8] += 2.0 * variables[i*varsPerVert + 8];
		grad[i*varsPerVert + 9] += 2.0 * variables[i*varsPerVert + 9];
#endif*/
	}

	if (enableSurfaceMatch) {
		surfaceErr *= surfaceMatchWeight;
		ret += surfaceErr;
	}

	if (enableSmoothness) {
		smoothErr *= smoothnessWeight;
		ret += smoothErr;
	}

	if (enableMarkerMatch) {
		for (i=0; i < baMin((int)markerRefs.size(), (int)markers->numMarkers); i++) {
			int n = markerRefs[i];
			if (n < 0 || markers->v(i).iszero())
				continue;

			Vec3d delta = curMesh->getPt(n) - markers->v(i);
			markerErr += delta.length2();
			addGradientVec(n, markerMatchWeight * 2.0 * delta);
		}

		markerErr *= markerMatchWeight;
		ret += markerErr;
	}

	cout << ret << " (surface: " << surfaceErr << "; smoothness: " << smoothErr << "; markers: " << markerErr << ")" << endl;

	lastErr = ret;
	return ret;
}

void LADeformationGoalFunction::evaluateGradient(Vecd& variables, Vecd& gradient) {
	gradient = grad;
}

void LADeformationGoalFunction::solverStep() {
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
}


// EdgeMatchGF ============================================

EdgeMatchGF::EdgeMatchGF(TriMesh *mesh, vector <int>* neigh, TriMesh *matchTM) : LADeformationGoalFunction(mesh, neigh) {
	int i;

	restrict = NULL;
	baryRecon = NULL;

	varsPerVert = 3;
	vars.resize(curMesh->numPts() * varsPerVert, true);
	grad.resize(curMesh->numPts() * varsPerVert);

	geodesics = new vector<TMNeigh>[curMesh->numPts()];
/*
	geodesics = findTMGeodesics(matchTM, 0.1);
	undefVerts = new Vec3d(matchTM->numPts());
	for (i=0; i < matchTM->numPts(); i++)
		undefVerts[i] = matchTM->getPt(i);
*/
	templateEdges = buildEdges(mesh, numEdges);
	for (i=0; i < numEdges; i++) {
		templateEdges[i].init(matchTM);
	}
	restAreas = new double[matchTM->numTris()];
	for (i=0; i < matchTM->numTris(); i++) {
		int v0 = curMesh->getTri(i, 0);
		int v1 = curMesh->getTri(i, 1);
		int v2 = curMesh->getTri(i, 2);

		Vec3d v = curMesh->getPt(v2)-curMesh->getPt(v0);
		Vec3d w = curMesh->getPt(v1)-curMesh->getPt(v0);
		restAreas[i] = 0.5 * (v ^ w).length();
	}
	enableSurfaceMatch = true;
	enableSmoothness = true;
	enableMarkerMatch = true;
	conf = new double[curMesh->numPts()];
}

void EdgeMatchGF::unfold(int reps) {
	int n = curMesh->numPts();
	Vec3d *newV = new Vec3d[n];
	int *newC = new int[n];
	int pt, rep, edge;
	Mat4d m;

	for (rep=0; rep < reps; rep++) {
		memset(newV, 0, sizeof(Vec3d)*n);
		memset(newC, 0, sizeof(int)*n);

		for (edge=0; edge < numEdges; edge++) {
			if (templateEdges[edge].f1 < 0)
				continue;

			templateEdges[edge].update(curMesh);

			Vec3d x0 = curMesh->getPt(templateEdges[edge].opp0);
			Vec3d x1 = curMesh->getPt(templateEdges[edge].v0);
			Vec3d x2 = curMesh->getPt(templateEdges[edge].v1);
			Vec3d x3 = curMesh->getPt(templateEdges[edge].opp1);
	
			Vec3d x20 = x0 - x2;
			m = QuatNorm(-(templateEdges[edge].eCurAngle - templateEdges[edge].eRestAngle) / 2.0,
				templateEdges[edge].e).toMatrixD();
			x20 = m * x20;

			Vec3d x23 = x3 - x2;
			m = QuatNorm((templateEdges[edge].eCurAngle - templateEdges[edge].eRestAngle) / 2.0,
				templateEdges[edge].e).toMatrixD();
			x23 = m * x23;

			Vec3d delta = 0.5 * (0.5 * (x0 - x2 + x3 - x2) - 0.5 * (x20 + x23));
//			Vec3d delta;

			newV[templateEdges[edge].opp0] += x2 + x20 + delta;
			newC[templateEdges[edge].opp0]++;
			newV[templateEdges[edge].v0] += x1 + delta;
			newC[templateEdges[edge].v0]++;
			newV[templateEdges[edge].v1] += x2 + delta;
			newC[templateEdges[edge].v1]++;
			newV[templateEdges[edge].opp1] += x2 + x23 + delta;
			newC[templateEdges[edge].opp1]++;
		}

		for (pt=0; pt < n; pt++) {
			if (newC[pt] > 0) {
				curMesh->getPt(pt) = (1.0 / newC[pt]) * newV[pt];
			}
		}
	}
	curMesh->calcNormals();

	for (pt=0; pt < n; pt++) {
		vars[pt*3+0] = curMesh->getPt(pt)[0];
		vars[pt*3+1] = curMesh->getPt(pt)[1];
		vars[pt*3+2] = curMesh->getPt(pt)[2];
	}

	delete []newV;
	delete []newC;
}

void EdgeMatchGF::applyDef(Vecd &variables) {
	int i;

	if (&variables != &vars)
		vars = variables;

	for (i=0; i < variables.size(); i += varsPerVert) {
//		Vec3d cp = origVerts[i/varsPerVert];
		curMesh->getPt(i/varsPerVert) = Vec3d(variables[i+0], variables[i+1], variables[i+2]);
	}
	curMesh->calcNormals();
}

void EdgeMatchGF::zeroDeformation() {
	if (vars.size() != curMesh->numPts() * varsPerVert) {
		vars.resize(curMesh->numPts() * varsPerVert);
		grad.resize(curMesh->numPts() * varsPerVert);
	}

	int i;
	for (i=0; i < curMesh->numPts(); i++) {
		vars[i*varsPerVert + 0] = origVerts[i][0];
		vars[i*varsPerVert + 1] = origVerts[i][1];
		vars[i*varsPerVert + 2] = origVerts[i][2];
//		vars[i*varsPerVert + 0] = origVerts[i][0] + boundedRand(-0.5, 0.5);
//		vars[i*varsPerVert + 1] = origVerts[i][1] + boundedRand(-0.5, 0.5);
//		vars[i*varsPerVert + 2] = origVerts[i][2] + boundedRand(-0.5, 0.5);
	}
}

void EdgeMatchGF::addGradientVec(int index, Vec3d v) {
//	Vec3d cp = origVerts[index];
	int ofs = index*varsPerVert;
	baAssert(varsPerVert == 3, "varsPerVert is wrong", false);

	grad[ofs + 0] += v[0];
	grad[ofs + 1] += v[1];
	grad[ofs + 2] += v[2];
}

double EdgeMatchGF::evaluateFunction(Vecd& variables) {
	double ret = 0;
	double surfaceErr = 0;
	double smoothErr = 0;
	double markerErr = 0;
	int i, j;

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

	for (i = 0; i < curMesh->numPts(); i++) {
		if (restrict && restrict[i])
			continue;

		// surface term
		if (surfaceMatchWeight > 0) {
			updateCurrent(i);

//			curMesh->getPtColor(i)[1] = curSurfaceWeight;
			conf[i] = curSurfaceWeight;
			curSurfaceWeight = baMin(1.0,CONF_MULT*curSurfaceWeight);

			// hack
//			if (lockShape)
//				curSurfaceWeight *= surfWeights[i];
//			curMesh->getPtColor(i)[1] = curSurfaceWeight;

//			curSurfaceWeight = baMin(1.0,CONF_MULT*curSurfaceWeight);

			surfaceErr += curSurfaceWeight * sqr(curDistance);

			addGradientVec(i, surfaceMatchWeight * curSurfaceWeight * 2.0 * curGradient);

//			if (smatchShowError) {
				//curMesh->getPtColor(i)[0] = 0.2 + baMin(0.8, sqr(curDistance) * 4000.0);
				curMesh->getPtColor(i)[1] = curSurfaceWeight;
//			}
		}
	}
/*
	// smoothness term
	if (enableSmoothness) {
		for (i=0; i < numEdges; i++) {
			MeshEdge &ei = templateEdges[i];
			ei.update(curMesh);

			if (restrict && (restrict[ei.v0] || restrict[ei.v1]))
				continue;

			Vec3d x1 = curMesh->getPt(ei.v0);
			Vec3d x2 = curMesh->getPt(ei.v1);

			double rest = ei.eRestLen * spandex;
			double delta = 1.0 - ei.eCurLen/rest;
			double weight = smoothnessWeight * neighWeights[ei.v0] * neighWeights[ei.v1] * 0.1;
			smoothErr += weight * sqr(delta) * rest;

			if (ei.eCurLen != 0) {
				grad[ei.v0*varsPerVert + 0] += weight * 
					-2.0*delta*(x1[0]-x2[0])/rest;
				grad[ei.v0*varsPerVert + 1] += weight * 
					-2.0*delta*(x1[1]-x2[1])/rest;
				grad[ei.v0*varsPerVert + 2] += weight * 
					-2.0*delta*(x1[2]-x2[2])/rest;
				grad[ei.v1*varsPerVert + 0] -= weight * 
					-2.0*delta*(x1[0]-x2[0])/rest;
				grad[ei.v1*varsPerVert + 1] -= weight * 
					-2.0*delta*(x1[1]-x2[1])/rest;
				grad[ei.v1*varsPerVert + 2] -= weight * 
					-2.0*delta*(x1[2]-x2[2])/rest;
			}

//			grad[ei.v1*varsPerVert + 2] -= smoothnessWeight * (newDist - ei.eRestLen) *
//				(1.0 / newDist) * 2.0 * (x1[2] - x2[2]);

			if (ei.f1 >= 0) {
				Vec3d x0 = curMesh->getPt(ei.opp0);
				Vec3d x3 = curMesh->getPt(ei.opp1);

				// bend term
				double restAngle = ei.eRestAngle;
//				restAngle = 0;
//				double aDelta = sin(ei.eCurAngle / 2.0) - restSin;
				double aDelta = angleComp(ei.eCurAngle, restAngle, true);
				double restSin = sin(restAngle / 2.0);
				int deltaSign = (aDelta > 0 ? 1 : -1);
				if (fabs(aDelta) < 0.001) continue;
//				smoothErr += bendWeight*fabs(aDelta) * ei.restBendFactor;
				double weight = bendWeight * ei.restBendFactor * neighWeights[ei.v0] * neighWeights[ei.v1];
//				weight *= surfWeights[ei.v0] * surfWeights[ei.v1];
				smoothErr += fabs(aDelta) * weight;
//				smoothErr += bendWeight*sin(ei.eCurAngle) * ei.restBendFactor;

				vColors[ei.v0] += sqr(aDelta) * ei.restBendFactor; vCounts[ei.v0]++;
				vColors[ei.v1] += sqr(aDelta) * ei.restBendFactor; vCounts[ei.v1]++;

#ifdef USE_FINITE_DIFF_DERIVS
				MeshEdge newEI = ei;
				const double DELTA = 0.000001;
				Vec3d deltas[3]={Vec3d(DELTA,0,0), Vec3d(0,DELTA,0), Vec3d(0,0,DELTA)};

				for (j=0; j < 3; j++) {
					newEI.update(curMesh->getPt(ei.opp0)+deltas[j], curMesh->getPt(ei.v0), curMesh->getPt(ei.v1), curMesh->getPt(ei.opp1));
					grad[ei.opp0*3+j] += deltaSign * (angleComp(newEI.eCurAngle, ei.eCurAngle) / DELTA) * weight;
					if (2.0*aDelta* (angleComp(newEI.eCurAngle, ei.eCurAngle) / DELTA) * weight > 1e6) {
						cout << "big deriv on pt " << i << "; newAngle = " << newEI.eCurAngle << "; old angle = " << ei.eCurAngle << endl;
					}
					newEI.update(curMesh->getPt(ei.opp0), curMesh->getPt(ei.v0)+deltas[j], curMesh->getPt(ei.v1), curMesh->getPt(ei.opp1));
					grad[ei.v0*3+j] += deltaSign * (angleComp(newEI.eCurAngle, ei.eCurAngle) / DELTA) * weight;
					if (2.0*aDelta* (angleComp(newEI.eCurAngle, ei.eCurAngle) / DELTA) * weight > 1e6) {
						cout << "big deriv on pt " << i << "; newAngle = " << newEI.eCurAngle << "; old angle = " << ei.eCurAngle << endl;
					}
					newEI.update(curMesh->getPt(ei.opp0), curMesh->getPt(ei.v0), curMesh->getPt(ei.v1)+deltas[j], curMesh->getPt(ei.opp1));
					grad[ei.v1*3+j] += deltaSign * (angleComp(newEI.eCurAngle, ei.eCurAngle) / DELTA) * weight;
					if (2.0*aDelta* (angleComp(newEI.eCurAngle, ei.eCurAngle) / DELTA) * weight > 1e6) {
						cout << "big deriv on pt " << i << "; newAngle = " << newEI.eCurAngle << "; old angle = " << ei.eCurAngle << endl;
					}
					newEI.update(curMesh->getPt(ei.opp0), curMesh->getPt(ei.v0), curMesh->getPt(ei.v1), curMesh->getPt(ei.opp1)+deltas[j]);
					grad[ei.opp1*3+j] += deltaSign * (angleComp(newEI.eCurAngle, ei.eCurAngle) / DELTA) * weight;
					if (2.0*aDelta* (angleComp(newEI.eCurAngle, ei.eCurAngle) / DELTA) * weight > 1e6) {
						cout << "big deriv on pt " << i << "; newAngle = " << newEI.eCurAngle << "; old angle = " << ei.eCurAngle << endl;
					}
				}

#else
				Vec3d nAGrad, nBGrad, eGrad;
				double cosGrad, sinGrad, eLenGrad, sin2Grad;
				double cosTheta = cos(ei.eCurAngle);
				double sinTheta = sin(ei.eCurAngle);
				double sin2Theta = sin(ei.eCurAngle/2.0);
	//			smoothErr += ei.eCurAngle;

				Vec3d qA[4], qB[4];
				qA[0] = x2-x1; qA[1] = x0-x2; qA[2] = x1-x0; qA[3] = 0;
				qB[0] = 0;     qB[1] = x2-x3; qB[2] = x3-x1; qB[3] = x1-x2;
				double e[4];
				e[0] = 0; e[1] = 1; e[2] = -1; e[3] = 0;

				Mat3d m0, m1, me;
				//int foo, bar;
				//for (foo=0; foo < 3; foo++)
				//	for (bar=0; bar < 3; bar++) {
				//		m0[foo][bar] -= ei.n0[foo] * ei.n0[bar];
				//		m1[foo][bar] -= ei.n1[foo] * ei.n1[bar];
				//		me[foo][bar] -= ei.e[foo] * ei.e[bar];
				//	}

//				double mult = bendWeight * ei.restBendFactor;
				double mult = weight * 2.0 * deltaSign; //aDelta
//				double mult = bendWeight * ei.restBendFactor * deltaSign;

				for (j=0; j < 4; j++) {
					int index;
					if (j == 0)
						index = ei.opp0*varsPerVert;
					else if (j == 1)
						index = ei.v0*varsPerVert;
					else if (j == 2)
						index = ei.v1*varsPerVert;
					else
						index = ei.opp1*varsPerVert;


					//if (index/varsPerVert == 21) {
					//	cout << "--21-- n0=" << ei.n0 << " " << ei.n0Len << " n1=" << ei.n1 << " " << ei.n1Len << " e=" << ei.e << " " << ei.eCurLen << endl;
					//	cout << ei.eCurAngle << " " << restAngle << endl;
					//}

					double datan;
					const double DATAN_CAP = 1e20; //100;

					int sign = ((ei.n0 ^ ei.n1) * ei.e > 0) ? 1 : -1;

					// x
					nAGrad = m0 * ((1.0/ei.n0Len) * Vec3d(0, -qA[j][2], qA[j][1]));
					nBGrad = m1 * ((1.0/ei.n1Len) * Vec3d(0, -qB[j][2], qB[j][1]));
					eGrad = me * ((1.0/ei.eCurLen) * Vec3d(e[j], 0, 0));
					eLenGrad = 2.0 * ei.e[0]*ei.eCurLen*e[j];
					cosGrad = nAGrad * ei.n1 + ei.n0 * nBGrad;
					sinGrad = ((nAGrad ^ ei.n1) + (ei.n0 ^ nBGrad)) * ei.e + (ei.n0 ^ ei.n1) * eGrad;
					sin2Grad = -sign*0.5*(1.0/sqrt(0.5*(1-ei.n0*ei.n1)))*0.5*(ei.n0*nBGrad + nAGrad*ei.n1);
	//				grad[index++] += smoothnessWeight * (cosTheta*sinGrad - sinTheta*cosGrad) / sqr(cosTheta);
					datan = (1.0 / (1.0 + sqr(tan(ei.eCurAngle)))) * (cosTheta*sinGrad - sinTheta*cosGrad) / sqr(cosTheta);
					if (!_finite(datan)) datan = 0;
					if (fabs(datan) > DATAN_CAP) datan = DATAN_CAP*DATAN_CAP/datan;
//					datan = (cosTheta*sinGrad - sinTheta*cosGrad);
//					grad[index++] += smoothnessWeight * 2.0 * ei.eCurAngle * datan;
					grad[index++] += mult * datan;
//					grad[index++] += mult * sin2Grad;
					// y
					nAGrad = m0 * ((1.0/ei.n0Len) * Vec3d(qA[j][2], 0, -qA[j][0]));
					nBGrad = m1 * ((1.0/ei.n1Len) * Vec3d(qB[j][2], 0, -qB[j][0]));
					eGrad = me * ((1.0/ei.eCurLen) * Vec3d(0, e[j], 0));
					eLenGrad = 2.0 * ei.e[1]*ei.eCurLen*e[j];
					cosGrad = nAGrad * ei.n1 + ei.n0 * nBGrad;
					sinGrad = ((nAGrad ^ ei.n1) + (ei.n0 ^ nBGrad)) * ei.e + (ei.n0 ^ ei.n1) * eGrad;
					sin2Grad = -sign*0.5*(1.0/sqrt(0.5*(1-ei.n0*ei.n1)))*0.5*(ei.n0*nBGrad + nAGrad*ei.n1);
					datan = (1.0 / (1.0 + sqr(tan(ei.eCurAngle)))) * (cosTheta*sinGrad - sinTheta*cosGrad) / sqr(cosTheta);
					if (!_finite(datan)) datan = 0;
					if (fabs(datan) > DATAN_CAP) datan = DATAN_CAP*DATAN_CAP/datan;
					grad[index++] += mult * datan;
					// z
					nAGrad = m0 * ((1.0/ei.n0Len) * Vec3d(-qA[j][1], qA[j][0], 0));
					nBGrad = m1 * ((1.0/ei.n1Len) * Vec3d(-qB[j][1], qB[j][0], 0));
					eGrad = me * ((1.0/ei.eCurLen) * Vec3d(0, 0, e[j]));
					eLenGrad = 2.0 * ei.e[2]*ei.eCurLen*e[j];
					cosGrad = nAGrad * ei.n1 + ei.n0 * nBGrad;
					sinGrad = ((nAGrad ^ ei.n1) + (ei.n0 ^ nBGrad)) * ei.e + (ei.n0 ^ ei.n1) * eGrad;
					sin2Grad = -sign*0.5*(1.0/sqrt(0.5*(1-ei.n0*ei.n1)))*0.5*(ei.n0*nBGrad + nAGrad*ei.n1);
					datan = (1.0 / (1.0 + sqr(tan(ei.eCurAngle)))) * (cosTheta*sinGrad - sinTheta*cosGrad) / sqr(cosTheta);
					if (!_finite(datan)) datan = 0;
					if (fabs(datan) > DATAN_CAP) datan = DATAN_CAP*DATAN_CAP/datan;
					grad[index++] += mult * datan;
				}
#endif
			}
		}

		// area term
		for (i=0; i < curMesh->numTris(); i++) {
			int v0 = curMesh->getTri(i, 0);
			int v1 = curMesh->getTri(i, 1);
			int v2 = curMesh->getTri(i, 2);
			if (restrict && (restrict[v0] || restrict[v1] || restrict[v2]))
				continue;
			Vec3d p0 = curMesh->getPt(v0);
			Vec3d p1 = curMesh->getPt(v1);
			Vec3d p2 = curMesh->getPt(v2);

			Vec3d v = p2-p0;
			Vec3d w = p1-p0;
			Vec3d cross = v^w;
			double area = 0.5 * cross.length();

			double rArea = restAreas[i] * spandex * spandex;
			double delta = (1.0 - area / rArea);
			smoothErr += smoothnessWeight * sqr(delta) * rArea;

			double mult = -smoothnessWeight * delta / cross.length();
			Vec3d q = p2-p1;
			grad[v0*varsPerVert + 0] += mult * (cross[1]*-q[2] + cross[2]*q[1]);
			grad[v0*varsPerVert + 1] += mult * (cross[0]*q[2] + cross[2]*-q[0]);
			grad[v0*varsPerVert + 2] += mult * (cross[0]*-q[1] + cross[1]*q[0]);
			q = p0-p2;
			grad[v1*varsPerVert + 0] += mult * (cross[1]*-q[2] + cross[2]*q[1]);
			grad[v1*varsPerVert + 1] += mult * (cross[0]*q[2] + cross[2]*-q[0]);
			grad[v1*varsPerVert + 2] += mult * (cross[0]*-q[1] + cross[1]*q[0]);
			q = p1-p0;
			grad[v2*varsPerVert + 0] += mult * (cross[1]*-q[2] + cross[2]*q[1]);
			grad[v2*varsPerVert + 1] += mult * (cross[0]*q[2] + cross[2]*-q[0]);
			grad[v2*varsPerVert + 2] += mult * (cross[0]*-q[1] + cross[1]*q[0]);
		}
	}
*/
	
	// geodesic-based springs
	for (i=0; i < curMesh->numPts(); i++) {
		for (j=0; j < geodesics[i].size(); j++) {
			Vec3d x1 = curMesh->getPt(i);
			Vec3d x2 = curMesh->getPt(geodesics[i][j].vert);
//			double restLen = (origVerts[i] - origVerts[geodesics[i][j].vert]).length();
			double restLen = geodesics[i][j].dist;
			double newLen = (x1 - x2).length();
			double delta = 1.0 - newLen/restLen;

//			double weight = ((0.10 - geodesics[i][j].dist)/0.10) * smoothnessWeight * neighWeights[i] * neighWeights[geodesics[i][j].vert];
			double weight = smoothnessWeight * geodesics[i][j].weight * neighWeights[i] * neighWeights[geodesics[i][j].vert];
			smoothErr += weight * sqr(delta) * restLen;

			if (newLen != 0) {
				grad[i*varsPerVert + 0] += weight * 
					-2.0*delta*(x1[0]-x2[0])/newLen;
				grad[i*varsPerVert + 1] += weight * 
					-2.0*delta*(x1[1]-x2[1])/newLen;
				grad[i*varsPerVert + 2] += weight * 
					-2.0*delta*(x1[2]-x2[2])/newLen;
				grad[geodesics[i][j].vert*varsPerVert + 0] -= weight * 
					-2.0*delta*(x1[0]-x2[0])/newLen;
				grad[geodesics[i][j].vert*varsPerVert + 1] -= weight * 
					-2.0*delta*(x1[1]-x2[1])/newLen;
				grad[geodesics[i][j].vert*varsPerVert + 2] -= weight * 
					-2.0*delta*(x1[2]-x2[2])/newLen;
			}
		}
	}
	

	if (enableSurfaceMatch) {
		surfaceErr *= surfaceMatchWeight;
		ret += surfaceErr;
	}

	if (enableSmoothness) {
		ret += smoothErr;
	}

	if (enableMarkerMatch && markerMatchWeight > 0) {
		for (i=0; i < baMin((int)markerRefs.size(), (int)markers->numMarkers); i++) {
			int n = markerRefs[i];
			if (n < 0 || markers->v(i).iszero())
				continue;
			if (restrict && restrict[n])
				continue;

			Vec3d delta = curMesh->getPt(n) - markers->v(i);
			markerErr += delta.length2();
			addGradientVec(n, markerMatchWeight * 2.0 * delta);
		}

		markerErr *= markerMatchWeight;
		ret += markerErr;
	}

	cout << ret << " (surface: " << surfaceErr << "; smoothness: " << smoothErr << "; markers: " << markerErr << ")" << endl;
/*
	for (i=0; i < curMesh->numPts(); i++) {
		if (vCounts[i] == 0) {
			curMesh->getPtColor(i) = Vec3d(0,0,0);
		}
		else {
			double d = max(0, 1.0 - 100* vColors[i] / vCounts[i]);
			curMesh->getPtColor(i) = Vec3d(1, d, d);
		}
	}*/


	double max = 0;
	for (i=0; i < grad.size(); i++) {
		if (!_finite(grad[i]))
			grad[i] = 0;
		if (fabs(grad[i]) > max)
			max = fabs(grad[i]);
	}
/*	for (i=0; i < grad.size(); i+=3) {
		curMesh->getPtColor(i/3) = Vec3d(1,1,1) - Vec3d(fabs(grad[i+0]), fabs(grad[i+1]), fabs(grad[i+2])) / max;
	}*/


	lastErr = ret;
	return ret;
}


#ifdef OLD_THING
// SkinMatchGF ============================================

SkinMatchGF::SkinMatchGF(int numEx, TriMesh *mesh, Skin *sk, Skeleton *mp) : LADeformationGoalFunction(mesh, NULL) {
	int i, index;

	numExamples = numEx;
	skin = sk;
	poses = mp;

	numVars = 0;
	varIndex = new int[skin->numPts];
	for (i=0; i < skin->numPts; i++) {
		varIndex[i] = numVars;
		numVars += skin->points[i].numVars;
	}
	cout << "number of variables: " << numVars << endl;

	vars.resize(numVars);
	grad.resize(numVars);

	index = 0;
	for (i=0; i < skin->numPts; i++) {
		skin->points[i].copyToVars(vars.n + index);
		index += skin->points[i].numVars;
	}

	templateEdgeList.buildFromTriMesh(*curMesh);
}


void SkinMatchGF::applyDef(Vecd &variables) {
	if (&variables != &vars)
		vars = variables;

	int index, i;

	index = 0;
	for (i=0; i < skin->numPts; i++) {
		skin->points[i].copyFromVars(vars.n + index);
		index += skin->points[i].numVars;
	}
}

void SkinMatchGF::zeroDeformation() {
	int index, i;

	index = 0;
	for (i=0; i < skin->numPts; i++) {
		skin->points[i].copyFromVars(vars.n + index);
		index += skin->points[i].numVars;
	}
}

void SkinMatchGF::prepareTriMesh(TriMesh *dtm) {
	meshArray = dtm;

	targetTM = dtm;
/*
	neighArray = new vector<int>*[numExamples];
	neighArray[0] = curMeshNeigh;
	for (i=1; i < numExamples; i++) {
		neighArray[1] = findTMNeighbors(&meshArray[i]);
	}*/

	edgeListArray = new EdgeList*[numExamples];
	vertListArray = new char*[numExamples];

	int i;
	for (i=0; i < numExamples; i++) {
		edgeListArray[i] = new EdgeList();
		edgeListArray[i]->buildFromTriMesh(meshArray[i]);
		vertListArray[i] = new char[edgeListArray[i]->numVerts];
		edgeListArray[i]->markVerts(vertListArray[i]);
	}
	edgeList = edgeListArray[0];
	vertList = vertListArray[0];
}

void SkinMatchGF::addGradientVec(int index, Vec3d v) {
/*	Vec3d cp = origVerts[index];
	int ofs = index*varsPerVert;

	grad[ofs + 1] += v[0];
	grad[ofs + 2] += v[1];
	grad[ofs + 3] += v[2];*/
}

double SkinMatchGF::evaluateFunction(Vecd& variables) {
	double ret = 0;
	double surfaceErr = 0;
	double smoothErr = 0;
	double markerErr = 0;
	int i, j, k,ex, pt, inf;

	uiWait();

	// update curMesh
	applyDef(variables);

	// zero out gradient
	grad.zeroElements();

	for (ex = 0; ex < numExamples; ex++) {
		skin->skel->copyVals(&poses[ex]);
		skin->skel->updateCoords();
		skin->updatePoints();
		skin->updateMesh(curMesh);

		targetTM = meshArray + ex;
		vertList = vertListArray[ex];
		edgeList = edgeListArray[ex];

		for (i = 0; i < curMesh->numPts(); i++) {
			// surface term
			if (enableSurfaceMatch  && surfaceMatchWeight > 0) {
				updateCurrent(i);

				// hack
				if (lockShape)
					curSurfaceWeight *= surfWeights[i];
				curSurfaceWeight *= surfaceMatchWeight;
	//			curMesh->getPtColor(i)[1] = curSurfaceWeight;

	//			curSurfaceWeight = baMin(1.0,CONF_MULT*curSurfaceWeight);

				if (curSurfaceWeight > 0)
					surfaceErr += skin->points[i].calcGrad(grad.n + varIndex[i], targetTM->closestPt, skin->skel, curSurfaceWeight);

				if (smatchShowError) {
					curMesh->getPtColor(i)[0] = 0.2 + baMin(0.8, sqr(curDistance) * 4000.0);
				}
			}

			// smoothness term
			if (enableSmoothness && smoothnessWeight > 0) {
				for (j=0; j < templateEdgeList.edges[i].size(); j++)  {
					EdgeInfo &ei = templateEdgeList.edges[i][j];
					double oldDist = (edgeMatchTM->getPt(i) - edgeMatchTM->getPt(ei.vert)).length();

					Vec3d vi = curMesh->getPt(i);
					Vec3d vei = curMesh->getPt(ei.vert);
					Vec3d newVec = vi - vei;
					double newDist = newVec.length();

// commented out because gcc wasn't happy with these lines...
//					smoothErr += skin->points[i].calcGrad(grad.n + varIndex[i], vi + (0.5*(newDist-oldDist)/newDist)*(vei-vi), skin->skel, smoothnessWeight);
//					smoothErr += skin->points[ei.vert].calcGrad(grad.n + varIndex[ei.vert], vi + (1.0 - 0.5*(newDist-oldDist)/newDist)*(vei-vi), skin->skel, smoothnessWeight);
				}
			}
		}

/*		if (enableMarkerMatch) {
			for (i=0; i < baMin((int)markerRefs.size(), (int)markers->numMarkers); i++) {
				int n = markerRefs[i];
				if (n < 0 || markers->v(i).iszero())
					continue;

//				Vec3d delta = curMesh->getPt(n) - markers->v(i);
				markerErr += skin->points[n].calcGrad(grad.n + varIndex[n], markers->v(i), skin->skel, markerMatchWeight);
			}

			ret += markerErr;
		}*/
	}
/*
	// smoothness term
	if (enableSmoothness && smoothnessWeight > 0) {
		for (pt=0; pt < skin->numPts; pt++) {
			for (inf=0; inf < skin->points[pt].numTInf; inf++) {
				int curTrans = skin->points[pt].tTransforms[inf];
				int ptInd = varIndex[pt] + 4*inf;

				for (i=0; i < curMeshNeigh[pt].size(); i++) {
					int n = curMeshNeigh[pt][i];

					for (j=0; j < skin->points[n].numTInf; j++) {
						if (skin->points[n].tTransforms[j] == curTrans) {
							break;
						}
					}
					if (j < skin->points[n].numTInf) {
						int nInd = varIndex[n] + 4*j;

						for (k=0; k < 4; k++) {
							double delta = vars[nInd+k] - vars[ptInd+k];
							smoothErr += smoothnessWeight * delta * delta;
							grad[nInd+k] += 2.0 * smoothnessWeight * delta;
							grad[ptInd+k] -= 2.0 * smoothnessWeight * delta;
						}
					}
					else {
						for (k=0; k < 4; k++) {
							double delta = vars[ptInd+k];
							smoothErr += smoothnessWeight * delta * delta;
							grad[ptInd+k] += 2.0 * smoothnessWeight * delta;
						}
					}
				}
			}
		}
	}
*/

	if (enableSurfaceMatch) {
		ret += surfaceErr;
	}

	if (enableSmoothness) {
		smoothErr *= smoothnessWeight;
		ret += smoothErr;
	}

	cout << ret << " (surface: " << surfaceErr << "; smoothness: " << smoothErr << "; markers: " << markerErr << ")" << endl;

	lastErr = ret;
	return ret;
}
#endif

// MeshEdge ===============================================

MeshEdge::MeshEdge() {
	v0 = 0;
	v1 = 0;
	f0 = -1;
	f1 = -1;
}

void MeshEdge::update(TriMesh *tm) {
	e = tm->getPt(v0) - tm->getPt(v1);
	eCurLen = e.length();
	e = (1.0 / eCurLen) * e;

	if (f1 >= 0) {
		n0 = (tm->getPt(v1) - tm->getPt(opp0)) ^ (tm->getPt(v0) - tm->getPt(opp0));
		n0Len = n0.length();
		n0 = (1.0/n0Len) * n0;

		n1 = (tm->getPt(v0) - tm->getPt(opp1)) ^ (tm->getPt(v1) - tm->getPt(opp1));
		n1Len = n1.length();
		n1 = (1.0/n1Len) * n1;

		eCurAngle = atan2((n0 ^ n1) * e, n0 * n1);
		// handle NaN/infinite cases:
		if (!_finite(eCurAngle))
			eCurAngle = 0;

		he0 = n0Len / eCurLen;
		he1 = n1Len / eCurLen;
	}
}

void MeshEdge::update(Vec3d &x0, Vec3d &x1, Vec3d &x2, Vec3d &x3) {
	e = x1 - x2;
	eCurLen = e.length();
	e = (1.0 / eCurLen) * e;

	if (f1 >= 0) {
		n0 = (x2 - x0) ^ (x1 - x0);
		n0Len = n0.length();
		n0 = (1.0/n0Len) * n0;

		n1 = (x1 - x3) ^ (x2 - x3);
		n1Len = n1.length();
		n1 = (1.0/n1Len) * n1;

		eCurAngle = atan2((n0 ^ n1) * e, n0 * n1);
		// handle NaN/infinite cases:
		if (!_finite(eCurAngle))
			eCurAngle = 0;

		he0 = n0Len / eCurLen;
		he1 = n1Len / eCurLen;
	}
}

void MeshEdge::init(TriMesh *tm) {
	update(tm);
	eRestLen = eCurLen;
	eRestAngle = eCurAngle;
	restBendFactor = 6.0 * sqr(eCurLen) / (n0Len + n1Len);
}

MeshEdge *buildEdges(TriMesh *tm, int &numEdges) {
	int i, j, k;
	MeshEdge *edges;

	// build edge list
	vector<vector<MeshEdge> > tempEdges;
	numEdges = 0;

	tempEdges.resize(tm->numPts());
	for (i=0; i < tm->numTris(); i++) {
		for (j=0; j < 3; j++) {
			int v0 = tm->getTri(i, j);
			int v1 = tm->getTri(i, (j+1)%3);
			int opp = tm->getTri(i, (j+2)%3);

			if (v1 < v0)
				swap(v0, v1);

			int n = (int)tempEdges[v0].size();
			for (k=0; k < n; k++) {
				if (tempEdges[v0][k].v1 == v1) {
					tempEdges[v0][k].f1 = i;
					tempEdges[v0][k].opp1 = opp;
					break;
				}
			}
			if (k == n) {
				tempEdges[v0].push_back(MeshEdge());
				tempEdges[v0][n].v0 = v0;
				tempEdges[v0][n].v1 = v1;
				tempEdges[v0][n].f0 = i;
				tempEdges[v0][n].opp0 = opp;
				numEdges++;
			}
		}
	}

	edges = new MeshEdge[numEdges];
	int index = 0;
	for (i=0; i < tm->numPts(); i++) {
		for (j=0; j < tempEdges[i].size(); j++) {
//			if (tempEdges[i][j].f1 < 0) {
				// oops, there's a cut... ignore this edge
//				cout << "bogus edge: " << tempEdges[i][j].opp1 << endl;
//				numEdges--;
//			}
//			else
				edges[index++] = tempEdges[i][j];
		}
	}

	return edges;
}