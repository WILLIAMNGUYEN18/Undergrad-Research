#include "doppel2.h"
#include "rotMatch.h"
#include "quatnorm.h"
#include "markers.h"
#include <float.h>
#include "vl/VLd.h"

static const double CONF_MULT = 5.0;

//#define IMPLICIT_ROT

Mat3d ptsRot(Vec3d *v0, Vec3d *v1, int numPts) {
	int i, j, m;
	Vec3d mu0, mu1;
	VLMatd k(3,3), u(3,3), v(3,3);
	VLVecd diag(3);

	// compute centroids
	for (i=0; i < numPts; i++) {
		mu0 += v0[i];
		mu1 += v1[i];
	}
	mu0 /= numPts;
	mu1 /= numPts;

	// compute k
	k.MakeZero();
	for (i=0; i < numPts; i++) {
		for (j=0; j < 3; j++)
			for (m=0; m < 3; m++)
				k[j][m] += (v0[i][j] - mu0[j]) * (v1[i][m] - mu1[m]);
	}

	// svd
	SVDFactorization(k, u, v, diag);
	k = v * trans(u);


	double det = 
		k[0][0]*(k[1][1]*k[2][2]-k[1][2]*k[2][1]) +
		k[1][0]*(k[2][1]*k[0][2]-k[2][2]*k[0][1]) +
		k[2][0]*(k[0][1]*k[1][2]-k[0][2]*k[1][1]);
	if (det < 0) {
		double minV = diag[0];
		int minInd = 0;
		for (i=1; i < 3; i++) {
			if (diag[i] < minV) {
				minV = diag[i];
				minInd = i;
			}
		}
		v[0][minInd] *= -1;
		v[1][minInd] *= -1;
		v[2][minInd] *= -1;

		k = v * trans(u);
	}

	return Mat3d(
		k[0][0], k[0][1], k[0][2], 
		k[1][0], k[1][1], k[1][2], 
		k[2][0], k[2][1], k[2][2]);
}


// RotMatchGF =============================================

RotMatchGF::RotMatchGF(TriMesh *mesh, vector <int>* neigh, TriMesh *matchTM) : LADeformationGoalFunction(mesh, neigh) {
	int i, j;

	matchVerts = new Vec3d[matchTM->numPts()];
	for (i=0; i < matchTM->numPts(); i++)
		matchVerts[i] = matchTM->getPt(i);

	restrict = NULL;
	baryRecon = NULL;
	ptRots = NULL;
	tmNeigh = NULL;
	targets = NULL;
	findTargets = true;

#ifdef IMPLICIT_ROT
	varsPerVert = 3;
#else
	varsPerVert = 7;
#endif
	vars.resize(curMesh->numPts() * varsPerVert, true);
	grad.resize(curMesh->numPts() * varsPerVert);

//	curMesh->calcNormals();

	geodesics = new vector<TMNeigh>[curMesh->numPts()];
/*
	geodesics = findTMGeodesics(matchTM, 0.1);
	undefVerts = new Vec3d(matchTM->numPts());
	for (i=0; i < matchTM->numPts(); i++)
		undefVerts[i] = matchTM->getPt(i);
*/
	templateEdges = buildEdges(mesh, numEdges);
	restAreas = new double[matchTM->numTris()];
	enableSurfaceMatch = true;
	enableSmoothness = true;
	enableMarkerMatch = true;
	conf = new double[curMesh->numPts()];

	newMatch(curMesh, matchTM);
}


const int MAX_NEIGH = 20;

void RotMatchGF::newMatch(TriMesh *mesh, TriMesh *matchTM) {
	int i, j;

	/*
	// test: use skinned mesh as base
	for (i=0; i < mesh->numPts(); i++) {
		matchVerts[i] = mesh->getPt(i);
	}*/

	for (i=0; i < curMesh->numPts(); i++) {
		vars[i*varsPerVert + 0] = 0;
		vars[i*varsPerVert + 1] = 0;
		vars[i*varsPerVert + 2] = 0;
		vars[i*varsPerVert + 3] = 0;
		vars[i*varsPerVert + 4] = 0;
		vars[i*varsPerVert + 5] = 0;
		vars[i*varsPerVert + 6] = 1;
	}

//#ifdef OLD_WAY
#ifndef IMPLICIT_ROT
	// initialize transformations
	if (tmNeigh)
		delete tmNeigh;
	tmNeigh = findTMNeighbors(curMesh);
/*	static vector<TMNeigh> *geodesics = NULL;
	if (!geodesics) {
		geodesics = findTMGeodesics(curMesh, 0.15);
	}*/
	
	for (i=0; i < matchTM->numPts(); i++) {
		int numInvolved = 0;
		Vec3d v0[MAX_NEIGH], v1[MAX_NEIGH];
		int vInd[MAX_NEIGH];
		vInd[numInvolved++] = i;
		for (j=0; j < (int)tmNeigh[i].size(); j++) {
			vInd[numInvolved++] = tmNeigh[i][j];
			if (numInvolved >= MAX_NEIGH)
				break;
		}
		for (j=1; j < numInvolved; j++) {
			int k;
			for (k=0; k < (int)tmNeigh[vInd[j]].size(); k++) {
				int neigh = tmNeigh[vInd[j]][k];
				int x;
				for (x=0; x < numInvolved; x++)
					if (vInd[x] == neigh) break;
				if (x == numInvolved)
					vInd[numInvolved++] = neigh;
				if (numInvolved >= MAX_NEIGH)
					break;
			}
			if (numInvolved >= MAX_NEIGH)
				break;
		}


		for (j=0; j < numInvolved; j++) {
			v0[j] = matchVerts[vInd[j]];
			v1[j] = curMesh->getPt(vInd[j]);
		}
/*

		int numInvolved = min(MAX_NEIGH, (int)tmNeigh[i].size() + 1);
//		int numInvolved = min(MAX_NEIGH, (int)geodesics[i].size() + 1);
		Vec3d v0[MAX_NEIGH], v1[MAX_NEIGH];

		v0[0] = matchVerts[i];
		v1[0] = curMesh->getPt(i);
		for (j=1; j < numInvolved; j++) {
			v0[j] = matchVerts[tmNeigh[i][j-1]];
			v1[j] = curMesh->getPt(tmNeigh[i][j-1]);
//			v0[j] = matchVerts[geodesics[i][j-1].vert];
//			v1[j] = curMesh->getPt(geodesics[i][j-1].vert);
		}
*/
		Mat3d m = ptsRot(v0, v1, numInvolved);
		QuatNorm q = matToQuat(m);

		Vec3d rotPt = q.toMatrixD() * v0[0];
/*		Vec3d check = m *v0[0];
		if ((check-rotPt).length() > 0.1) {
			cout << "bad conversion: (" << (check-rotPt).length() << ")" << endl << m << endl << q << endl << q.toMatrixD() << endl;
			Vec3d c0(m[0][0], m[0][1], m[0][2]);
			Vec3d c1(m[1][0], m[1][1], m[1][2]);
			Vec3d c2(m[2][0], m[2][1], m[2][2]);
			cout << c0*c0 << " " << c1*c1 << " " << c2*c2 << endl;
			cout << c0*c1 << " " << c0*c2 << " " << c1*c2 << endl;
		}*/
		vars[i*varsPerVert + 0] = v1[0][0] - rotPt[0];
		vars[i*varsPerVert + 1] = v1[0][1] - rotPt[1];
		vars[i*varsPerVert + 2] = v1[0][2] - rotPt[2];
		vars[i*varsPerVert + 3] = q.x;
		vars[i*varsPerVert + 4] = q.y;
		vars[i*varsPerVert + 5] = q.z;
		vars[i*varsPerVert + 6] = q.w;
	}
#else
	for (i=0; i < curMesh->numPts(); i++) {
		vars[i*varsPerVert + 0] = curMesh->getPt(i)[0];
		vars[i*varsPerVert + 1] = curMesh->getPt(i)[1];
		vars[i*varsPerVert + 2] = curMesh->getPt(i)[2];
	}
#endif
//#endif

	for (i=0; i < mesh->numPts(); i++) {
		origVerts[i] = mesh->getPt(i);
	}

	for (i=0; i < numEdges; i++) {
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

	// save deformation
	{
		Vecd v(mesh->numPts()*12);
		for (i=0; i < mesh->numPts(); i++) {
			Mat4d m = QuatNorm(vars[i*7+3],vars[i*7+4],vars[i*7+5],vars[i*7+6]).toMatrixD();
			m[0][3] = vars[i*7+0];
			m[1][3] = vars[i*7+1];
			m[2][3] = vars[i*7+2];
			int j;
			for (j=0;j < 12; j++) {
				v[i*12+j] = m.n[j];
			}
		}
		FILE *f;
		if (!openFile(&f, "temp.disp", "wb", "displacements")) {
			return;
		}
		char version = '0';
		fwrite(&version, sizeof(char), 1, f);
		int size = v.size();
		fwrite(&size, sizeof(int), 1, f);
		fwrite(v.n, sizeof(double), size, f);
		fclose(f);
	}

	applyDef(vars);
}

void RotMatchGF::unfold(int reps) {
}

void RotMatchGF::applyDef(Vecd &variables) {
	int i;

	if (&variables != &vars)
		vars = variables;


	for (i=0; i < variables.size(); i += varsPerVert) {
#ifdef IMPLICIT_ROT
		curMesh->getPt(i/varsPerVert) = Vec3d(variables[i+0], variables[i+1], variables[i+2]);
#else
		Vec3d cp = matchVerts[i/varsPerVert];

		cp = QuatNorm(vars[i+3],vars[i+4],vars[i+5],vars[i+6]).toMatrixD() * cp;
		cp += Vec3d(vars[i+0], vars[i+1], vars[i+2]);
		curMesh->getPt(i/varsPerVert) = cp;
#endif
	}
	curMesh->calcNormals();
}

void RotMatchGF::zeroDeformation() {
/*	if (vars.size() != curMesh->numPts() * varsPerVert) {
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
	}*/
}

void RotMatchGF::addGradientVec(int index, Vec3d v) {
//	Vec3d cp = origVerts[index];
	int ofs = index*varsPerVert;

	grad[ofs + 0] += v[0];
	grad[ofs + 1] += v[1];
	grad[ofs + 2] += v[2];

#ifndef IMPLICIT_ROT
	Mat4d mat, derivs[4];
	QuatNorm q(vars[ofs+3], vars[ofs+4], vars[ofs+5], vars[ofs+6]);
	q.getMatrices(mat, derivs);

	grad[ofs + 3]+= (derivs[0] * matchVerts[index]) * v;
	grad[ofs + 4]+= (derivs[1] * matchVerts[index]) * v;
	grad[ofs + 5]+= (derivs[2] * matchVerts[index]) * v;
	grad[ofs + 6]+= (derivs[3] * matchVerts[index]) * v;
#endif
}

//static const double FD_DELTA = 0.000001;
static const double FD_DELTA = 1e-8;

double qDistDeriv(QuatNorm q0, QuatNorm q0d, QuatNorm q1, QuatNorm q1d) {
	double ret = 0;
	ret += (q0.x - q1.x) * ((q0d.x-q0.x)/FD_DELTA - ((q1d.x-q1.x)/FD_DELTA));
	ret += (q0.y - q1.y) * ((q0d.y-q0.y)/FD_DELTA - ((q1d.y-q1.y)/FD_DELTA));
	ret += (q0.z - q1.z) * ((q0d.z-q0.z)/FD_DELTA - ((q1d.z-q1.z)/FD_DELTA));
	ret += (q0.w - q1.w) * ((q0d.w-q0.w)/FD_DELTA - ((q1d.w-q1.w)/FD_DELTA));
	return ret*2;
}

double RotMatchGF::evaluateFunction(Vecd& variables) {
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


#ifdef IMPLICIT_ROT
	// calculate point rotations
	bool initPR = false;
	if (ptRots == NULL) {
		ptRots = new QuatNorm[curMesh->numPts()];
		ptRotDerivs = new QuatNorm*[curMesh->numPts()];
		tmNeigh = findTMNeighbors(curMesh);
		initPR = true;
	}
	for (i = 0; i < curMesh->numPts(); i++) {
		int numInvolved = min(8, (int)tmNeigh[i].size() + 1);
		Vec3d v0[8], v1[8];

		v0[0] = matchVerts[i];
		v1[0] = curMesh->getPt(i);
		for (j=1; j < numInvolved; j++) {
			v0[j] = matchVerts[tmNeigh[i][j-1]];
			v1[j] = curMesh->getPt(tmNeigh[i][j-1]);
		}

		ptRots[i] = matToQuat(ptsRot(v0, v1, numInvolved));
		

		// now, the derivatives
		if (initPR) {
			ptRotDerivs[i] = new QuatNorm[numInvolved*3];
		}
		for (j=0; j < numInvolved; j++) {
			for (k=0; k < 3; k++) {
				double oldVal = v1[j][k];
				v1[j][k] += FD_DELTA;
				ptRotDerivs[i][j*3+k] = matToQuat(ptsRot(v0, v1, numInvolved));
				v1[j][k] = oldVal;
			}
		}
	}
#endif

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

			Vec3d x1 = curMesh->getPt(ei.v0);
			Vec3d x2 = curMesh->getPt(ei.v1);

			double rest = ei.eRestLen;
			double delta = ei.eCurLen - rest;
			double weight = smoothnessWeight * neighWeights[ei.v0] * neighWeights[ei.v1];
			if (weight < 0) {
				cout << "bad smoothness weight. v0 = " << ei.v0 << " v1 = " << ei.v1 
					<< ", w0 = " <<
					neighWeights[ei.v0] << " w1 = " << neighWeights[ei.v1] << endl;
			}
			smoothErr += weight * sqr(delta);

			if (ei.eCurLen != 0) {
				addGradientVec(ei.v0, weight * 2.0 * delta * (1.0 / ei.eCurLen) * Vec3d(
					x1[0]-x2[0],
					x1[1]-x2[1],
					x1[2]-x2[2]));
				addGradientVec(ei.v1, -weight * 2.0 * delta * (1.0 / ei.eCurLen) * Vec3d(
					x1[0]-x2[0],
					x1[1]-x2[1],
					x1[2]-x2[2]));
			}


#ifdef IMPLICIT_ROT
			// rotation term
			weight = bendWeight * neighWeights[ei.v0] * neighWeights[ei.v1];
			QuatNorm q0 = ptRots[ei.v0];
			QuatNorm q1 = ptRots[ei.v1];
			double dist = sqr(q0.x-q1.x) + sqr(q0.y-q1.y) + sqr(q0.z-q1.z) + sqr(q0.w-q1.w);
//			double dist = quatDist(q0, q1);
			smoothErr += weight * dist;

			// derivative of rotations term
			// (these are tricky because some points are influenced by v0 and v1
			bool seenVert[8];
			memset(seenVert, 0, sizeof(bool)*8);
			int numInvolved0 = min(8, (int)tmNeigh[ei.v0].size() + 1);
			int numInvolved1 = min(8, (int)tmNeigh[ei.v1].size() + 1);
			for (j = 0; j < numInvolved0; j++) {
				int pt0 = (j == 0) ? ei.v0 : tmNeigh[ei.v0][j-1];

				for (k=0; k < numInvolved1; k++) {
					int pt1 = (k == 0) ? ei.v1 : tmNeigh[ei.v1][k-1];

					if (pt0 == pt1) {
						// case 1: pt0 and pt1 affect both v0 and v1
						//grad[pt0*3+0] += weight * (quatDist(ptRotDerivs[ei.v0][j*3+0], ptRotDerivs[ei.v1][k*3+0]) - dist) / FD_DELTA;
						//grad[pt0*3+1] += weight * (quatDist(ptRotDerivs[ei.v0][j*3+1], ptRotDerivs[ei.v1][k*3+1]) - dist) / FD_DELTA;
						//grad[pt0*3+2] += weight * (quatDist(ptRotDerivs[ei.v0][j*3+2], ptRotDerivs[ei.v1][k*3+2]) - dist) / FD_DELTA;
						grad[pt0*3+0] += weight * qDistDeriv(q0, ptRotDerivs[ei.v0][j*3+0], q1, ptRotDerivs[ei.v1][k*3+0]);
						grad[pt0*3+1] += weight * qDistDeriv(q0, ptRotDerivs[ei.v0][j*3+1], q1, ptRotDerivs[ei.v1][k*3+1]);
						grad[pt0*3+2] += weight * qDistDeriv(q0, ptRotDerivs[ei.v0][j*3+2], q1, ptRotDerivs[ei.v1][k*3+2]);
						seenVert[k] = true;
						break;
					}
				}

				if (k == numInvolved1) {
					// case 2: pt0 affects only v0
					//grad[pt0*3+0] += weight * (quatDist(ptRotDerivs[ei.v0][j*3+0], ptRots[ei.v1]) - dist) / FD_DELTA;
					//grad[pt0*3+1] += weight * (quatDist(ptRotDerivs[ei.v0][j*3+1], ptRots[ei.v1]) - dist) / FD_DELTA;
					//grad[pt0*3+2] += weight * (quatDist(ptRotDerivs[ei.v0][j*3+2], ptRots[ei.v1]) - dist) / FD_DELTA;
					grad[pt0*3+0] += weight * qDistDeriv(q0, ptRotDerivs[ei.v0][j*3+0], q1, q1);
					grad[pt0*3+1] += weight * qDistDeriv(q0, ptRotDerivs[ei.v0][j*3+1], q1, q1);
					grad[pt0*3+2] += weight * qDistDeriv(q0, ptRotDerivs[ei.v0][j*3+2], q1, q1);
				}
			}
			for (k=0; k < numInvolved1; k++) {
				int pt1 = (k == 0) ? ei.v1 : tmNeigh[ei.v1][k-1];
				if (!seenVert[k]) {
					// case 3: pt1 affects only v1
					//grad[pt1*3+0] += weight * (quatDist(ptRots[ei.v0], ptRotDerivs[ei.v1][k*3+0]) - dist) / FD_DELTA;
					//grad[pt1*3+1] += weight * (quatDist(ptRots[ei.v0], ptRotDerivs[ei.v1][k*3+1]) - dist) / FD_DELTA;
					//grad[pt1*3+2] += weight * (quatDist(ptRots[ei.v0], ptRotDerivs[ei.v1][k*3+2]) - dist) / FD_DELTA;
					grad[pt1*3+0] += weight * qDistDeriv(q0, q0, q1, ptRotDerivs[ei.v1][k*3+0]);
					grad[pt1*3+1] += weight * qDistDeriv(q0, q0, q1, ptRotDerivs[ei.v1][k*3+1]);
					grad[pt1*3+2] += weight * qDistDeriv(q0, q0, q1, ptRotDerivs[ei.v1][k*3+2]);
				}
			}
#else
			// rotation term
			weight = bendWeight * neighWeights[ei.v0] * neighWeights[ei.v1];

			Vec4d q0(vars[ei.v0*varsPerVert + 3], vars[ei.v0*varsPerVert + 4], vars[ei.v0*varsPerVert + 5], vars[ei.v0*varsPerVert + 6]);
			Vec4d q1(vars[ei.v1*varsPerVert + 3], vars[ei.v1*varsPerVert + 4], vars[ei.v1*varsPerVert + 5], vars[ei.v1*varsPerVert + 6]);
			int sign = 1;
			if ((q0+q1).length2() < (q0-q1).length2()) {
				sign = -1;
			}

			const double EXPONENT = 2;
			for (j=0; j < varsPerVert; j++) {
//				if (j == 3)
//					weight = weight * 0.01;
				if (j >= 3 && sign < 0) {
					double delta = vars[ei.v0*varsPerVert + j] + vars[ei.v1*varsPerVert + j];
					
					if (weight < 0) {
						cout << "bad bend weight. v0 = " << ei.v0 << " v1 = " 
							<< ei.v1 << ", w0 = " << neighWeights[ei.v0] 
							<< " w1 = " << neighWeights[ei.v1] << endl;
					}
					smoothErr += weight * log(1.0 + sqr(delta));
					grad[ei.v0*varsPerVert + j] += weight * (1.0/(1.0+sqr(delta))) * 2.0 * delta;
					grad[ei.v1*varsPerVert + j] += weight * (1.0/(1.0+sqr(delta))) * 2.0 * delta;

					/*int dSign = 1;
					if (delta < 0) {
						dSign = -1;
						delta = -delta;
					}

					smoothErr += weight * pow(delta, EXPONENT);
					if (delta > 1e-7) {
						grad[ei.v0*varsPerVert + j] += EXPONENT * weight * dSign * pow(delta, EXPONENT-1);
						grad[ei.v1*varsPerVert + j] += EXPONENT * weight * dSign * pow(delta, EXPONENT-1);
					}*/
					//smoothErr += weight * sqr(delta);
					//grad[ei.v0*varsPerVert + j] += 2.0 * weight * delta;
					//grad[ei.v1*varsPerVert + j] += 2.0 * weight * delta;
				}
				else {
					double delta = vars[ei.v0*varsPerVert + j] - vars[ei.v1*varsPerVert + j];
					smoothErr += weight * log(1 + sqr(delta));
					grad[ei.v0*varsPerVert + j] += weight * (1.0/(1.0+sqr(delta))) * 2.0 * delta;
					grad[ei.v1*varsPerVert + j] -= weight * (1.0/(1.0+sqr(delta))) * 2.0 * delta;

					
/*					int dSign = 1;
					if (delta < 0) {
						dSign = -1;
						delta = -delta;
					}


					smoothErr += weight * pow(delta, EXPONENT);
					if (delta > 1e-7) {
						grad[ei.v0*varsPerVert + j] += EXPONENT * weight * dSign * pow(delta, EXPONENT-1);
						grad[ei.v1*varsPerVert + j] -= EXPONENT * weight * dSign * pow
(delta, EXPONENT-1);
					}*/
/*					smoothErr += weight * sqr(delta);
					grad[ei.v0*varsPerVert + j] += 2.0 * weight * delta;
					grad[ei.v1*varsPerVert + j] -= 2.0 * weight * delta;*/
				}
			}

/*			QuatNorm q0(vars[ei.v0*varsPerVert + 3], vars[ei.v0*varsPerVert + 4], vars[ei.v0*varsPerVert + 5], vars[ei.v0*varsPerVert + 6]);
			QuatNorm q1(vars[ei.v1*varsPerVert + 3], vars[ei.v1*varsPerVert + 4], vars[ei.v1*varsPerVert + 5], vars[ei.v1*varsPerVert + 6]);

			smoothErr += weight * (sqr(q0.x - q1.x) + sqr(q0.y - q1.y) + sqr(q0.z - q1.z) + sqr(q0.w - q1.w));
			grad[ei.v0*varsPerVert + 3] += weight * 2.0 * (q0.x - q1.x);
			grad[ei.v0*varsPerVert + 4] += weight * 2.0 * (q0.y - q1.y);
			grad[ei.v0*varsPerVert + 5] += weight * 2.0 * (q0.z - q1.z);
			grad[ei.v0*varsPerVert + 6] += weight * 2.0 * (q0.w - q1.w);
			grad[ei.v1*varsPerVert + 3] -= weight * 2.0 * (q0.x - q1.x);
			grad[ei.v1*varsPerVert + 4] -= weight * 2.0 * (q0.y - q1.y);
			grad[ei.v1*varsPerVert + 5] -= weight * 2.0 * (q0.z - q1.z);
			grad[ei.v1*varsPerVert + 6] -= weight * 2.0 * (q0.w - q1.w);*/
/*
			double dist = quatDist(q0, q1);
			smoothErr += weight * dist;
			
			q0.x += FD_DELTA;
			grad[ei.v0*varsPerVert + 3] += weight * (quatDist(q0, q1) - dist) / FD_DELTA;
			q0.x = vars[ei.v0*varsPerVert + 3];
			q0.y += FD_DELTA;
			grad[ei.v0*varsPerVert + 4] += weight * (quatDist(q0, q1) - dist) / FD_DELTA;
			q0.y = vars[ei.v0*varsPerVert + 4];
			q0.z += FD_DELTA;
			grad[ei.v0*varsPerVert + 5] += weight * (quatDist(q0, q1) - dist) / FD_DELTA;
			q0.z = vars[ei.v0*varsPerVert + 5];
			q0.w += FD_DELTA;
			grad[ei.v0*varsPerVert + 6] += weight * (quatDist(q0, q1) - dist) / FD_DELTA;
			q0.w = vars[ei.v0*varsPerVert + 6];
			q1.x += FD_DELTA;
			grad[ei.v1*varsPerVert + 3] += weight * (quatDist(q0, q1) - dist) / FD_DELTA;
			q1.x = vars[ei.v1*varsPerVert + 3];
			q1.y += FD_DELTA;
			grad[ei.v1*varsPerVert + 4] += weight * (quatDist(q0, q1) - dist) / FD_DELTA;
			q1.y = vars[ei.v1*varsPerVert + 4];
			q1.z += FD_DELTA;
			grad[ei.v1*varsPerVert + 5] += weight * (quatDist(q0, q1) - dist) / FD_DELTA;
			q1.z = vars[ei.v1*varsPerVert + 5];
			q1.w += FD_DELTA;
			grad[ei.v1*varsPerVert + 6] += weight * (quatDist(q0, q1) - dist) / FD_DELTA;
			q1.w = vars[ei.v1*varsPerVert + 6];*/

#endif

		}
	}

	/*
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
	*/

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
