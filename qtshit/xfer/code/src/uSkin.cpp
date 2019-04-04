#include "uSkin.h"
#include "uSolver.h"
#include "skeleton.h"
#include "trimesh.h"
#include "normalMap.h"
#include "vl/VLd.h"
#include <float.h>

const bool USE_LBS = true;
//const bool USE_LBS = false;

// from uMaster
extern Vec3d *normalMap;


// UPoseDepDef ======================================================

Vec3d *UPoseDepDef::load(USkin *iSkin, TriMesh *mesh, ifstream &f) {
	char s[80];
	int i;

	s[0] = 0;
	f >> s;
	if (strlen(s) < 1) {
		cout << "end of file" << endl;
		return NULL;
	}
	skin = iSkin;
	transInd = skin->skel->transforms.lookupName(s);
	if (transInd < 0) {
		cout << "WARNING: unknown transform: " << s << endl;
		return NULL;
	}
	transform = skin->skel->transforms.getT(transInd);

	f >> numSamples;
	double *samples = new double[4 * numSamples];
	for (i = 0; i < numSamples*4; i++) {
		f >> samples[i];
	}

	rbf = new RBF(numSamples, 4, true);
	rbf->init(samples);

	// load colors
	f >> s;
	FILE *f2;
	if (!openFile(&f2, s, "rb", "color map")) {
		return false;
	}
	char version;
	fread(&version, sizeof(char), 1, f2);
	if (version != '0') {
		cout << "unknown version; treating as raw data" << endl;
		fseek(f2, 0, SEEK_SET);
	}
	int size;
	fread(&size, sizeof(int), 1, f2);
	if (size != mesh->numPts() * 3) {
		cout << "size mismatch in " << s << endl;
		fclose(f2);
		return false;
	}
	Vec3d *colorData = new Vec3d[mesh->numPts()];
	fread(colorData, sizeof(double), size, f2);
	fclose(f2);

	return colorData;
}

void UPoseDepDef::mirrorFrom(UPoseDepDef *origPDD, int iTransInd) {
	int i;

	transInd = iTransInd;
	skin = origPDD->skin;
	transform = skin->skel->transforms.getT(transInd);
	numSamples = origPDD->numSamples;

	double *samples = new double[4 * numSamples];
	for (i=0; i < numSamples; i++) {
		samples[i*4 + 0] = origPDD->rbf->samplePts[i][0];
		samples[i*4 + 1] = -origPDD->rbf->samplePts[i][1];
		samples[i*4 + 2] = origPDD->rbf->samplePts[i][2];
		samples[i*4 + 3] = -origPDD->rbf->samplePts[i][3];
	}

	rbf = new RBF(numSamples, 4, true);
	rbf->init(samples);

	isMirror = true;
}


// USkin ============================================================

USkin::USkin() {
	int i;

	skel = NULL;
	numPts = 0;
	dressPts = NULL;
	pddDressPts = NULL;
	dressJoints = NULL;
	numTransInit = 0;
	tiFrames = NULL;
	tiMarkers = NULL;
	maxInf = 0;
	infJoints = NULL;
	infWeights = NULL;
	curMats = NULL;
	curJointQuat = NULL;
	curJointPos = NULL;
	curPts = NULL;
	curJointDerivs = NULL;

	pddPtKeys = NULL;
	pddNMKeys = NULL;

	baseNormalMap = new Vec3d[NORMAL_MAP_SIZE];
	for (i=0; i < NORMAL_MAP_SIZE; i++)
		baseNormalMap[i] = Vec3d(0, 0, -1);
}

void USkin::init(Skeleton *iSkel, int iNumPts, int iMaxInf) {
	int i;

	skel = iSkel;
	numPts = iNumPts;
	maxInf = iMaxInf;

	// allocate variables
	dressPts = new Vec3d[numPts];
	pddDressPts = new Vec3d[numPts];
	dressJoints = new Vec3d[skel->transforms.size()];
	ptNames = new char80[numPts];
	pddNMIndex = new vector<int>[NORMAL_MAP_SIZE];
	pddPtIndex = new vector<int>[numPts];

	infJoints = new int[numPts * maxInf];
	for (i=0; i < numPts * maxInf; i++) infJoints[i] = -1;
	infWeights = new double[numPts * maxInf];

	curMats = new Mat4d[skel->transforms.size()];
	curJointQuat = new QuatNorm[skel->transforms.size()];
	curJointPos = new Vec3d[skel->transforms.size()];
	curPts = new Vec3d[numPts];
	curJointDerivs = new Vec4d[skel->transforms.size() * skel->numDofs];

	dofIndex = new int[skel->transforms.size()];
	int pDof = 0;
	int iDof = 0;
	for (i=0; i < skel->transforms.size(); i++) {
		SkelTransform *curTrans = skel->transforms.getT(i);
		if (curTrans->isIntrinsic) {
			dofIndex[i] = iDof;
			iDof += curTrans->numDofs();
		}
		else {
			dofIndex[i] = pDof;
			pDof += curTrans->numDofs();
		}
	}
}

void USkin::transInit(UExample *ex) {
	// initialize the skeleton translations based on tiMarkers
	int i;

//	skel->zero();
	skel->updateCoords();
	for (i=0; i < numTransInit; i++) {
		SkelTranslation *curTrans = (SkelTranslation*)skel->transforms.getT(tiFrames[i]);
		if (!curTrans) {
			cout << "WARNING: don't recognize transform " << tiFrames[i] << endl;
			continue;
		}

		Vec3d v0, v1;
		ex->getPt(tiMarkers[i*2], &v0);
		if (tiMarkers[i*2+1] >= 0) {
			ex->getPt(tiMarkers[i*2+1], &v1);
			curTrans->curVal = 0.5 * (v0 + v1) - curTrans->globalCoord.v;
		}
		else {
			curTrans->curVal = v0 - curTrans->globalCoord.v;
		}
		skel->updateCoords();
	}
}

int USkin::getPddJoint(int pt, int inf) {
	int joint = infJoints[pt * maxInf + inf];

	// hack: use masks from single-axis joint on combined joints
	if (joint == 6)	// lClavicle
		joint = 8;
	else if (joint == 10) // lShoulder
		joint = 13;
	else if (joint == 24) // rClavicle
		joint = 26;
	else if (joint == 28) // rShoulder
		joint = 31;

	return joint;
}

void USkin::getSkinMat(int pt, Mat3d &m, Vec3d &base) {
	int inf, j, joint;

	m = Mat3d(0,0,0, 0,0,0, 0,0,0);
	base = Vec3d();

	if (USE_LBS) {
		// LBS
		for (inf=0; inf < maxInf; inf++) {
			int joint = infJoints[pt*maxInf + inf];
			double weight = infWeights[pt*maxInf + inf];
			if (joint >= 0 && weight != 0) {
				Mat4d &cm = curMats[joint];
				Mat3d cm3(cm[0][0],cm[0][1],cm[0][2],cm[1][0],cm[1][1],cm[1][2],cm[2][0],cm[2][1],cm[2][2]);
				m += weight * cm3;
				base += weight * (Vec3d(cm[0][3],cm[1][3],cm[2][3]) - cm3 * dressJoints[joint]);
			}
		}
	}
	else {
		// rotation-interpolated skinning 
		QuatNorm q(0,0,0,0);
		base = Vec3d();
		Vec3d ojc, njc;
		for (inf=0; inf < maxInf; inf++) {
			int joint = infJoints[pt*maxInf + inf];
			double weight = infWeights[pt*maxInf + inf];
			if (joint >= 0 && weight != 0) {
				Mat4d &cm = curMats[joint];
				Mat3d cm3(cm[0][0],cm[0][1],cm[0][2],cm[1][0],cm[1][1],cm[1][2],cm[2][0],cm[2][1],cm[2][2]);
				base += weight * (Vec3d(cm[0][3],cm[1][3],cm[2][3]) - cm3 * dressJoints[joint]);
				ojc += weight * dressJoints[joint];
				njc += weight * curJointPos[joint];
				if (q.x != 0 || q.y != 0 || q.z != 0 || q.w != 0) {
					if (q.x * curJointQuat[joint].x + 
						q.y * curJointQuat[joint].y +
						q.z * curJointQuat[joint].z +
						q.w * curJointQuat[joint].w < 0) {
						curJointQuat[joint].x *= -1;
						curJointQuat[joint].y *= -1;
						curJointQuat[joint].z *= -1;
						curJointQuat[joint].w *= -1;
					}
				}
				q += weight * curJointQuat[joint];
			}
		}
		q.normalize();
		Mat4d cm = q.toMatrixD();
		m = Mat3d(cm[0][0],cm[0][1],cm[0][2],cm[1][0],cm[1][1],cm[1][2],cm[2][0],cm[2][1],cm[2][2]);
		base = njc - m * ojc;
	}

/*
	// spherical blend skinning
	int numInf = 0;
	for (inf=0; inf < maxInf; inf++) {
		if (infJoints[pt*maxInf + inf] >= 0)
			numInf++;
	}

	if (numInf == 1) {
		// if only one influence, just do LBS
		for (inf=0; inf < maxInf; inf++) {
			int joint = infJoints[pt*maxInf + inf];
			double weight = infWeights[pt*maxInf + inf];
			if (joint >= 0 && weight != 0) {
				Mat4d &cm = curMats[joint];
				Mat3d cm3(cm[0][0],cm[0][1],cm[0][2],cm[1][0],cm[1][1],cm[1][2],cm[2][0],cm[2][1],cm[2][2]);
				m += weight * cm3;
				base += weight * (Vec3d(cm[0][3],cm[1][3],cm[2][3]) - cm3 * dressJoints[joint]);
			}
		}
		return;
	}

	int numCombo = (numInf * (numInf-1)) / 2;
	TMat dMat(numCombo * 3, 3);
	TVec eVec(numCombo * 3);
	int ind = 0, j0, j1;
	for (j0 = 0; j0 < maxInf; j0++) {
		if (infJoints[pt*maxInf + j0] < 0)
			continue;
		Mat4d &m0 = curMats[infJoints[pt*maxInf + j0]];
		for (j1 = j0+1; j1 < maxInf; j1++) {
			if (infJoints[pt*maxInf + j1] < 0)
				continue;
			Mat4d &m1 = curMats[infJoints[pt*maxInf + j1]];

			dMat[ind*3+0][0] = m0[0][0] - m1[0][0];
			dMat[ind*3+0][1] = m0[0][1] - m1[0][1];
			dMat[ind*3+0][2] = m0[0][2] - m1[0][2];
			dMat[ind*3+1][0] = m0[1][0] - m1[1][0];
			dMat[ind*3+1][1] = m0[1][1] - m1[1][1];
			dMat[ind*3+1][2] = m0[1][2] - m1[1][2];
			dMat[ind*3+2][0] = m0[2][0] - m1[2][0];
			dMat[ind*3+2][1] = m0[2][1] - m1[2][1];
			dMat[ind*3+2][2] = m0[2][2] - m1[2][2];
			eVec[ind*3+0] = m1[0][3] - m0[0][3];
			eVec[ind*3+1] = m1[1][3] - m0[1][3];
			eVec[ind*3+2] = m1[2][3] - m0[2][3];

			ind++;
		}
	}

	// calculate pseudoinverse
	TMat U(numCombo * 3, 3), V(3, 3);
	TVec diagonal(3);
	SVDFactorization(dMat, U, V, diagonal);
	TMat pinv(3, numCombo * 3), diag(3, 3);
	diag = vl_0;
	for (j = 0; j < 3; j++) {
		if (fabs(diagonal[j]) > 0.00001)
			diag.Elt(j, j) = 1.0 / diagonal[j];
	}
	pinv = V * diag * trans(U);
	diagonal = pinv * eVec;
	Vec3d rc(diagonal[0], diagonal[1], diagonal[2]);

	QuatNorm q(0,0,0,0), qPivot(0,0,0,0);

	for (j=0; j < maxInf; j++) {
		joint = infJoints[pt*maxInf + j];
		if (joint >= 0 && infWeights[pt*maxInf + j] != 0) {
			QuatNorm q2 = skel->transforms.getT(joint)->globalCoord.q;
			if (qPivot.w == 0) {
				qPivot = q2;
			}
			else {
				if (q2.x*qPivot.x + q2.y*qPivot.y + q2.z*qPivot.z + q2.w*qPivot.w < 0)
					q2 = -1 * q2;
			}
			q += infWeights[pt*maxInf + j] * q2;
			base += infWeights[pt*maxInf + j] * (curMats[joint] * (rc - dressJoints[joint]));
		}
	}
	q.normalize();
	Mat4d m4 = q.toMatrixD();
	m[0][0] += m4[0][0];
	m[0][1] += m4[0][1];
	m[0][2] += m4[0][2];
	m[1][0] += m4[1][0];
	m[1][1] += m4[1][1];
	m[1][2] += m4[1][2];
	m[2][0] += m4[2][0];
	m[2][1] += m4[2][1];
	m[2][2] += m4[2][2];
	base -= m4 * rc;*/
}

void USkin::updateMats() {
	int i;
	for (i=0; i < skel->transforms.size(); i++) {
		curMats[i] = skel->transforms.getT(i)->globalCoord.mat;
		curJointQuat[i] = skel->transforms.getT(i)->globalCoord.q;
		curJointQuat[i].normalize();
		curJointPos[i] = skel->transforms.getT(i)->globalCoord.v;
	}
}

void USkin::updateJoints() {
	int i, j;
	int index = 0;

	static Skeleton *tempSkel = NULL;
	if (!tempSkel) {
		tempSkel = new Skeleton();
		tempSkel->copyFrom(skel);

		tempSkel->allocDerivs(tempSkel->numDofs);
	}

	// copy only the intrinsic dofs into a temporary skeleton
	for (i=0; i < tempSkel->transforms.size(); i++) {
		SkelTransform *curTrans = tempSkel->transforms.getT(i);
		if (curTrans->isIntrinsic && strstr(curTrans->name, "Carry") == NULL)
			curTrans->copyVal(skel->transforms.getT(i));
		else
			curTrans->zero();
	}
	tempSkel->updateCoords();
	tempSkel->updateDerivs();
	tempSkel->updateGlobalDerivs();

	// calculate the joint positions and derivatives
	for (i=0; i < tempSkel->transforms.size(); i++) {
		SkelTransform *curTrans = tempSkel->transforms.getT(i);
		dressJoints[i] = curTrans->globalCoord.v;
		for (j=0; j < tempSkel->numDofs; j++)
			curJointDerivs[index++] = curTrans->globalDerivs[j] * Vec4d(0,0,0,1);
	}
}

void USkin::updateRBFs() {
	int pdd;

	for (pdd=0; pdd < pdds.size(); pdd++) {
		QuatNorm q = pdds[pdd]->transform->curCoord.q;
		Vec4d v(q.x, q.y, q.z, q.w);
		pdds[pdd]->rbf->evalDerivs(v.n);
	}
}

void USkin::updatePts(int minPt, int maxPt, bool updateNM) {
	int i, j, w;
	int joint;
	if (maxPt < 0)
		maxPt = numPts;

	updateRBFs();

	memcpy(pddDressPts, dressPts, sizeof(Vec3d)*numPts);
	for (i=0; i < numPts; i++) {
		for (j=0; j < pddPtIndex[i].size(); j += 2) {
			RBF *rbf = pdds[pddPtIndex[i][j]]->rbf;
			int ofs = pddPtIndex[i][j+1];
			double yMult = 1.0;
			if (pdds[pddPtIndex[i][j]]->isMirror)
				yMult = -1.0;
			for (w=1; w < rbf->mN; w++) {
				pddDressPts[i] += rbf->curWeights[w] * prod(pddPtKeys[ofs + w-1], Vec3d(1, yMult, 1));
			}
		}
	}

	if (normalMap && updateNM) {
		memcpy(normalMap, baseNormalMap, sizeof(Vec3d)*NORMAL_MAP_SIZE);

		for (i=0; i < NORMAL_MAP_SIZE; i++) {
			for (j=0; j < pddNMIndex[i].size(); j += 2) {
				RBF *rbf = pdds[pddNMIndex[i][j]]->rbf;
				int ofs = pddNMIndex[i][j+1];
				double yMult = 1.0;
				if (pdds[pddNMIndex[i][j]]->isMirror)
					yMult = -1.0;
				for (w=1; w < rbf->mN; w++) {
					normalMap[i] += rbf->curWeights[w] * prod(pddNMKeys[ofs + w-1], Vec3d(1, yMult, 1));
				}
			}
		}
	}

/*
	// spherical blend skinning
	for (i=minPt; i < maxPt; i++) {
		curPts[i] = Vec3d();

		int numInf = 0;
		for (j=0; j < maxInf; j++) {
			if (infJoints[i*maxInf + j] >= 0)
				numInf++;
		}

		if (numInf == 1) {
			// if only one influence, just do regular skinning
			curPts[i] = Vec3d();
			for (j=0; j < maxInf; j++) {
				joint = infJoints[i*maxInf + j];
				if (joint >= 0 && infWeights[i*maxInf + j] != 0) {
					curPts[i] += infWeights[i*maxInf + j] * (curMats[joint] * (pddDressPts[i] - dressJoints[joint]));
				}
			}
			continue;
		}

		int numCombo = (numInf * (numInf-1)) / 2;
		TMat dMat(numCombo * 3, 3);
		TVec eVec(numCombo * 3);
		int ind = 0, j0, j1;
		for (j0 = 0; j0 < maxInf; j0++) {
			if (infJoints[i*maxInf + j0] < 0)
				continue;
			Mat4d &m0 = curMats[infJoints[i*maxInf + j0]];
			for (j1 = j0+1; j1 < maxInf; j1++) {
				if (infJoints[i*maxInf + j1] < 0)
					continue;
				Mat4d &m1 = curMats[infJoints[i*maxInf + j1]];

				dMat[ind*3+0][0] = m0[0][0] - m1[0][0];
				dMat[ind*3+0][1] = m0[0][1] - m1[0][1];
				dMat[ind*3+0][2] = m0[0][2] - m1[0][2];
				dMat[ind*3+1][0] = m0[1][0] - m1[1][0];
				dMat[ind*3+1][1] = m0[1][1] - m1[1][1];
				dMat[ind*3+1][2] = m0[1][2] - m1[1][2];
				dMat[ind*3+2][0] = m0[2][0] - m1[2][0];
				dMat[ind*3+2][1] = m0[2][1] - m1[2][1];
				dMat[ind*3+2][2] = m0[2][2] - m1[2][2];
				eVec[ind*3+0] = m1[0][3] - m0[0][3];
				eVec[ind*3+1] = m1[1][3] - m0[1][3];
				eVec[ind*3+2] = m1[2][3] - m0[2][3];

				ind++;
			}
		}

		// calculate pseudoinverse
		TMat U(numCombo * 3, 3), V(3, 3);
		TVec diagonal(3);
		SVDFactorization(dMat, U, V, diagonal);
		TMat pinv(3, numCombo * 3), diag(3, 3);
		diag = vl_0;
		for (j = 0; j < 3; j++) {
			if (fabs(diagonal[j]) > 0.00001)
				diag.Elt(j, j) = 1.0 / diagonal[j];
		}
		pinv = V * diag * trans(U);
		diagonal = pinv * eVec;
		Vec3d rc(diagonal[0], diagonal[1], diagonal[2]);

		QuatNorm q(0,0,0,0), qPivot(0,0,0,0);

		for (j=0; j < maxInf; j++) {
			joint = infJoints[i*maxInf + j];
			if (joint >= 0 && infWeights[i*maxInf + j] != 0) {
				QuatNorm q2 = skel->transforms.getT(joint)->globalCoord.q;
				if (qPivot.w == 0) {
					qPivot = q2;
				}
				else {
					if (q2.x*qPivot.x + q2.y*qPivot.y + q2.z*qPivot.z + q2.w*qPivot.w < 0)
						q2 = -1 * q2;
				}
				q += infWeights[i*maxInf + j] * q2;
				curPts[i] += infWeights[i*maxInf + j] * (curMats[joint] * (rc - dressJoints[joint]));
			}
		}
		q.normalize();
		curPts[i] += q.toMatrixD() * (pddDressPts[i] - rc);
	}
*/

	if (USE_LBS) {
		// linear blend skinning
		for (i=minPt; i < maxPt; i++) {
			curPts[i] = Vec3d();
			for (j=0; j < maxInf; j++) {
				joint = infJoints[i*maxInf + j];
				if (joint >= 0 && infWeights[i*maxInf + j] != 0) {
					curPts[i] += infWeights[i*maxInf + j] * (curMats[joint] * (pddDressPts[i] - dressJoints[joint]));
				}
			}
		}
	}
	else {
		// rotation-interpolated skinning
		for (i=minPt; i < maxPt; i++) {
/*			Vec3d origPos, targetPos;
			QuatNorm targetRot(0,0,0,0);
			
			for (j=0; j < maxInf; j++) {
				joint = infJoints[i*maxInf + j];
				if (joint >= 0 && infWeights[i*maxInf + j] != 0) {
					origPos += infWeights[i*maxInf + j] * dressJoints[joint];
					if (targetRot.x != 0 || targetRot.y != 0 || targetRot.z != 0 || targetRot.w != 0) {
						if (targetRot.x * curJointQuat[joint].x + 
							targetRot.y * curJointQuat[joint].y +
							targetRot.z * curJointQuat[joint].z +
							targetRot.w * curJointQuat[joint].w < 0) {
							curJointQuat[joint].x *= -1;
							curJointQuat[joint].y *= -1;
							curJointQuat[joint].z *= -1;
							curJointQuat[joint].w *= -1;
						}
					}
					targetPos += infWeights[i*maxInf + j] * curJointPos[joint];
					targetRot += infWeights[i*maxInf + j] * curJointQuat[joint];
				}
			}

			curPts[i] = targetPos + (targetRot.toMatrixD() * (pddDressPts[i] - origPos));
*/
			Mat3d m;
			Vec3d base;
			getSkinMat(i, m, base);
			curPts[i] = base + m * pddDressPts[i];
		}
	}
}

void USkin::calcVecGrad(int pt, 
					 Vec3d *dressGrad, Vec3d *intGrad,
					 Vec3d *poseGrad, Vec3d *weightGrad, Vec3d *pddGrad) {
	int inf, pddInd, samp;
	Vec3d dv;
	int joint;

	for (inf=0; inf < maxInf; inf++) {
		if (infJoints[pt*maxInf + inf] >= 0) {
			double w = infWeights[pt*maxInf + inf];
			joint = infJoints[pt*maxInf + inf];
			Vec3d local = pddDressPts[pt] - dressJoints[joint];
			SkelTransform *curTrans = skel->transforms.getT(joint);

			if (pddGrad /*|| poseGrad*/) {
				// pose-dependent deformation variables
				for (pddInd = 0; pddInd < pddPtIndex[pt].size(); pddInd += 2) {
					UPoseDepDef *pdd = pdds[pddPtIndex[pt][pddInd]];
					int ofs = pddPtIndex[pt][pddInd+1];
					double yMult = (pdd->isMirror ? -1 : 1);

					for (samp = 0; samp < pdd->rbf->mN - 1; samp++) {
						double w2 = pdd->rbf->curWeights[samp + 1];

						if (pddGrad && w2 != 0) {
							pddGrad[ofs*3 + samp*3 + 0] += w * w2 *
								vec4to3(curMats[joint] * Vec4d(1, 0, 0, 0));
							pddGrad[ofs*3 + samp*3 + 1] += yMult * w * w2 *
								vec4to3(curMats[joint] * Vec4d(0, 1, 0, 0));
							pddGrad[ofs*3 + samp*3 + 2] += w * w2 *
								vec4to3(curMats[joint] * Vec4d(0, 0, 1, 0));
						}

						/*if (poseGrad) {
							Vec4d key = Vec4d(pddPtKeys[ofs + samp][0], 
								pddPtKeys[ofs + samp][1],
								pddPtKeys[ofs + samp][2], 0);
							if (strcmp(pdd->transform->className, "SkelQuatRotation") == 0) {
								int dof;
								for (dof = 0; dof < 4; dof++) {
									poseGrad[pdd->transform->dofInd + dof] += 
										pdd->rbf->curDerivs[dof*pdd->rbf->mN + samp + 1] * w * 
										vec4to3(curMats[joint] * key);
								}
							}
							else if (strcmp(pdd->transform->className, "SkelEulerRotation") == 0) {
								// chain rule from euler angle to quaternion
								int axis = ((SkelEulerRotation*)pdd->transform)->axis;
								double angle = ((SkelEulerRotation*)pdd->transform)->curAngle;
								poseGrad[pdd->transform->dofInd] += 
									0.5 * cos(0.5 * angle) *
									pdd->rbf->curDerivs[axis*pdd->rbf->mN + samp + 1] * w * 
									vec4to3(curMats[joint] * key);
								poseGrad[pdd->transform->dofInd] += 
									0.5 * sin(0.5 * angle) *
									pdd->rbf->curDerivs[3*pdd->rbf->mN + samp + 1] * w * 
									vec4to3(curMats[joint] * key);
							}
							else if (strcmp(pdd->transform->name+1, "ShoulderCR") == 0) {
								SkelCombinedTransform *shoulder = 
									(SkelCombinedTransform*)pdd->transform;
								SkelEulerRotation *euler[3];
								euler[0] = (SkelEulerRotation*)shoulder->orig[0];
								euler[1] = (SkelEulerRotation*)shoulder->orig[1];
								euler[2] = (SkelEulerRotation*)shoulder->orig[2];
								double cx = cos(euler[0]->curAngle/2);
								double sx = sin(euler[0]->curAngle/2);
								double cy = cos(euler[2]->curAngle/2);
								double sy = sin(euler[2]->curAngle/2);
								double cz = cos(euler[1]->curAngle/2);
								double sz = sin(euler[1]->curAngle/2);
								Vec4d dqdx;
								// ShoulderXR
								dqdx = Vec4d(
									 0.5*cy*cz*cx + 0.5*sy*sz*sx,
									-0.5*sy*cz*sx - 0.5*cy*sz*cx,
									-0.5*cy*sz*sx + 0.5*sy*cz*cx,
									 0.5*cy*cz*sx - 0.5*sy*sz*cx);
								poseGrad[euler[0]->dofInd] += 
									(dqdx[0] * pdd->rbf->curDerivs[0*pdd->rbf->mN + samp + 1] +
									dqdx[1] * pdd->rbf->curDerivs[1*pdd->rbf->mN + samp + 1] +
									dqdx[2] * pdd->rbf->curDerivs[2*pdd->rbf->mN + samp + 1] +
									dqdx[3] * pdd->rbf->curDerivs[3*pdd->rbf->mN + samp + 1]) 
									* w * vec4to3(curMats[joint] * key);
								// ShoulderZR
								dqdx = Vec4d(
									 -0.5*cy*sz*sx - 0.5*sy*cz*cx,
									 -0.5*sy*sz*cx - 0.5*cy*cz*sx,
									  0.5*cy*cz*cx - 0.5*sy*sz*sx,
									  0.5*cy*sz*cx - 0.5*sy*cz*sx);
								poseGrad[euler[1]->dofInd] += 
									(dqdx[0] * pdd->rbf->curDerivs[0*pdd->rbf->mN + samp + 1] +
									dqdx[1] * pdd->rbf->curDerivs[1*pdd->rbf->mN + samp + 1] +
									dqdx[2] * pdd->rbf->curDerivs[2*pdd->rbf->mN + samp + 1] +
									dqdx[3] * pdd->rbf->curDerivs[3*pdd->rbf->mN + samp + 1])
									* w * vec4to3(curMats[joint] * key);
								// ShoulderYR
								dqdx = Vec4d(
									 -0.5*sy*cz*sx - 0.5*cy*sz*cx,
									  0.5*cy*cz*cx + 0.5*sy*sz*sx,
									 -0.5*sy*sz*cx + 0.5*cy*cz*sx,
									 +0.5*sy*cz*cx - 0.5*cy*sz*sx);
								poseGrad[euler[2]->dofInd] += 
									(dqdx[0] * pdd->rbf->curDerivs[0*pdd->rbf->mN + samp + 1] +
									dqdx[1] * pdd->rbf->curDerivs[1*pdd->rbf->mN + samp + 1] +
									dqdx[2] * pdd->rbf->curDerivs[2*pdd->rbf->mN + samp + 1] +
									dqdx[3] * pdd->rbf->curDerivs[3*pdd->rbf->mN + samp + 1])
									* w * vec4to3(curMats[joint] * key);
							}
							else if (strcmp(pdd->transform->name+1, "ClavicleCQ") == 0) {
								SkelCombinedTransform *clav = 
									(SkelCombinedTransform*)pdd->transform;
								SkelEulerRotation *euler[2];
								euler[0] = (SkelEulerRotation*)clav->orig[0];
								euler[1] = (SkelEulerRotation*)clav->orig[1];
								double cx = cos(euler[1]->curAngle/2);
								double sx = sin(euler[1]->curAngle/2);
								double cz = cos(euler[0]->curAngle/2);
								double sz = sin(euler[0]->curAngle/2);
								Vec4d dqdx;
								// ClavicleZR
								dqdx = Vec4d(0.5*sx*sz, -0.5*sx*cz, -0.5*cx*cz, -0.5*cx*sz);
								poseGrad[euler[0]->dofInd] += 
									(dqdx[0] * pdd->rbf->curDerivs[0*pdd->rbf->mN + samp + 1] +
									dqdx[1] * pdd->rbf->curDerivs[1*pdd->rbf->mN + samp + 1] +
									dqdx[2] * pdd->rbf->curDerivs[2*pdd->rbf->mN + samp + 1] +
									dqdx[3] * pdd->rbf->curDerivs[3*pdd->rbf->mN + samp + 1])
									* w * vec4to3(curMats[joint] * key);
								// ClavicleXR
								dqdx = Vec4d(-0.5*cx*cz, -0.5*cx*sz, 0.5*sx*sz, -0.5*sx*cz);
								poseGrad[euler[1]->dofInd] += 
									(dqdx[0] * pdd->rbf->curDerivs[0*pdd->rbf->mN + samp + 1] +
									dqdx[1] * pdd->rbf->curDerivs[1*pdd->rbf->mN + samp + 1] +
									dqdx[2] * pdd->rbf->curDerivs[2*pdd->rbf->mN + samp + 1] +
									dqdx[3] * pdd->rbf->curDerivs[3*pdd->rbf->mN + samp + 1])
									* w * vec4to3(curMats[joint] * key);
							}
							else {
								cout << "WARNING: unsupported class type " << pdd->transform->className << " for joint " << pdd->transform->name << endl;
							}
						}*/
					}

				}
			}

			if (w != 0) {
				if (dressGrad) {
					// dressGrad
					dv = vec4to3(curMats[joint] * Vec4d(1, 0, 0, 0));
					dressGrad[0] += w * dv;
					dv = vec4to3(curMats[joint] * Vec4d(0, 1, 0, 0));
					dressGrad[1] += w * dv;
					dv = vec4to3(curMats[joint] * Vec4d(0, 0, 1, 0));
					dressGrad[2] += w * dv;
				}

				// pose-dependent gradients
				int dof = 0;
				int ip = 0;
				int pp = 0;
				for (dof = 0; dof < skel->numDofs; dof++) {
					if (skel->dofIntrins[dof]) {
						if (intGrad) 
							intGrad[ip++] += w * (curTrans->globalDerivs[dof] * local - vec4to3(curMats[joint] * curJointDerivs[joint*skel->numDofs + dof]));
					}
					else {
						if (poseGrad)
							poseGrad[pp++] += w * (curTrans->globalDerivs[dof] * local);
					}
				}
			}

			if (weightGrad)
				weightGrad[inf] += (curMats[joint] * local);
		}
	}
}

					 /*
double USkin::calcGradExp(int pt, Vec3d v,
					 double *dressGrad, double *intGrad,
					 double *poseGrad, double *weightGrad, double *pddGrad,
					 VLMatd &phi, USolver *uSolver) {
	int i;
	int joint, inf;
	int comp, comp2, tr;
	int numComponents = phi.Rows();
	double ret = 0;
	int ip, pp, dof;

	static Vec3d *aMat = NULL, *daDress, *daInt, *daPose, *daWeight, *daPDD;
	if (aMat == NULL) {
		aMat = new Vec3d[numComponents];
		daDress = new Vec3d[3];
		daInt = new Vec3d[skel->numIntrinsicDofs * (numComponents+1)];
		daPose = new Vec3d[skel->numPoseDofs * (numComponents+1)];
		daWeight = new Vec3d[(maxInf-1) * (numComponents+1)];
	}

	memset(aMat, 0, sizeof(Vec3d) * numComponents);
	memset(daDress, 0, sizeof(Vec3d) * 3);
	memset(daInt, 0, sizeof(Vec3d) * skel->numIntrinsicDofs * (numComponents+1));
	memset(daPose, 0, sizeof(Vec3d) * skel->numPoseDofs * (numComponents+1));
	memset(daWeight, 0, sizeof(Vec3d) * (maxInf-1) * (numComponents+1));

	Vec3d trans;
	for (inf=0; inf < maxInf; inf++) {
		ip = 0;
		pp = 0;

		if (infJoints[pt*maxInf + inf] >= 0) {
			double w = infWeights[pt*maxInf + inf];
			joint = infJoints[pt*maxInf + inf];
			SkelTransform *curTrans = skel->transforms.getT(joint);

			for (comp = 0; comp < numComponents; comp++) {
				Vec4d dress = vec3to4z(uSolver->varDressPt(pt, comp));
				aMat[comp] += w * vec4to3(curMats[joint] * dress);

				// intrinsic and pose derivatives
				for (dof = 0; dof < skel->numDofs; dof++) {
					if (skel->dofIntrins[dof]) {
						daInt[ip++] += w * vec4to3(curTrans->globalDerivs[dof] * dress);
					}
					else {
						daPose[pp++] += w * vec4to3(curTrans->globalDerivs[dof] * dress);
					}
				}
				// weight derivatives
				if (inf == 0) {
					for (i=0; i < maxInf-1; i++) {
						if (infJoints[pt*maxInf + i + 1] >= 0) {
							daWeight[(maxInf-1) * comp + i] -= vec4to3(curMats[joint] * dress);
						}
					}
				}
				else
					daWeight[(maxInf-1) * comp + inf - 1] += vec4to3(curMats[joint] * dress);
			}
			// dress derivatives
			daDress[0] += w * vec4to3(curMats[joint] * Vec4d(1, 0, 0, 0));
			daDress[1] += w * vec4to3(curMats[joint] * Vec4d(0, 1, 0, 0));
			daDress[2] += w * vec4to3(curMats[joint] * Vec4d(0, 0, 1, 0));


			trans += w * (curMats[joint]*Vec3d() - vec4to3(curMats[joint] * vec3to4z(dressJoints[joint])));

			// intrinsic and pose derivatives
			for (dof = 0; dof < skel->numDofs; dof++) {
				if (skel->dofIntrins[dof]) {
					daInt[ip++] += w *
						(curTrans->globalDerivs[dof]*Vec3d() - 
						vec4to3(curTrans->globalDerivs[dof] * vec3to4z(dressJoints[joint]) +
						curMats[joint] * curJointDerivs[joint*skel->numDofs + dof]));
				}
				else {
					daPose[pp++] += w *
						(curTrans->globalDerivs[dof]*Vec3d() - 
						vec4to3(curTrans->globalDerivs[dof] * vec3to4z(dressJoints[joint])));
				}
			}
			// weight derivatives
			if (inf == 0) {
				for (i=0; i < maxInf-1; i++) {
					if (infJoints[pt*maxInf + i + 1] >= 0) {
						daWeight[(maxInf-1) * numComponents + i] -= (curMats[joint]*Vec3d() - vec4to3(curMats[joint] * vec3to4z(dressJoints[joint])));
					}
				}
			}
			else
				daWeight[(maxInf-1) * numComponents + inf - 1] += (curMats[joint]*Vec3d() - vec4to3(curMats[joint] * vec3to4z(dressJoints[joint])));
		}
	}

	// actual calculation starts here ---------------------

	// first part: (v-t)^2
	ret = (v-trans) * (v-trans);
	if (intGrad) {
		for (ip = 0; ip < skel->numIntrinsicDofs; ip++) {
			intGrad[ip] += -2.0 * daInt[numComponents*skel->numIntrinsicDofs + ip] * (v-trans);
		}
	}
	if (poseGrad) {
		for (pp = 0; pp < skel->numPoseDofs; pp++) {
			poseGrad[pp] += -2.0 * daPose[numComponents*skel->numPoseDofs + pp] * (v-trans);
		}
	}
	if (weightGrad) {
		for (i=0; i < maxInf-1; i++) {
			weightGrad[i] += -2.0 * daWeight[(maxInf-1) * numComponents + i] * (v-trans);
		}
	}

	for (comp = 0; comp < numComponents; comp++) {
		// second part: -2 (v-t) A mu
		ret += -2.0 * (v - trans) * aMat[comp] * phi[0][comp];
		if (dressGrad) {
			dressGrad[(comp * numPts + pt) * 3 + 0] += -2.0 * ((v-trans) * daDress[0]) * phi[0][comp];
			dressGrad[(comp * numPts + pt) * 3 + 1] += -2.0 * ((v-trans) * daDress[1]) * phi[0][comp];
			dressGrad[(comp * numPts + pt) * 3 + 2] += -2.0 * ((v-trans) * daDress[2]) * phi[0][comp];
		}
		if (intGrad) {
			for (ip = 0; ip < skel->numIntrinsicDofs; ip++) {
				intGrad[ip] += -2.0 * (
					-daInt[numComponents*skel->numIntrinsicDofs + ip] * aMat[comp] +
					(v - trans) * daInt[comp*skel->numIntrinsicDofs + ip]
					) * phi[0][comp];
			}
		}
		if (poseGrad) {
			for (pp = 0; pp < skel->numPoseDofs; pp++) {
				poseGrad[pp] += -2.0 * (
					-daPose[numComponents*skel->numPoseDofs + pp] * aMat[comp] +
					(v - trans) * daPose[comp*skel->numPoseDofs + pp]
					) * phi[0][comp];
			}
		}
		if (weightGrad) {
			for (i=0; i < maxInf-1; i++) {
				weightGrad[i] += -2.0 * (
					-daWeight[(maxInf-1) * numComponents + i] * aMat[comp] +
					(v - trans) * daWeight[(maxInf-1) * comp + i]) * phi[0][comp];
			}
		}

		for (comp2 = 0; comp2 < numComponents; comp2++) {
			// third part: tr(A'A phi)
			ret += aMat[comp] * aMat[comp2] * phi[comp][comp2];
			if (dressGrad) {
				dressGrad[(comp * numPts + pt) * 3 + 0]  += daDress[0] * aMat[comp2] * phi[comp][comp2];
				dressGrad[(comp2 * numPts + pt) * 3 + 0] += daDress[0] * aMat[comp] * phi[comp][comp2];
				dressGrad[(comp * numPts + pt) * 3 + 1]  += daDress[1] * aMat[comp2] * phi[comp][comp2];
				dressGrad[(comp2 * numPts + pt) * 3 + 1] += daDress[1] * aMat[comp] * phi[comp][comp2];
				dressGrad[(comp * numPts + pt) * 3 + 2]  += daDress[2] * aMat[comp2] * phi[comp][comp2];
				dressGrad[(comp2 * numPts + pt) * 3 + 2] += daDress[2] * aMat[comp] * phi[comp][comp2];
			}
			if (intGrad) {
				for (ip = 0; ip < skel->numIntrinsicDofs; ip++) {
					intGrad[ip] += (daInt[comp*skel->numIntrinsicDofs + ip] * aMat[comp2] +
						aMat[comp] * daInt[comp2*skel->numIntrinsicDofs + ip]) *
						phi[comp][comp2];
				}
			}
			if (poseGrad) {
				for (pp = 0; pp < skel->numIntrinsicDofs; pp++) {
					poseGrad[pp] += (daPose[comp*skel->numPoseDofs + pp] * aMat[comp2] +
						aMat[comp] * daPose[comp2*skel->numPoseDofs + pp]) *
						phi[comp][comp2];
				}
			}
			if (weightGrad) {
				for (i=0; i < maxInf-1; i++) {
					weightGrad[i] += (
						aMat[comp] * daWeight[(maxInf-1) * comp + i] +
						daWeight[(maxInf-1) * comp2 + i] * aMat[comp2]) * phi[0][comp];
				}
			}
		}
	}

	return ret;
}

void USkin::calcGrad(int pt, Vec3d v, 
					 double *dressGrad, double *intGrad,
					 double *poseGrad, double *weightGrad, double *pddGrad) {
	int inf, inf2, samp;
	Vec3d dv;
	int joint, pddJoint;
	int firstPdd = ptPddOfs[pt*maxInf]*3;

	for (inf=0; inf < maxInf; inf++) {
		if (infJoints[pt*maxInf + inf] >= 0) {
			double w = infWeights[pt*maxInf + inf];
			joint = infJoints[pt*maxInf + inf];
			pddJoint = getPddJoint(pt, inf);
			Vec3d local = pddDressPts[pt] - dressJoints[joint];
			SkelTransform *curTrans = skel->transforms.getT(joint);

			if (pddGrad && transRBFs && transRBFs[pddJoint]) {
				RBF *rbf = transRBFs[pddJoint];

				// pose-dependent deformation variables
				for (inf2 = 0; inf2 < maxInf; inf2++) {
					double w2 = infWeights[pt*maxInf + inf2];
					int joint2 = infJoints[pt*maxInf + inf2];
					if (joint2 == -1 || w2 == 0)
						continue;

					for (samp = 0; samp < rbf->mN-1; samp++) {
						if (rbf->curWeights[samp+1] == 0) continue;

						pddGrad[ptPddOfs[pt*maxInf + inf]*3 + samp*3 + 0 - firstPdd] += 
							v * vec4to3(curMats[joint2] * 
							Vec4d(rbf->curWeights[samp+1]*w2, 0, 0, 0));
						pddGrad[ptPddOfs[pt*maxInf + inf]*3 + samp*3 + 1 - firstPdd] += 
							v * vec4to3(curMats[joint2] * 
							Vec4d(0, rbf->curWeights[samp+1]*w2, 0, 0));
						pddGrad[ptPddOfs[pt*maxInf + inf]*3 + samp*3 + 2 - firstPdd] += 
							v * vec4to3(curMats[joint2] * 
							Vec4d(0, 0, rbf->curWeights[samp+1]*w2, 0));
					}

					//// how joint angles affect this pdd
					//if (poseGrad) {
					//	int qInd;
					//	for (qInd = 0; qInd < 4; qInd++) {
					//		Vec3d deriv;
					//		for (samp = 0; samp < rbf->mD-1; samp++) {
					//			deriv += rbf->curDerivs[(samp+1)*rbf->mD + qInd] * 
					//				pddKeys[ptPddOfs[pt*maxInf + inf2]];
					//		}

					//		dv = vec4to3(curMats[joint2] * 
					//			Vec4d(deriv[0], deriv[1], deriv[2], 0));
					//		poseGrad[dofIndex[infJoints[pt*maxInf + inf]] + qInd] += 
					//			w2 * v * dv;
					//	}
					//}
				}
			}

			if (w != 0) {
				if (dressGrad) {
					// dressGrad
					dv = vec4to3(curMats[joint] * Vec4d(1, 0, 0, 0));
					dressGrad[0] += w * v * dv;
					dv = vec4to3(curMats[joint] * Vec4d(0, 1, 0, 0));
					dressGrad[1] += w * v * dv;
					dv = vec4to3(curMats[joint] * Vec4d(0, 0, 1, 0));
					dressGrad[2] += w * v * dv;
				}

				// pose-dependent gradients
				int dof = 0;
				int ip = 0;
				int pp = 0;
				for (dof = 0; dof < skel->numDofs; dof++) {
					if (skel->dofIntrins[dof]) {
						if (intGrad)
							intGrad[ip++] += w * v * (curTrans->globalDerivs[dof] * local - vec4to3(curMats[joint] * curJointDerivs[joint*skel->numDofs + dof]));
					}
					else {
						if (poseGrad)
							poseGrad[pp++] += w * v * (curTrans->globalDerivs[dof] * local);
					}
				}
			}

			if (weightGrad)
				weightGrad[inf] += v * (curMats[joint] * local);
		}
	}
}
*/