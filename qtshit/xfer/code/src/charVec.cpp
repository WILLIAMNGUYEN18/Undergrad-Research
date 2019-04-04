#include "charVec.h"
#include "skeleton.h"
#include "trimesh.h"
#include "trimesh_render.h"
#include "trimesh_util.h"
#include "ba.h"
#include "skinCalc.h"

extern Skin scSkin;

#define PIECEWISE_LINEAR
const int CV_DEGREE = 3;
const int CV_PIECES = 4;

// CharVecPt ==============================================

static double multipliers[1024];


CharVecPt::CharVecPt() {
	numInfluences = 0;
	influences = NULL;
	curComponents = NULL;
	cv = NULL;
}

CharVecPt::~CharVecPt() {
	if (influences)
		delete []influences;
	if (curComponents)
		delete []curComponents;
}

void CharVecPt::init(CharVec *iCV, int nInf) {
	cv = iCV;

	numInfluences = nInf;
	influences = new int[numInfluences];
	infWeights = new double[numInfluences];

	curComponents = new Vec3d[cv->numComponents];
	
#ifdef PIECEWISE_LINEAR
	dataParts = numInfluences * CV_PIECES;
#else
	dataParts = numInfluences * CV_DEGREE;
#endif

	data.resize(3 * cv->numComponents * dataParts);
	backup.resize(3 * cv->numComponents * dataParts);
}

void CharVecPt::initData() { 
	int index = 0;

#ifdef PIECEWISE_LINEAR
	int comp, inf, piece;
	Vec3d v;

	for (comp = 0; comp < cv->numComponents; comp++) {
		for (inf = 0; inf < numInfluences; inf++) {
			if (comp == 0 && inf == 0)
				v = localPos;
			else
				v = Vec3d();
//				v = Vec3d(boundedRand(-0.001, 0.001), boundedRand(-0.001, 0.001), boundedRand(-0.001, 0.001));

			for (piece = 0; piece < CV_PIECES; piece++) {
				data[index++] = v[0];
				data[index++] = v[1];
				data[index++] = v[2];
			}
		}
	}
	baAssert(index == data.size(), "size mismatch when initializing data!", true);
#else
	data[index++] = localPos[0];
	data[index++] = localPos[1];
	data[index++] = localPos[2];

	for (; index < CV_DEGREE*3; index++) {
		data[index] = 0;
	}

	for (; index < cv->numComponents*CV_DEGREE*3; index++) {
		// seed W with random values of +/- 1 mm
		data[index] = boundedRand(-0.001, 0.001);
	}
	for (; index < data.size(); index++) {
		data[index] = 0;
	}
#endif
}

double *CharVecPt::calcMultipliers(double *jointDofs) {
	int comp, inf, degree;
	int index = 0;

#ifdef PIECEWISE_LINEAR
	memset(multipliers, 0, sizeof(double) * dataParts);
	if (!cv->dofMin)
		return multipliers;
	for (inf = 0; inf < numInfluences; inf++) {
		int dofInd = influences[inf];
		double dof = jointDofs[dofInd];
		if (dof <= cv->dofMin[dofInd])
			multipliers[index + 0] = infWeights[inf];
		else if (dof >= cv->dofMax[dofInd])
			multipliers[index + CV_PIECES - 1] = infWeights[inf];
		else {
			double delta = cv->dofMax[dofInd] - cv->dofMin[dofInd];
			double pos = (CV_PIECES-1) * (dof - cv->dofMin[dofInd]) / delta;
			int first = pos;
			pos -= first;
			multipliers[index + first] = infWeights[inf] * (1.0 - pos);
			multipliers[index + first + 1] = infWeights[inf] * pos;
		}
		index += CV_PIECES;
	}
#else
	for (inf = 0; inf < numInfluences; inf++) {
		double dof = jointDofs[influences[inf]];
		double mult = dof;
		dof *= infWeights[inf];

		for (degree = 0; degree < CV_DEGREE; degree++) {
			multipliers[index] = dof;
			dof *= mult;
			index++;
		}
	}
#endif
	return multipliers;
}

void CharVecPt::updateCurComponents(double *jointDofs) {
	if (jointDofs)
		calcMultipliers(jointDofs);

	memset(curComponents, 0, sizeof(Vec3d) * cv->numComponents);

	int i, comp;
	int n = data.size() / 3 / cv->numComponents;
//	cout << n << endl;
	int index = 0;

	for (i=0; i < n; i++) {
		for (comp = 0; comp < cv->numComponents; comp++) {
			curComponents[comp] += multipliers[i] * Vec3d(data[index+0], data[index+1], data[index+2]);
			index += 3;
		}
	}
}

void CharVecPt::updateLocalPos(double *n, int nSize) {
	int comp;

	if (nSize < 0)
		nSize = cv->numComponents;
	if (nSize > cv->numComponents)
		nSize = cv->numComponents;
	if (!n)
		nSize = 1;

	localPos = curComponents[0];	// assume n[0] = 1
	if (!n) return;
	for (comp = 1; comp < nSize; comp++) {
		localPos += curComponents[comp] * n[comp];
	}
}

int CharVecPt::getInfIndex(int inf) {
	int i;
	for (i=0; i < numInfluences; i++)
		if (influences[i] == inf)
			return i;
	return -1;
}

/*
Vec3d CharVecPt::calcDeriv(int index) {
	switch (index) {
	case 0:
		return vec4to3(coord * Vec4d(1, 0, 0, 0));
		break;
	case 1:
		return vec4to3(coord * Vec4d(0, 1, 0, 0));
		break;
	default:
		return vec4to3(coord * Vec4d(0, 0, 1, 0));
		break;
	}

}
*/

// CharVec ================================================

void CharVec::createCylArm(int lSegs, int rSegs, double rad) {
/*	tm = new TriMesh();
	numPts = lSegs * rSegs;
	cvPts = new CharVecPt[numPts];
	tm->init(numPts, 2 * (lSegs - 1) * rSegs);

	vis = false;

	double elbowBlend = 0.2;
	double wristBlend = 0.2;

	int elbowPos = lSegs * 0.4;
	int wristPos = lSegs * 0.9;

	int shoulderQ = skel->transforms.lookupName("lShoulderQ");
	int upperArmT = skel->transforms.lookupName("lUpperArmT");
	int elbowA = skel->transforms.lookupName("lElbowA");
	int forearmA = skel->transforms.lookupName("lForearmA");
	int forearmT = skel->transforms.lookupName("lForearmT");
	int wristA = skel->transforms.lookupName("lWristA");
	int handT = skel->transforms.lookupName("lHandT");

	int x, y, i;
	int index = 0;
	for (x=0; x < lSegs; x++) {
		int joints[2];
		int jointMin = -lSegs * 0.1;
		int jointMax = lSegs * 0.4;

		if (x >= wristPos) {
			joints[0] = wristA;
			joints[1] = handT;
			jointMin = wristPos;
			jointMax = lSegs;
		}
		else if (x >= elbowPos) {
			joints[0] = elbowA;
			joints[1] = forearmT;
			jointMin = elbowPos;
			jointMax = wristPos;
		}
		else {
			joints[0] = shoulderQ;
			joints[1] = upperArmT;
		}

		double relPos = (1.0 * x - jointMin) / (1.0 * jointMax - jointMin);
		if (joints[0] == wristA)
			relPos /= 1.5;

		for (y=0; y < rSegs; y++) {
			cvPts[index].init(2);
			cvPts[index].pInfluences[0] = joints[0];
			cvPts[index].pInfluences[1] = joints[1];
			cvPts[index].skinPosWeights[0] = 1.0 - relPos;
			cvPts[index].skinPosWeights[1] = relPos;

			if (joints[0] == shoulderQ) {
				if (relPos > 1.0 - elbowBlend) {
					cvPts[index].rInfluences[0] = shoulderQ;
					cvPts[index].rInfluences[1] = elbowA;
					cvPts[index].skinRotWeights[0] = (relPos - (1.0 - elbowBlend)) * (0.5 / elbowBlend);
				}
				else {
					cvPts[index].rInfluences[0] = shoulderQ;
					cvPts[index].rInfluences[1] = shoulderQ;
					cvPts[index].skinRotWeights[0] = 1;
				}
			}
			else if (joints[0] == elbowA) {
				if (relPos < elbowBlend) {
					cvPts[index].rInfluences[0] = shoulderQ;
					cvPts[index].rInfluences[1] = elbowA;
					cvPts[index].skinRotWeights[0] = 0.5 + relPos * (0.5 / elbowBlend);
				}
				else if (relPos < 1.0 - wristBlend) {
					cvPts[index].rInfluences[0] = elbowA;
					cvPts[index].rInfluences[1] = forearmA;
					cvPts[index].skinRotWeights[0] = 1.0 - ((1.0 - wristBlend - relPos) / (1.0 - wristBlend - elbowBlend));
				}
				else {
					cvPts[index].rInfluences[0] = forearmA;
					cvPts[index].rInfluences[1] = wristA;
					cvPts[index].skinRotWeights[0] = (relPos - (1.0 - wristBlend)) * (0.5 / wristBlend);
				}
			}
			else {
				if (relPos < wristBlend * 1.5) {
					cvPts[index].rInfluences[0] = forearmA;
					cvPts[index].rInfluences[1] = wristA;
					cvPts[index].skinRotWeights[0] = 0.5 + relPos * (0.5 / (wristBlend*1.5));
				}
				else {
					cvPts[index].rInfluences[0] = handT;
					cvPts[index].rInfluences[1] = handT;
					cvPts[index].skinRotWeights[0] = 1.0;
				}
			}

			cvPts[index].axisRot = -2.0 * PI * y / rSegs;
			cvPts[index].pValues[0] = Vec3d(0, rad, 0); //Vec3d(0, sin(2.0 * PI * y / rSegs) * rad, cos(2.0 * PI * y / rSegs) * rad);
			for (i = 1; i < NUM_P_VALUES; i++)
				cvPts[index].pValues[i] = Vec3d();

			index++;
		}
	}

	// now let's make some triangles
	index = 0;
	for (x=0; x < lSegs-1; x++) {
		for (y=0; y < rSegs; y++) {
			tm->addTri(x*rSegs + y, (x+1)*rSegs + ((y+1) % rSegs), x*rSegs + ((y+1) % rSegs));
			tm->addTri(x*rSegs + y, (x+1)*rSegs + y, (x+1)*rSegs + ((y+1) % rSegs));
		}
	}
	*/
}

void CharVec::createFromSkin(Skin &skin, int iNumComponents) {
	int i;
	int numInf, inf[50];
	double infWeights[50];
	double w;

	numComponents = iNumComponents;

	cvPts = new CharVecPt[skin.numPts];
	numPts = skin.numPts;
	vis = false;

	numPreMirror = skin.numPreMirror;
	mirrorMap = skin.mirrorMap;

	int abdomenT = skin.skel->transforms.lookupName("abdomenT");
	int chestT = skin.skel->transforms.lookupName("chestT");
	int lClavicleT = skin.skel->transforms.lookupName("lClavicleT");
	int lUpperArmT = skin.skel->transforms.lookupName("lUpperArmT");
	int lForearm1T = skin.skel->transforms.lookupName("lForearm1T");
	int lForearm2T = skin.skel->transforms.lookupName("lForearm2T");
	int lHandT = skin.skel->transforms.lookupName("lHandT");
	int rClavicleT = skin.skel->transforms.lookupName("rClavicleT");
	int rUpperArmT = skin.skel->transforms.lookupName("rUpperArmT");
	int rForearm1T = skin.skel->transforms.lookupName("rForearm1T");
	int rForearm2T = skin.skel->transforms.lookupName("rForearm2T");
	int rHandT = skin.skel->transforms.lookupName("rHandT");
	int lPelvisT = skin.skel->transforms.lookupName("lPelvisT");
	int lThighT = skin.skel->transforms.lookupName("lThighT");
	int lShinT = skin.skel->transforms.lookupName("lShinT");
	int lFootT = skin.skel->transforms.lookupName("lFootT");
	int lToeT = skin.skel->transforms.lookupName("lToeT");
	int rPelvisT = skin.skel->transforms.lookupName("rPelvisT");
	int rThighT = skin.skel->transforms.lookupName("rThighT");
	int rShinT = skin.skel->transforms.lookupName("rShinT");
	int rFootT = skin.skel->transforms.lookupName("rFootT");
	int rToeT = skin.skel->transforms.lookupName("rToeT");
	int neckT = skin.skel->transforms.lookupName("neckT");

	const double MIN_WEIGHT = 0.0;;

	int index;
	for (index = 0; index < skin.numPts; index++) {
		numInf = 1;
		inf[0] = 0;
		infWeights[0] = 1;

//		if (index == 7683)
//			cout << "foo";
		// waistA
		w = skin.points[index].getTransWeight(abdomenT) + skin.points[index].getTransWeight(lPelvisT) + skin.points[index].getTransWeight(rPelvisT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 1;
			infWeights[numInf] = w; inf[numInf++] = 2;
			infWeights[numInf] = w; inf[numInf++] = 3;
		}
		// abdomenA
		w = skin.points[index].getTransWeight(abdomenT) + skin.points[index].getTransWeight(chestT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 4;
			infWeights[numInf] = w; inf[numInf++] = 5;
			infWeights[numInf] = w; inf[numInf++] = 6;
		}
		// lClavicleA
		w = skin.points[index].getTransWeight(lClavicleT) + skin.points[index].getTransWeight(chestT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 7;
			infWeights[numInf] = w; inf[numInf++] = 8;
			infWeights[numInf] = w; inf[numInf++] = 9;
		}
		// lShoulderA
		w = skin.points[index].getTransWeight(lClavicleT) + skin.points[index].getTransWeight(lUpperArmT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 10;
			infWeights[numInf] = w; inf[numInf++] = 11;
			infWeights[numInf] = w; inf[numInf++] = 12;
		}
		// lElbowA
		w = skin.points[index].getTransWeight(lForearm1T) + skin.points[index].getTransWeight(lUpperArmT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 13;
		}
		// lForearmA
		w = skin.points[index].getTransWeight(lForearm2T) + skin.points[index].getTransWeight(lForearm1T);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 14;
		}
		// lHandT
		w = skin.points[index].getTransWeight(lForearm2T) + skin.points[index].getTransWeight(lHandT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 15;
		}
		// rClavicleA
		w = skin.points[index].getTransWeight(chestT) + skin.points[index].getTransWeight(rClavicleT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 16;
			infWeights[numInf] = w; inf[numInf++] = 17;
			infWeights[numInf] = w; inf[numInf++] = 18;
		}
		// rShoulderA
		w = skin.points[index].getTransWeight(rClavicleT) + skin.points[index].getTransWeight(rUpperArmT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 19;
			infWeights[numInf] = w; inf[numInf++] = 20;
			infWeights[numInf] = w; inf[numInf++] = 21;
		}
		// rElbowA
		w = skin.points[index].getTransWeight(rForearm1T) + skin.points[index].getTransWeight(rUpperArmT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 22;
		}
		// rForearmA
		w = skin.points[index].getTransWeight(rForearm2T) + skin.points[index].getTransWeight(rForearm1T);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 23;
		}
		// rHandT
		w = skin.points[index].getTransWeight(rForearm2T) + skin.points[index].getTransWeight(rHandT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 24;
		}
		// lHipQ
		w = skin.points[index].getTransWeight(lPelvisT) + skin.points[index].getTransWeight(lThighT) + skin.points[index].getTransWeight(rPelvisT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 25;
			infWeights[numInf] = w; inf[numInf++] = 26;
			infWeights[numInf] = w; inf[numInf++] = 27;
		}
		// lKneeA and lShinA
		w = skin.points[index].getTransWeight(lThighT) + skin.points[index].getTransWeight(lShinT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 28;
			infWeights[numInf] = w; inf[numInf++] = 29;
		}
		// lAnkleA
		w = skin.points[index].getTransWeight(lShinT) + skin.points[index].getTransWeight(lFootT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 30;
		}
		// lFootA
		w = skin.points[index].getTransWeight(lFootT) + skin.points[index].getTransWeight(lToeT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 31;
		}
		// rHipQ
		w = skin.points[index].getTransWeight(rPelvisT) + skin.points[index].getTransWeight(rThighT) + skin.points[index].getTransWeight(lPelvisT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 32;
			infWeights[numInf] = w; inf[numInf++] = 33;
			infWeights[numInf] = w; inf[numInf++] = 34;
		}
		// rKneeA and rShinA
		w = skin.points[index].getTransWeight(rThighT) + skin.points[index].getTransWeight(rShinT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 35;
			infWeights[numInf] = w; inf[numInf++] = 36;
		}
		// rAnkleA
		w = skin.points[index].getTransWeight(rShinT) + skin.points[index].getTransWeight(rFootT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 37;
		}
		// rFootA
		w = skin.points[index].getTransWeight(rFootT) + skin.points[index].getTransWeight(rToeT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 38;
		}
		// neck
		w = skin.points[index].getTransWeight(chestT) + skin.points[index].getTransWeight(neckT);
		if (w > MIN_WEIGHT) {
			infWeights[numInf] = w; inf[numInf++] = 39;
			infWeights[numInf] = w; inf[numInf++] = 40;
			infWeights[numInf] = w; inf[numInf++] = 41;
		}
		
		cvPts[index].init(this, numInf);
		memcpy(cvPts[index].influences, inf, numInf*sizeof(int));
		memcpy(cvPts[index].infWeights, infWeights, numInf*sizeof(double));

#ifdef USE_POSITIONS
		cvPts[index].templateLocalPos = skin.points[index].localPos;
#else
		cvPts[index].templateLocalPos = Vec3d();
#endif
		cvPts[index].localPos = cvPts[index].templateLocalPos;

		cvPts[index].initData();
	}

	neighbors = findTMNeighbors(tm);
}

void CharVec::createFromMesh(TriMesh *mesh) {
	tm = mesh;
	numPts = tm->numPts();
	cvPts = new CharVecPt[numPts];
	vis = false;

	int index;
	for (index = 0; index < numPts; index++) {
		cvPts[index].init(this, 0);
	}

	neighbors = findTMNeighbors(tm);
}

void CharVec::setFromLocalTM(Mat4d *localIFrames) {
	int i;
	for (i=0; i < numPts; i++) {
		cvPts[i].localPos = localIFrames[i] * tm->getPt(i);
	}
}

void CharVec::backupData() {
	int pt;
	for (pt=0; pt < numPts; pt++)
		cvPts[pt].backup = cvPts[pt].data;
}

void CharVec::updateCurComponents(double *jointDofs) {
	int pt;
	for (pt=0; pt < numPts; pt++)
		cvPts[pt].updateCurComponents(jointDofs);
}

void CharVec::updateLocalPos(double *n, int nSize) {
	int pt;
	for (pt=0; pt < numPts; pt++)
		cvPts[pt].updateLocalPos(n, nSize);
}

void CharVec::updateTM(Mat4d *localFrames) {
	int pt;
	Vec3d curPos;

	baAssert(numPts == tm->numPts(), "mismatched number of points updating tm", false);

	for (pt=0; pt < numPts; pt++) {
		if (localFrames)
			curPos = localFrames[pt] * cvPts[pt].localPos;
		else
#ifdef USE_POSITIONS
			curPos = scSkin.points[pt].localMat * cvPts[pt].localPos;
#else
			curPos = Vec3d();
#endif
		tm->getPt(pt) = curPos;

//		if (tm->gsPtColors()) {
//			tm->getPtColor(pt) = Vec3d(0.8, 0.8, 0.8);
//		}
	}
	tm->calcNormals();
}

void CharVec::render(int viewMode, Vec3d bkg) {
	if (vis && tm) {
		renderTriMesh(tm, viewMode, bkg);
//		skel->drawGL();
	}
}

void CharVec::dofsFromSkel(Skeleton *skel, double *dofs) {
	int index = 0;
	Vec3d euler;

	dofs[index++] = 1;

	// torso
	euler = matToEuler(((SkelQuatRotation*)skel->transforms.getT("waistQ"))->curCoord.mat);
	dofs[index++] = euler[0];
	dofs[index++] = euler[1];
	dofs[index++] = euler[2];
	euler = matToEuler(((SkelQuatRotation*)skel->transforms.getT("abdomenQ"))->curCoord.mat);
	dofs[index++] = euler[0];
	dofs[index++] = euler[1];
	dofs[index++] = euler[2];

	// left arm
	euler = matToEuler(((SkelQuatRotation*)skel->transforms.getT("lClavicleQ"))->curCoord.mat);
	dofs[index++] = euler[0];
	dofs[index++] = euler[1];
	dofs[index++] = euler[2];
	euler = matToEuler(((SkelQuatRotation*)skel->transforms.getT("lShoulderQ"))->curCoord.mat);
	dofs[index++] = euler[0];
	dofs[index++] = euler[1];
	dofs[index++] = euler[2];
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("lElbowA"))->curAngle;
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("lForearmA"))->curAngle;
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("lWristA"))->curAngle;

	// right arm
	euler = matToEuler(((SkelQuatRotation*)skel->transforms.getT("rClavicleQ"))->curCoord.mat);
	dofs[index++] = euler[0];
	dofs[index++] = euler[1];
	dofs[index++] = euler[2];
	euler = matToEuler(((SkelQuatRotation*)skel->transforms.getT("rShoulderQ"))->curCoord.mat);
	dofs[index++] = euler[0];
	dofs[index++] = euler[1];
	dofs[index++] = euler[2];
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("rElbowA"))->curAngle;
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("rForearmA"))->curAngle;
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("rWristA"))->curAngle;

	// left leg
	euler = matToEuler(((SkelQuatRotation*)skel->transforms.getT("lHipQ"))->curCoord.mat);
	dofs[index++] = euler[0];
	dofs[index++] = euler[1];
	dofs[index++] = euler[2];
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("lKneeA"))->curAngle;
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("lShinA"))->curAngle;
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("lAnkleA"))->curAngle;
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("lFootA"))->curAngle;

	// right leg
	euler = matToEuler(((SkelQuatRotation*)skel->transforms.getT("rHipQ"))->curCoord.mat);
	dofs[index++] = euler[0];
	dofs[index++] = euler[1];
	dofs[index++] = euler[2];
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("rKneeA"))->curAngle;
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("rShinA"))->curAngle;
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("rAnkleA"))->curAngle;
//	dofs[index++] = ((SkelEulerRotation*)skel->transforms.getT("rFootA"))->curAngle;

	// neck
	euler = matToEuler(((SkelQuatRotation*)skel->transforms.getT("neckR"))->curCoord.mat);
	dofs[index++] = euler[0];
	dofs[index++] = euler[1];
	dofs[index++] = euler[2];
}

void swapNeg(double &a, double &b) {
	double temp = a;
	a = -b;
	b = -temp;
}

void CharVec::mirrorDofs(double *dofs) {
	int index = 0;
	Vec3d euler;

	dofs[index++] = 1;

	// torso
	dofs[index++] *= -1;
	dofs[index++] *= 1;
	dofs[index++] *= -1;
	dofs[index++] *= -1;
	dofs[index++] *= 1;
	dofs[index++] *= -1;

	// left arm
	// right arm
	swapNeg(dofs[index+0], dofs[index+9]);
	swap(dofs[index+1], dofs[index+10]);
	swapNeg(dofs[index+2], dofs[index+11]);
	swapNeg(dofs[index+3], dofs[index+12]);
	swap(dofs[index+4], dofs[index+13]);
	swapNeg(dofs[index+5], dofs[index+14]);
	swapNeg(dofs[index+6], dofs[index+15]);
	swapNeg(dofs[index+7], dofs[index+16]);
	swap(dofs[index+8], dofs[index+17]);
	index += 18;

	// left leg
	// right leg
	swap(dofs[index+0], dofs[index+8]);
	swapNeg(dofs[index+1], dofs[index+9]);
	swap(dofs[index+2], dofs[index+10]);
	swapNeg(dofs[index+3], dofs[index+11]);
	swap(dofs[index+4], dofs[index+12]);
	swapNeg(dofs[index+5], dofs[index+13]);
	swap(dofs[index+6], dofs[index+14]);
	swap(dofs[index+7], dofs[index+15]);
	index += 16;

	// neck
	dofs[index++] *= -1;
	dofs[index++] *= 1;
	dofs[index++] *= -1;
}