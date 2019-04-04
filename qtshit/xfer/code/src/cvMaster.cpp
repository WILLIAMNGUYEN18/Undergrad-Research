#include <fstream>
using namespace std;

#include "doppel2.h"
#include <GL/glu.h>
#include "doppel2.h"
#include "cvMaster.h"
#include "cli.h"
#include "dof.h"
#include "skeleton.h"
#include "charVec.h"
#include "trimesh.h"
#include "cvgf.h"
#include "skinCalc.h"
#include "surfaceDef.h"
#include "solver.h"
#include "jactest.h"
#include "meshmap.h"
#include "mocap.h"
#include "main_win.h"
#include "vl/VLd.h"

extern Skin scSkin;
extern TriMesh *scMesh;

CharVec charVec;

//int numMatches, numChars;
Skeleton *matchPoses;
TriMesh *matchMeshes;
//int *charMap;
CVGoalFunction *cvGF = NULL;
CVSkinningGF *cvSGF = NULL;
//SkinMatchGF *skinMatchGF = NULL;

static LBFGSSolver *solver = NULL;
LADeformationGoalFunction *matchGF = NULL;
EdgeMatchGF *edgeMatchGF = NULL;

MocapData mcData;

double brushInt = 1;
double brushRad = 0.02;

Mat4d transformFromPts(vector<Vec3d> &orig, vector<Vec3d> &target) {
	int numPts = (int)orig.size();
	int i;

	int r = 3, c = numPts;
	TMat oMat(r, c), tMat(3, c);
	TMat U(r, c), V(c, c);
	TVec diagonal(c);

	Vec3d oAvg, tAvg;
	for (i=0; i < numPts; i++) {
		oAvg += orig[i];
		tAvg += target[i];
	}
	oAvg /= numPts;
	tAvg /= numPts;

	for (i=0; i < numPts; i++) {
		oMat.Elt(0, i) = orig[i][0] - oAvg[0];
		oMat.Elt(1, i) = orig[i][1] - oAvg[1];
		oMat.Elt(2, i) = orig[i][2] - oAvg[2];
		tMat.Elt(0, i) = target[i][0] - oAvg[0];
		tMat.Elt(1, i) = target[i][1] - oAvg[1];
		tMat.Elt(2, i) = target[i][2] - oAvg[2];
	}

	// calculate pseudoinverse
	SVDFactorization(oMat, U, V, diagonal);
	TMat pinv(c, r), diag(c, c);
	diag = vl_0;
	for (i = 0; i < c; i++) {
		if (fabs(diagonal[i]) > 0.00001)
			diag.Elt(i, i) = 1.0 / diagonal[i];
	}
	pinv = V * diag * trans(U);

	TMat affine(3, 3);
	affine = tMat * pinv;

/*	if (fabs(affine.Elt(0,0)) > 10) {
		cout << oMat << endl << tMat << endl;
		cout << affine << endl;
		exit(0);
	}*/

	Mat4d m = Mat4d(0,0,0,oAvg[0],0,0,0,oAvg[1],0,0,0,oAvg[2],0,0,0,1);
	Mat4d ret = Mat4d(
		affine.Elt(0,0), affine.Elt(0,1), affine.Elt(0,2), 0,
		affine.Elt(1,0), affine.Elt(1,1), affine.Elt(1,2), 0,
		affine.Elt(2,0), affine.Elt(2,1), affine.Elt(2,2), 0,
		0, 0, 0, 1);
	ret = (m) * ret * (-m);
//	ret[0][3] += tAvg[0] - oAvg[0];
//	ret[1][3] += tAvg[1] - oAvg[1];
//	ret[2][3] += tAvg[2] - oAvg[2];
	return ret;
/*	return Mat4d(
		affine.Elt(0,0), affine.Elt(0,1), affine.Elt(0,2), tAvg[0] - oAvg[0],
		affine.Elt(1,0), affine.Elt(1,1), affine.Elt(1,2), tAvg[1] - oAvg[1],
		affine.Elt(2,0), affine.Elt(2,1), affine.Elt(2,2), tAvg[2] - oAvg[2],
		0, 0, 0, 1);*/
}

void cvBrushPt(int x, int y) {
	Vec3d orig, direction;
	y = mainWin->viewer->viewPort[3] - y - 1;
	gluUnProject(GLdouble(x), GLdouble(y), GLdouble(0), mainWin->viewer->modelMatrix, 
		mainWin->viewer->projMatrix, mainWin->viewer->viewPort, &orig[0], &orig[1], &orig[2]);
	gluUnProject(GLdouble(x), GLdouble(y), GLdouble(1), mainWin->viewer->modelMatrix, 
		mainWin->viewer->projMatrix, mainWin->viewer->viewPort, &direction[0], &direction[1], &direction[2]);
	direction = -(direction - orig);
	direction.normalize();

	if (cvGF->cv->tm->calcRayIntersection(orig, direction)){
		//if it hits the positive side
		if(cvGF->cv->tm->hitNeg){
			double time = cvGF->cv->tm->tNeg;
			Vec3d intersectPt = direction * time + orig;
			double nearestDist = -1;
			Vec3d nearestPt;
			int selectedMeshIndex;
			//Check each point in the mesh to see if it is closer to the intersection than any others seen so far.
			for(int i = 0; i < cvGF->cv->tm->numPts(); i++){
				double dist = (cvGF->cv->tm->getPt(i) - intersectPt).length();
				//calculate absolute value of dist.
				if(dist < 0){
					dist = -dist;
				}
				if (dist < brushRad) {
//					double c = (1.0 - brushInt) + (dist / brushRad) * brushInt;
					double c = 1.0 - brushInt;
					Vec3d &v = cvGF->cv->tm->getPtColor(i);
					c = min(c, v[1]);
					v = Vec3d(1, c, c);
					if (cvGF->cv->mirrorMap && cvGF->cv->mirrorMap[i] > -1)
						cvGF->cv->tm->getPtColor(cvGF->cv->mirrorMap[i]) = v;
				}
				if(dist < nearestDist || nearestDist < 0){
					nearestDist = dist;
					nearestPt = cvGF->cv->tm->getPt(i);
					selectedMeshIndex = i;
				}
			}
			cout << selectedMeshIndex << endl;
		}
	}
	redrawV();
}

void cvClick(int x, int y, int button) {
//	if (!cvGF || !cvGF->cv->tm)
//		return;
//	cvBrushPt(x, y);
}

void cvDrag(int x, int y) {
	if (!cvGF || !cvGF->cv->tm)
		return;
	cvBrushPt(x, y);
}

bool cvStartDrag(GLuint *nameBuf, int x, int y) {
	if (!cvGF || !cvGF->cv->tm)
		return false;
	cvBrushPt(x, y);
	return true;
}

void cvBrush(const char *params) {
	params = extractDouble(params, &brushInt);
	params = extractDouble(params, &brushRad);
}

void cvSaveColor(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);

	FILE *f;
	int i;

	if (!openFile(&f, fname, "wb", "texture")) {
		return;
	}

	char ver = '0';
	fwrite(&ver, sizeof(char), 1, f);

	i = cvGF->cv->tm->numPts() * 3;
	fwrite(&i, sizeof(int), 1, f);
	for (i=0; i < cvGF->cv->tm->numPts(); i++) {
		fwrite(cvGF->cv->tm->getPtColor(i).n, sizeof(double), 3, f);
	}
	fclose(f);

	/* ascii version
	ofstream os;
	int i;

	if (!openOFStream(&os, fname, "texture")) {
		return;
	}

	os << (curMesh->numEvalPts * 3) << endl;
	for (i=0; i < curMesh->numEvalPts; i++) {
		os << curMesh->evalPts[i].color[0] << endl;
		os << curMesh->evalPts[i].color[1] << endl;
		os << curMesh->evalPts[i].color[2] << endl;
	}
	os.close();
	*/
}

void cvLoadColor(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);

	FILE *f;
	if (!openFile(&f, fname, "rb", "points"))
		return;

	char version;
	fread(&version, sizeof(char), 1, f);

	if (version == '0') {
		int i;
		fread(&i, sizeof(int), 1, f);

		if (i != cvGF->cv->tm->numPts() * 3) {
			cout << "point size mismatch; expected " << cvGF->cv->tm->numPts()*3 << "; read " << i << endl;
		}
		else {
			for (i=0; i < cvGF->cv->tm->numPts(); i++) {
				fread(cvGF->cv->tm->getPtColor(i).n, sizeof(double), 3, f);
			}
		}
		fclose(f);
	}
	else {
		// text version
		fclose(f);
		cout << "loading text format" << endl;

		ifstream is;
		if (!openIFStream(&is, fname, "texture")) {
			return;
		}
		int i;
		is >> i;

		if (i != cvGF->cv->tm->numPts() * 3) {
			cout << "size mismatch in deformation" << endl;
			is.close();
			return;
		}
		for (i=0; i < cvGF->cv->tm->numPts(); i++) {
			is >> cvGF->cv->tm->getPtColor(i);
		}
		is.close();
	}
}

void cvLoadNeighWeights(const char *params) {
	Vec3d v;
	char fname[80];
	params = extractString(params, fname, 80);

	FILE *f;
	if (!openFile(&f, fname, "rb", "points"))
		return;

	char version;
	fread(&version, sizeof(char), 1, f);

	if (version == '0') {
		int i;
		fread(&i, sizeof(int), 1, f);

		if (i != cvGF->cv->numPts * 3) {
			cout << "point size mismatch; expected " << cvGF->cv->numPts*3 << "; read " << i << endl;
		}
		else {
			for (i=0; i < cvGF->cv->numPts; i++) {
				fread(v.n, sizeof(double), 3, f);
				if (matchGF)
					matchGF->neighWeights[i] = 0.2 + 0.8 * v[1];
				if (edgeMatchGF)
					edgeMatchGF->neighWeights[i] = 0.1 + 0.9 * v[1];
			}
		}
		fclose(f);
	}
	else
		cout << "format not supported" << endl;
}

void cvLoadSurfWeights(const char *params) {
	Vec3d v;
	char fname[80];
	params = extractString(params, fname, 80);

	FILE *f;
	if (!openFile(&f, fname, "rb", "points"))
		return;

	char version;
	fread(&version, sizeof(char), 1, f);

	if (version == '0') {
		int i;
		fread(&i, sizeof(int), 1, f);

		if (i != cvGF->cv->numPts * 3) {
			cout << "point size mismatch; expected " << cvGF->cv->numPts*3 << "; read " << i << endl;
		}
		else {
			for (i=0; i < cvGF->cv->numPts; i++) {
				fread(v.n, sizeof(double), 3, f);
				if (matchGF)
					matchGF->surfWeights[i] = v[1];
				if (edgeMatchGF)
					edgeMatchGF->surfWeights[i] = v[1];
			}
		}
		fclose(f);
	}
	else
		cout << "format not supported" << endl;
}

void cvInit(const char *params) {
	cvGF = new CVGoalFunction();
	cvGF->cv = &charVec;

	charVec.vis = true;
//	charVec.skel = Skeleton::load("data/body.sk.txt");
	charVec.skel = Skeleton::load("data/james.sk.txt");

//	charVec.createCylArm(40, 20, 0.05);
	charVec.createCylArm(80, 40, 0.05);
//	charVec.createCylArm(100, 50, 0.05);
	charVec.updateTM();
	charVec.vis = true;
}

void cvInitFromSkin(const char *params) {
	int numComponents;
	params = extractInt(params, &numComponents);

	cvGF = new CVGoalFunction();
	cvGF->cv = &charVec;

	charVec.vis = true;
	charVec.skel = scSkin.skel; //Skeleton::load("data/body.sk.txt");

	charVec.tm = scMesh;
	charVec.createFromSkin(scSkin, numComponents);
//	charVec.updateTM();
	charVec.vis = true;
}

/*
void cvInitFromMesh(const char *params) {
	char skelFName[256], meshFName[256];
	params = extractString(params, skelFName, 256);
	params = extractString(params, meshFName, 256);

	cvGF = new CVGoalFunction();
	cvGF->cv = &charVec;

	charVec.vis = true;
	charVec.skel = Skeleton::load(skelFName); //"data/body.sk.txt");

	charVec.tm = new TriMesh();
	charVec.tm->loadFile(meshFName);
	charVec.createFromSkin(scSkin);
//	charVec.updateTM();
	charVec.vis = true;
}*/

void cvAutoUpdate(const char *params) {
	params = extractBool(params, &charVec.update);
}

void cvBuildExampleSet(const char *params) {
	char setName[80], temp[80];
	params = extractString(params, setName, 80);
	ifstream setIn, in;
	int i;

	if (!openIFStream(&setIn, setName, "example set"))
		return;

	setIn >> cvGF->numExamples >> cvGF->numChars;
	matchPoses = new Skeleton[cvGF->numExamples];
	matchMeshes = new TriMesh[cvGF->numExamples];
	cvGF->examples = new CVExample[cvGF->numExamples];

	for (i=0; i < cvGF->numExamples; i++) {
		setIn >> cvGF->examples[i].charID;

		matchPoses[i].copyFrom(charVec.skel);
		setIn >> temp;
		if (!openIFStream(&in, temp, "pose"))
			return;
		matchPoses[i].loadPose(in);
		in.close();

		setIn >> temp;
		matchMeshes[i].loadPly(temp);
		matchMeshes[i].calcNormals();
		matchMeshes[i].showColor = false;
		matchMeshes[i].solidColor = Vec3d(0.8, 0.8, 0.8);
		matchMeshes[i].alpha = 0.8;

		scSkin.skel->copyVals(&matchPoses[i]);
		scSkin.skel->updateCoords();
		scSkin.updatePoints();
		cvGF->examples[i].buildExample(&matchMeshes[i], &scSkin);

		setIn >> cvGF->examples[i].fname;
//		cvGF->examples[i].save(temp);
	}

	cvGF->init(cvGF->numExamples, cvGF->numChars);

	cvGF->cv->updateCurComponents(cvGF->examples[cvGF->numExamples-1].dofs);
//	cvGF->cv->updateLocalPos(NULL);
//	cvGF->cv->updateTM();
}

static void swapV(Vec3d &v0, Vec3d &v1) {
	Vec3d temp = v0;
	v0 = v1;
	v1 = temp;
}

static void swapM(Mat4d &v0, Mat4d &v1) {
	Mat4d temp = v0;
	v0 = v1;
	v1 = temp;
}

static void swapQ(QuatNorm &v0, QuatNorm &v1) {
	QuatNorm temp = v0;
	v0 = v1;
	v1 = temp;
}

void reflectPose(Skeleton *skel) {
	int i;
	char mName[255];
	SkelTransform *trans, *mTrans;

	for (i=0; i < skel->transforms.size(); i++) {
		trans = skel->transforms.getT(i);

		if (strcmp(trans->name, "baseT") == 0) {
			((SkelTranslation*)trans)->curVal[1] *= -1;
		} 
		else if (strcmp(trans->className, "SkelQuatRotation") == 0) {
			QuatNorm &q = ((SkelQuatRotation*)trans)->curQuat;
			q.w *= -1;
			q.y *= -1;
		}
		else if (strcmp(trans->className, "SkelEulerRotation") == 0) {
			if (((SkelEulerRotation*)trans)->axis != 1)
				((SkelEulerRotation*)trans)->curAngle *= -1;
		}
	}

	for (i=0; i < skel->transforms.size(); i++) {
		trans = skel->transforms.getT(i); 
		if (trans->name[0] == 'l' && trans->name[1] != 'e') {
			strcpy(mName, trans->name);
			mName[0] = 'r';
			mTrans = skel->transforms.getT(mName);
			if (!mTrans) {
				cout << "warning: can't find " << mName << endl;
				continue;
			}

			if (strcmp(trans->className, "SkelQuatRotation") == 0) {
				swapQ(((SkelQuatRotation*)trans)->curQuat, ((SkelQuatRotation*)mTrans)->curQuat);
			}
			else if (strcmp(trans->className, "SkelEulerRotation") == 0) {
				swap(((SkelEulerRotation*)trans)->curAngle, ((SkelEulerRotation*)mTrans)->curAngle);
			}
		}
	}
}

void cvLoadExampleSet(const char *params) {
	char setName[80], temp[80];
	params = extractString(params, setName, 80);
	bool mirrorSet = false;
	params = extractBool(params, &mirrorSet);
	bool localCoords = true;
	params = extractBool(params, &localCoords);
	ifstream setIn, in;
	int i, j;

	if (!openIFStream(&setIn, setName, "example set"))
		return;

	setIn >> cvGF->numExamples >> cvGF->numChars;
	if (mirrorSet)
		cvGF->numExamples *= 2;
	matchPoses = new Skeleton[cvGF->numExamples];
	cvGF->examples = new CVExample[cvGF->numExamples];

	for (i=0; i < cvGF->numExamples; i++) {
		setIn >> cvGF->examples[i].charID;

		matchPoses[i].copyFrom(charVec.skel);
		setIn >> temp;
		if (!openIFStream(&in, temp, "pose"))
			return;
		matchPoses[i].loadPose(in);
		in.close();

		setIn >> temp;

		setIn >> cvGF->examples[i].fname;
		cvGF->examples[i].load(cvGF->examples[i].fname);

		for (j=0; j < cvGF->examples[i].numPoints; j++)
			cvGF->examples[i].conf[j] = min(1.0, 0.01+0.99*cvGF->examples[i].conf[j]);

		if (cvGF->cv->mirrorMap && mirrorSet) {
			int j, n;

			i++;

			matchPoses[i].copyFrom(&matchPoses[i-1]);
			reflectPose(&matchPoses[i]);
			matchPoses[i].updateCoords();
			scSkin.skel->copyVals(&matchPoses[i]);
			scSkin.skel->updateCoords();
			scSkin.updatePoints();

			cvGF->examples[i].charID = cvGF->examples[i-1].charID;
			sprintf(cvGF->examples[i].fname, "%sm", cvGF->examples[i-1].fname);
			n = cvGF->examples[i-1].numPoints;
			cvGF->examples[i].init(n);
			for (j=0; j < n; j++) {
				cvGF->examples[i].points[j] = cvGF->examples[i-1].points[j];
				cvGF->examples[i].points[j][1] *= -1;
				cvGF->examples[i].conf[j] = cvGF->examples[i-1].conf[j];
#ifdef USE_POSITIONS
				cvGF->examples[i].trans[j] = scSkin.points[j].localMat;
#endif
				cvGF->examples[i].iTrans[j] = cvGF->examples[i].trans[j].inverse();
			}

/*				cvGF->examples[i].trans[j] = cvGF->examples[i-1].trans[j];
				cvGF->examples[i].trans[j][1][0] *= -1;
//				cvGF->examples[i].trans[j][1][1] *= -1;
				cvGF->examples[i].trans[j][1][2] *= -1;
				cvGF->examples[i].trans[j][0][1] *= -1;
				cvGF->examples[i].trans[j][0][2] *= -1;
				cvGF->examples[i].trans[j][1][3] *= -1;
				cvGF->examples[i].iTrans[j] = cvGF->examples[i].trans[j].inverse();
			}
*/			// now swap!
			for (j = 0; j < charVec.numPreMirror; j++) {
				if (charVec.mirrorMap[j] != j) {
					swapV(cvGF->examples[i].points[j], cvGF->examples[i].points[charVec.mirrorMap[j]]);
					swap(cvGF->examples[i].conf[j], cvGF->examples[i].conf[charVec.mirrorMap[j]]);
//					swapM(cvGF->examples[i].trans[j], cvGF->examples[i].trans[charVec.mirrorMap[j]]);
//					swapM(cvGF->examples[i].iTrans[j], cvGF->examples[i].iTrans[charVec.mirrorMap[j]]);
				}
			}
			cvGF->examples[i].numDofs = cvGF->examples[i-1].numDofs;
			cvGF->examples[i].dofs = new double[cvGF->examples[i].numDofs];
			CharVec::dofsFromSkel(&matchPoses[i], cvGF->examples[i].dofs);
//			memcpy(cvGF->examples[i].dofs, cvGF->examples[i-1].dofs, sizeof(double)*cvGF->examples[i].numDofs);
//			CharVec::mirrorDofs(cvGF->examples[i].dofs);
		}

		if (localCoords) {
			if (cvGF->examples[i].numPoints == cvGF->examples[i].numTPoints) {
				// put point in local coords
				int j;
				for (j=0; j < cvGF->examples[i].numPoints; j++)
					cvGF->examples[i].points[j] = cvGF->examples[i].iTrans[j] * cvGF->examples[i].points[j];
			}
			if (mirrorSet) {
				if (cvGF->examples[i-1].numPoints == cvGF->examples[i-1].numTPoints) {
					// put point in local coords
					int j;
					for (j=0; j < cvGF->examples[i-1].numPoints; j++)
						cvGF->examples[i-1].points[j] = cvGF->examples[i-1].iTrans[j] * cvGF->examples[i-1].points[j];
				}

			}
		}
	}

	cvGF->init(cvGF->numExamples, cvGF->numChars);

	cvGF->cv->updateCurComponents(cvGF->examples[cvGF->numExamples-1].dofs);
	cvGF->cv->updateTM(cvGF->examples[i-1].trans);

	// calculate min/max
	int numDofs = cvGF->examples[0].numDofs;
	charVec.dofMin = new double[numDofs];
	charVec.dofMax = new double[numDofs];
	for (j=0; j < numDofs; j++) {
		charVec.dofMin[j] = 1e6;
		charVec.dofMax[j] = -1e6;
	}
	for (i=0; i < cvGF->numExamples; i++) {
		for (j=0; j < numDofs; j++) {
			if (cvGF->examples[i].dofs[j] < charVec.dofMin[j])
				charVec.dofMin[j] = cvGF->examples[i].dofs[j];
			if (cvGF->examples[i].dofs[j] > charVec.dofMax[j])
				charVec.dofMax[j] = cvGF->examples[i].dofs[j];
		}
	}
}

//void cvInitSolve(const char *params) {

//}

void cvRunEM(const char *params) {
	int numIterations = 5;
	params = extractInt(params, &numIterations);

	cvGF->runEM(numIterations);
}

void cvStopEM(const char *params) {
	if (!cvGF) {
		cout << "no goal function created" << endl;
		return;
	}
	if (!cvGF->solver) {
		cout << "solver not running" << endl;
		return;
	}

	cvGF->solver->stopNow = true;
}

void cvElbowAnim(const char *params) {
	char fname[80], str[80];
	str[0] = 0;
	params = extractString(params, str, 80);
	if (str[0] == 0)
		strcpy(str, "elbow");

	double *elbowA = &(((SkelPolarAxisRotation*)charVec.skel->transforms.getT("lElbowA"))->curTheta);

	int i;
	for (i=0; i < 60; i++) {
		if (i <= 30)
			*elbowA = -i * 2.3 / 30.0;
		else
			*elbowA = -(60-i) * 2.3 / 30.0;

		sprintf(fname, "anim/%s-%03d.tga", str, i);
		redrawVNow();
		uiScreenshot(fname);
	}
}

void cvCalcRMS(const char *params) {
/*	double err = 0;
	int count = 0;

	int ex, pt;
	for (ex=0; ex < cvGF->numExamples; ex++) {
		charVec.skel->copyVals(&cvGF->poses[ex]);
		charVec.skel->updateCoords();
		charVec.updateTM();

		charVec.tm->calcHBB(100);

		for (pt=0; pt < cvGF->meshes[ex].numPts(); pt++) {
			double dist = 0;

			charVec.tm->closestRestrictNormal = false;
			if (!charVec.tm->calcClosestPoint(cvGF->meshes[ex].getPt(pt), 1.0)) {
				cout << "warning: no closest point " << pt << endl;
			}
			else {
				count++;
				err += sqr(charVec.tm->closestDist);
			}
		}
	}

	cout << "RMS error = " << sqrt(err / count) << endl;*/
}

void cvSave(const char *params) {
	char fname[80];
	fname[0] = 0;
	params = extractString(params, fname, 80);

	FILE *f;
	if (!openFile(&f, fname, "wb", "charvec file")) {
		return;
	}

	int i, n;

	i = charVec.numPts;
	fwrite(&i, sizeof(int), 1, f);
	i = charVec.numComponents;
	fwrite(&i, sizeof(int), 1, f);
//	i = CV_DEGREE;
//	fwrite(&i, sizeof(int), 1, f);

	for (i=0; i < charVec.numPts; i++) {
		n = charVec.cvPts[i].numInfluences;
		fwrite(&n, sizeof(int), 1, f);
		fwrite(charVec.cvPts[i].influences, sizeof(int), n, f);
		n = charVec.cvPts[i].data.size();
		fwrite(&n, sizeof(int), 1, f);
		fwrite(charVec.cvPts[i].data.n, sizeof(double), n, f);
	}

	fclose(f);
}

void cvLoad(const char *params) {
	char fname[80], mapFN[80];
	fname[0] = 0;
	mapFN[0] = 0;
	params = extractString(params, fname, 80);
	params = extractString(params, mapFN, 80);

	FILE *f;
	if (!openFile(&f, fname, "rb", "charvec file")) {
		return;
	}

	int i, j, k, x, n, loadPts;

	fread(&loadPts, sizeof(int), 1, f);
	fread(&i, sizeof(int), 1, f); // numComponents
//	i = CV_DEGREE;
//	fread(&i, sizeof(int), 1, f); // degree

	if (loadPts == charVec.numPts) {
		for (i=0; i < charVec.numPts; i++) {
			fread(&n, sizeof(int), 1, f);
			if (n != charVec.cvPts[i].numInfluences) {
				cout << "# of influences mismatch on point " << i << "; aborting!" << endl;
				return;
			}
			fread(charVec.cvPts[i].influences, sizeof(int), n, f);
			fread(&n, sizeof(int), 1, f);
			if (n != charVec.cvPts[i].data.size()) {
				cout << "# of data points mismatch on point " << i << "; aborting!" << endl;
				return;
			}
			fread(charVec.cvPts[i].data.n, sizeof(double), n, f);
		}
	}
	else if (mapFN[0] == '0') {
		cout << "point count mismatch; expected " << charVec.numPts << "; found " << i << endl;
		return;
	}
	else {
		cout << " mapping not fixed yet!!" << endl;
		/*
		// open map file
		FILE *mapF = NULL;
		int mapSmall, mapBig;
		if (!openFile(&mapF, mapFN, "rb", "mapping"))
			return;
		fread(&mapSmall, sizeof(int), 1, mapF);
		fread(&mapBig, sizeof(int), 1, mapF);
		if (mapBig != charVec.numPts) {
			cout << "invalid map file target size; expected " << charVec.numPts << "; found " << mapBig << endl;
			fclose(mapF);
			return;
		}
		if (mapSmall != loadPts) {
			cout << "invalid map file; expected " << loadPts << "; found " << mapSmall << endl;
			fclose(mapF);
			return;
		}

		// load original data
		CharVecPt *tempCV = new CharVecPt[loadPts];
		for (i=0; i < loadPts; i++) {
			fread(&(tempCV[i].numInfluences), sizeof(int), 1, f);
			tempCV[i].influences = new int[tempCV[i].numInfluences];
			fread(tempCV[i].influences, sizeof(int), tempCV[i].numInfluences, f);
			
			fread(&n), sizeof(int), 1, f);
			tempCV[i].data.resize(n);
			fread(tempCV[i].data.n, sizeof(double), n, f);
		}

		// apply mapping
		cout << "applying mapping..." << endl;
		for (i=0; i < mapBig; i++) {
			int numSources;
			fread(&numSources, sizeof(int), 1, mapF);

			charVec.cvPts[i].data.zeroElements();

			for (j=0; j < numSources; j++) {
				int curV;
				float curW;
				fread(&curV, sizeof(int), 1, mapF);
				fread(&curW, sizeof(float), 1, mapF);

				for (k=0; k < charVec.cvPts[i].numInfluences; k++) {
					int sInf = tempCV[curV].getInfIndex(charVec.cvPts[i].influences[k]);
					if (sInf > -1) {
						for (x=0; x < charVec.cvPts[i].data.size(); x++)
// FIX ME!
							charVec.cvPts[i].data[k * charVec.numComponents * 3 * CV_DEGREE + x] +=
								curW * tempCV[curV].data[sInf * charVec.numComponents * 3 * CV_DEGREE + x];
					}
				}
			}
		}

		// free original data
		delete []tempCV;

		fclose(mapF);
*/	}
}

void cvLoadWeighted(const char *params) {
	cout << "cvLoadWeighted is under construction." << endl;
	/*
	char fname[80];
	fname[0] = 0;
	double weight = 1.0;
	bool zero = false;
	params = extractString(params, fname, 80);
	params = extractDouble(params, &weight);
	params = extractBool(params, &zero);

	FILE *f;
	if (!openFile(&f, fname, "rb", "charvec file")) {
		return;
	}

	int i, j, n, loadPts;

	fread(&loadPts, sizeof(int), 1, f);
	fread(&i, sizeof(int), 1, f); // numComponents
//	i = CV_DEGREE;
//	fread(&i, sizeof(int), 1, f); // degree

	if (loadPts == charVec.numPts) {
		for (i=0; i < charVec.numPts; i++) {
			fread(&n, sizeof(int), 1, f);
			if (n != charVec.cvPts[i].numInfluences) {
				cout << "# of influences mismatch on point " << i << "; aborting!" << endl;
				return;
			}
			for (j=0; j < charVec.numComponents*CV_DEGREE*3*n; j++) {
				double d;
				fread(&d, sizeof(double), 1, f);
				d *= weight;
				if (zero)
					charVec.cvPts[i].data[j] = d;
				else
					charVec.cvPts[i].data[j] += d;
			}
		}
	}
	else {
		cout << "point count mismatch; expected " << charVec.numPts << "; found " << i << endl;
		return;
	}

	fclose(f);
	*/
}

void cvDumpPt(const char *params) {
	int pt = 0;
	params = extractInt(params, &pt);

	int i, j;

	ofstream out("dump.txt");
	out << "examples" << endl;
	for (i=0; i < cvGF->numExamples; i++)
		out << i << " " << cvGF->examples[i].points[pt] << endl;

	out << endl << "dofs" << endl;
	for (i=0; i < cvGF->numExamples; i++) {
		out << i << " ";
		for (j=0; j < cvGF->examples[i].numDofs; j++) {
			out << cvGF->examples[i].dofs[j] << " ";
		}
		out << endl;
	}

	out << endl << "fit" << endl;
	out << charVec.cvPts[pt].numInfluences << endl;
	for (i=0; i < charVec.cvPts[pt].numInfluences; i++)
		out << charVec.cvPts[pt].influences[i] << " ";
	out << endl;
	for (i=0; i < charVec.cvPts[pt].data.size(); i++) {
		out << charVec.cvPts[pt].data[i] << " ";
	}
	out << endl;
	out.close();
}

void cvDumpDOFs(const char *params) {
	ofstream out("dofs.txt");

	int i, j;
	for (i=0; i < cvGF->numExamples; i++) {
		out << cvGF->examples[i].fname;
		for (j=0; j < cvGF->examples[i].numDofs; j++) {
			out << ", " << cvGF->examples[i].dofs[j];
		}
		out << endl;
	}
	out.close();
}

void cvShowInf(const char *params) {
	int inf = 0;
	params = extractInt(params, &inf);

	int i;
	for (i=0; i < charVec.numPts; i++) {
		int ind = charVec.cvPts[i].getInfIndex(inf);

		if (ind > 0)
			charVec.tm->getPtColor(i) = Vec3d(1.0, 1.0-charVec.cvPts[i].infWeights[ind], 1.0-charVec.cvPts[i].infWeights[ind]);
		else
			charVec.tm->getPtColor(i) = Vec3d(1, 1, 1);
	}
}

void cvSaveSet(const char *params) {
	char fname[80];
	fname[0] = 0;
	params = extractString(params, fname, 80);

	FILE *f;
	if (!openFile(&f, fname, "wb", "charvec set"))
		return;

	// version:
	int i = 0;
	fwrite(&i, sizeof(int), 1, f);
	// # of points:
	fwrite(&cvGF->cv->numPts, sizeof(int), 1, f);
	// # of components:
	fwrite(&cvGF->cv->numComponents, sizeof(int), 1, f);
	// # of P values:
//	i = NUM_P_VALUES;
//	fwrite(&i, sizeof(int), 1, f);

	// data:
//	fwrite(cvGF->curVars.n, sizeof(double), cvGF->curVars.size(), f);

	fclose(f);
}

void cvLoadSet(const char *params) {
	char fname[80];
	fname[0] = 0;
	params = extractString(params, fname, 80);

	// ...
}

void cvMakeMovie(const char *params) {
	double lo, hi;
	int numFrames;
	params = extractInt(params, &numFrames);
	params = extractDouble(params, &lo);
	params = extractDouble(params, &hi);
	char fName[80];

	int i;
	for (i=0; i < numFrames; i++) {
		double comps[1];
		comps[0] = lo + (hi-lo) * (1.0 * i / (numFrames-1));
//		cvGF->constructCV(comps, 1);
		charVec.updateTM();

		redrawVNow();
		sprintf(fName, "anim\\sc%04d.tga", i);
		uiScreenshot(fName);
	}
}

void cvSavePly(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);

	charVec.tm->savePly(fname);
}

void cvLoadTexture(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);
	char mapping[80];
	mapping[0] = 0;
	params = extractString(params, mapping, 80);

	int expectedSize = charVec.numPts * 3;

	FILE *mapF = NULL;
	int mapSmall, mapBig;
	if (mapping[0]) {
		if (!openFile(&mapF, mapping, "rb", "mapping"))
			return;
		fread(&mapSmall, sizeof(int), 1, mapF);
		fread(&mapBig, sizeof(int), 1, mapF);
		expectedSize = mapSmall * 3;

		if (mapBig != charVec.numPts) {
			cout << "invalid map file; expected " << charVec.numPts << "; found " << mapBig << endl;
			fclose(mapF);
			return;
		}
	}

	FILE *f;
	if (!openFile(&f, fname, "rb", "points"))
		return;

	char version;
	fread(&version, sizeof(char), 1, f);

	if (version == '0') {
		int i, j;
		fread(&i, sizeof(int), 1, f);

		if (i != expectedSize) {
			cout << "point size mismatch; expected " << expectedSize << "; read " << i << endl;
			// temp hack:
			for (i=0; i < cvGF->cv->numPts; i++) {
				cvGF->cv->tm->getPtColor(i) = Vec3d(1,1,1);
			}
		}
		else {
			if (mapping[0]) {
				Vec3d *fileColors = new Vec3d[mapSmall];
				fread(fileColors, sizeof(Vec3d), mapSmall, f);

				for (i=0; i < mapBig; i++) {
					int numSources;
					fread(&numSources, sizeof(int), 1, mapF);

					Vec3d &v = cvGF->cv->tm->getPtColor(i);
					v = Vec3d();

					for (j=0; j < numSources; j++) {
						int curV;
						float curW;
						fread(&curV, sizeof(int), 1, mapF);
						fread(&curW, sizeof(float), 1, mapF);
						v += curW * fileColors[curV];
					}
				}
				fclose(mapF);
			}
			else {
				for (i=0; i < cvGF->cv->numPts; i++) {
					Vec3d *v = &(cvGF->cv->tm->getPtColor(i));
					fread(v, sizeof(double), 3, f);
				}
			}
		}
		fclose(f);
	}
	else {
		cout << "warning in cvLoadTexture: unsupported format" << endl;
	}
}

void cvSetLocalFromTM(const char *params) {
	int ex = 0;
	params = extractInt(params, &ex);
	charVec.setFromLocalTM(cvGF->examples[ex].iTrans);
}

void cvSaveLocal(const char *params) {
	char fname[255];
	params = extractString(params, fname, 255);

	FILE *f;
	if (!openFile(&f, fname, "wb", "local coordinate file"))
		return;

	fwrite(&charVec.numPts, sizeof(int), 1, f);
	int i;
	for (i=0; i < charVec.numPts; i++) {
		fwrite(&(charVec.cvPts[i].localPos), sizeof(Vec3d), 1, f);
	}

	fclose(f);
}

void cvLoadLocal(const char *params) {
	char fname[255];
	params = extractString(params, fname, 255);

	FILE *f;
	if (!openFile(&f, fname, "rb", "local coordinate file"))
		return;

	int i;
	fread(&i, sizeof(int), 1, f);
	if (i != charVec.numPts) {
		cout << "warning in cvLoadLocal: expected " << charVec.numPts << " points; found " << i << endl;
		fclose(f);
		return;
	}
	for (i=0; i < charVec.numPts; i++) {
		fread(&(charVec.cvPts[i].localPos), sizeof(Vec3d), 1, f);
	}

	fclose(f);
}

void cvSaveExample(const char *params) {
	int ex;
	params = extractInt(params, &ex);

	CVExample newEx = cvGF->examples[ex];
	newEx.numPoints = charVec.numPts;
	newEx.points = new Vec3d[newEx.numPoints];
	newEx.conf = new double[newEx.numPoints];
	int i;
	for (i=0; i < newEx.numPoints; i++) {
		newEx.points[i] = charVec.tm->getPt(i);
		newEx.conf[i] = charVec.tm->getPtColor(i)[1];
	}
	newEx.save(cvGF->examples[ex].fname);
	delete []newEx.points;
	delete []newEx.conf;
}

void cvSetCVFromExample(const char *params) {
	int ex = 0;
	params = extractInt(params, &ex);


	charVec.updateCurComponents(cvGF->examples[ex].dofs);
	charVec.updateLocalPos(NULL);

	int i;
	for (i=0; i < charVec.numPts; i++) {
		Vec3d delta = cvGF->examples[ex].points[i] - charVec.cvPts[i].localPos;
		charVec.cvPts[i].data[0] += delta[0];
		charVec.cvPts[i].data[1] += delta[1];
		charVec.cvPts[i].data[2] += delta[2];
	}

	charVec.updateCurComponents(cvGF->examples[ex].dofs);
	charVec.updateLocalPos(NULL);
	charVec.updateTM();
	redrawV();
}

void cvSaveExMatrix(const char *params) {
	char fname[256];
	params = extractString(params, fname, 256);

	ofstream f;
	if (!openOFStream(&f, fname, "matrix file"))
		return;

	int i, j, k;
	for (i=0; i < charVec.numPts; i++) {
		for (j=0; j < 3; j++) {
			for (k=0; k < cvGF->numExamples; k++) {
				cout << cvGF->examples[k].points[i][j] << " ";
			}
			cout << endl;
		}
	}

	f.close();
}


void cvLoadPose(const char *params) {
	char fname[256];
	params = extractString(params, fname, 256);

	ifstream in;
	if (!openIFStream(&in, fname, "pose file"))
		return;
	scSkin.skel->loadPose(in);
	in.close();

	scSkin.skel->updateCoords();
	scSkin.updatePoints();
}

void cvInitSurfaceMatch(const char *params) {
	int i;
	int example = 0;
	char markerName[255], mRefsName[255];

	params = extractInt(params, &example);
	params = extractString(params, markerName, 255);
	params = extractString(params, mRefsName, 255);

	// set global coordinates (0th example)
//	cvGF->cv->updateTM(cvGF->examples[0].trans);

	// load lock shape
//	cvLoadTexture("data/lock-hi.disp");
//	int i;
//	for (i=0; i < cvGF->cv->numPts; i++) {
//		cvGF->cv->tm->getPtColor(i) = Vec3d(1, 1, 1);
//	}

	// initialize goal function
	if (!matchGF) {
		matchGF = new LADeformationGoalFunction(cvGF->cv->tm, cvGF->cv->neighbors);	

		// load marker refs
		ifstream is;
		int i, n;
		if (!openIFStream(&is, mRefsName, "marker refs"))
			return;
		is >> n;
		matchGF->markerRefs.resize(n);
		for (i=0; i < n; i++) {
			is >> matchGF->markerRefs[i];
		}
		is.close();
		cout << "loaded " << n << " marker refs" << endl;
	}
	else {
		int i;
		for (i=0; i < cvGF->cv->tm->numPts(); i++) {
			matchGF->origVerts[i] = cvGF->cv->tm->getPt(i);
		}
	}

	// load markers
	matchGF->markers = new MarkerSet;
	if (markerName[strlen(markerName)-3] == 'm')
		matchGF->markers->loadFromMkr(markerName);
	else {
		cout << "loading from landmarks..." << endl;
		matchGF->markers->loadFromLandmarks(markerName);
		// wipe out extraneous markers
		for (i=73; i < 85; i++) {
			matchGF->markers->markers[i].pos = Vec3d();
		}
	}
	
	// ignore patella
//	matchGF->markers->markers[73].pos = Vec3d();
//	matchGF->markers->markers[74].pos = Vec3d();

	matchMeshes[example].calcNormals();
	matchMeshes[example].calcHBB(16);
	matchGF->prepareTriMesh(&(matchMeshes[example]));
	matchGF->zeroDeformation();
/*
	// apply skinning
	if (matchGF->varsPerVert == 10) {
		int j, index = 0;
		for (i=0; i < cvGF->cv->numPts; i++) {
			Mat4d m = cvGF->examples[0].trans[i];
			m = m.inverse();
			m = cvGF->examples[example].trans[i] * m;

			QuatNorm q = matToQuat(Mat3d(m.n[0],m.n[1],m.n[2],m.n[4],m.n[5],m.n[6],m.n[8],m.n[9],m.n[10]));
			matchGF->vars[i*10+0] = q.x;
			matchGF->vars[i*10+1] = q.y;
			matchGF->vars[i*10+2] = q.z;
			matchGF->vars[i*10+3] = q.w;
			matchGF->vars[i*10+4] = 1;
			matchGF->vars[i*10+5] = 1;
			matchGF->vars[i*10+6] = 1;
			matchGF->vars[i*10+7] = m.n[3];
			matchGF->vars[i*10+8] = m.n[7];
			matchGF->vars[i*10+9] = m.n[11];
		}
		matchGF->applyDef(matchGF->vars);
	}
	else if (matchGF->varsPerVert == 12) {
		int j, index = 0;
		for (i=0; i < cvGF->cv->numPts; i++) {
			Mat4d m = cvGF->examples[0].trans[i];
			m = m.inverse();
			m = cvGF->examples[example].trans[i] * m;
			for (j=0; j < 12; j++)
				matchGF->vars[index++] = m.n[j];
		}
		matchGF->applyDef(matchGF->vars);
	}
	else 
		cout << "wrong number of components for matchGF; deformation not applied" << endl;
		*/
}

void cvRunSurfaceMatch(const char *params) {
	if (!matchGF) {
		cout << "warning: goal function not initialized in cvRunSurfaceMatch" << endl;
		return;
	}
	if (solver) {
		cout << "warning: solver already running in cvRunSurfaceMatch" << endl;
		return;
	}

	params = extractDouble(params, &(matchGF->surfaceMatchWeight));
	params = extractDouble(params, &(matchGF->smoothnessWeight));
	params = extractDouble(params, &(matchGF->markerMatchWeight));
	int maxIter = 1000;
	params = extractInt(params, &maxIter);

	solver = new LBFGSSolver(matchGF);
	solver->solve(1e+3, 1e-5, matchGF->vars, maxIter);

	delete solver;
	solver = NULL;
}

void cvStop(const char *params) {
	if (solver)
		solver->stopNow = true;
}

void cvDumpVariance(const char *params) {
	int i, j, k, n;
	double sumX, sumX2;

	for (i=1; i < charVec.numComponents; i++) {
		sumX = 0;
		sumX2 = 0;
		n = 0;

		for (j=0; j < charVec.numPts; j++) {
			for (k=i*3; k < charVec.cvPts[j].data.size(); k += charVec.numComponents*3) {
				sumX += charVec.cvPts[j].data[k+0];
				sumX2 += sqr(charVec.cvPts[j].data[k+0]);
				sumX += charVec.cvPts[j].data[k+1];
				sumX2 += sqr(charVec.cvPts[j].data[k+1]);
				sumX += charVec.cvPts[j].data[k+2];
				sumX2 += sqr(charVec.cvPts[j].data[k+2]);
				n += 3;
			}
		}
		cout << "component " << i << ": " << (1.0/(n-1.0))*(sumX2 - (2.0/n)*sqr(sumX) + n*sqr(sumX/n)) << endl;
	}
}

void cvAnimPoses(const char *params) {
	vector<int> poses;
	vector<int> frames;
	int i, j;
	char fname[256];
	bool updateLocal = true;

	params = extractBool(params, &updateLocal);
	
	while (1) {
		i = -1;
		params = extractInt(params, &i);
		if (i < 0)
			break;
		poses.push_back(i);
		i = -1;
		params = extractInt(params, &i);
		if (i < 0)
			break;
		frames.push_back(i);
	}

	// wipe color
	for (i=0; i < charVec.tm->numPts(); i++)
		charVec.tm->getPtColor(i) = Vec3d(0.8, 0.8, 0.8);

	int frame = 0;
	int numDofs = cvGF->examples[0].numDofs;
	double *dofs = new double[numDofs];
	for (i=0; i < frames.size(); i++) {
		for (j=0; j < frames[i]; j++) {
			double interp = 1.0 * j / (frames[i]);
			
			charVec.skel->interpVals(&matchPoses[poses[i]], &matchPoses[poses[i+1]], interp);
			charVec.skel->updateCoords();
			scSkin.skel->copyVals(charVec.skel);
			scSkin.skel->updateCoords();
			scSkin.updatePoints();
			charVec.dofsFromSkel(charVec.skel, dofs);
			charVec.updateCurComponents(dofs);
			if (updateLocal)
				cvGF->cv->updateLocalPos(NULL);
			cvGF->cv->updateTM();
			redrawVNow();

			sprintf(fname, "anim/frame%04d.tga", frame);
			uiScreenshot(fname);
			uiWait();
			frame++;
		}
	}
	delete []dofs;
}

void cvCreateMapping(const char *params) {
	char fnameS[80], fnameB[80], fnameM[80]; //, fnameMI[80];
	params = extractString(params, fnameS, 80);
	params = extractString(params, fnameB, 80);
	params = extractString(params, fnameM, 80);

	createMapping(fnameS, fnameB, fnameM);
}

void cvSaveDisp(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);

	FILE *f;
	if (!openFile(&f, fname, "wb", "displacements")) {
		return;
	}

	char version = '0';
	fwrite(&version, sizeof(char), 1, f);
	if (matchGF) {
		int size = matchGF->vars.size();
		fwrite(&size, sizeof(int), 1, f);
		fwrite(matchGF->vars.n, sizeof(double), size, f);
	}
	else if (edgeMatchGF) {
		int size = edgeMatchGF->vars.size();
		fwrite(&size, sizeof(int), 1, f);
		fwrite(edgeMatchGF->vars.n, sizeof(double), size, f);
	}
	fclose(f);
}

void cvLoadDisp(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);
	char mapping[80];
	mapping[0] = 0;
	params = extractString(params, mapping, 80);

	int i, j; 
	int ofs = 0;
	int expectedSize = matchGF->vars.size();

	FILE *mapF = NULL;
	int mapSmall, mapBig;
	if (mapping) {
		if (!openFile(&mapF, mapping, "rb", "mapping"))
			return;
		fread(&mapSmall, sizeof(int), 1, mapF);
		fread(&mapBig, sizeof(int), 1, mapF);
		expectedSize = mapSmall * matchGF->varsPerVert;

		if (mapBig != charVec.numPts) {
			cout << "invalid map file; expected " << charVec.numPts << "; found " << mapBig << endl;
			fclose(mapF);
			return;
		}
	}

	FILE *f;
	if (!openFile(&f, fname, "rb", "displacements"))
		return;

	char version;
	fread(&version, sizeof(char), 1, f);

	if (version != '0') {
		cout << "unsupported version" << endl;
		return;
	}

	int fileSize;
	fread(&fileSize, sizeof(int), 1, f);

	double *fileDisps = new double[fileSize];
	baAssert(fileDisps != NULL, "can't allocate memory for displacements");
	fread(fileDisps, sizeof(double), fileSize, f);

	if (fileSize != expectedSize) {
		cout << "file has " << fileSize << " points; expected " << expectedSize << endl;
		fclose(f);
		return;
	}

	if (mapF) {
		cout << "applying mapping from " << mapSmall << " to " << mapBig << endl;
		for (i=0; i < mapBig; i++) {
			int numSources;
			fread(&numSources, sizeof(int), 1, mapF);

			int curOfs = ofs + i * matchGF->varsPerVert;
			int k;
			for (k=0; k < matchGF->varsPerVert; k++)
				matchGF->vars[curOfs + k] = 0;

			for (j=0; j < numSources; j++) {
				int curV;
				float curW;
				fread(&curV, sizeof(int), 1, mapF);
				fread(&curW, sizeof(float), 1, mapF);
				for (k=0; k < matchGF->varsPerVert; k++)
					matchGF->vars[curOfs + k] += curW * fileDisps[curV * matchGF->varsPerVert + k];
			}
		}
	}
	else {
		for (i=0; i < expectedSize; i++) {
			matchGF->vars[ofs + i] = fileDisps[i];
		}
	}

	delete []fileDisps;
	fclose(f);

	if (mapF)
		fclose(mapF);

	matchGF->applyDef(matchGF->vars);
}

void cvSaveMirrorMap(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);

	FILE *f;
	if (!openFile(&f, fname, "wb", "mirror map"))
		return;

	fwrite(&charVec.numPts, sizeof(int), 1, f);
	fwrite(charVec.mirrorMap, sizeof(int), charVec.numPts, f);
	fclose(f);
}

void cvLoadMocap(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);
	mcData.load(fname);
}

void cvMocapFrame(const char *params) {
	int frame = 0;
	params = extractInt(params, &frame);
	mcData.toSkel(frame, scSkin.skel);

	scSkin.skel->updateCoords();
	scSkin.updatePoints();
	double dofs[100];
	CharVec::dofsFromSkel(scSkin.skel, dofs);
	cvGF->cv->updateCurComponents(dofs);
	cvGF->cv->updateLocalPos(NULL);
	cvGF->cv->updateTM();
	redrawV();
}

void cvMocapAnim(const char *params) {
	int minFrame = 0;
	int maxFrame = 100;
	int frameStep = 4;
	bool save = false;
	params = extractInt(params, &minFrame);
	params = extractInt(params, &maxFrame);
	params = extractInt(params, &frameStep);
	params = extractBool(params, &save);

	int i, frame = 0;
	char s[256];
	for (i = minFrame; i < maxFrame; i += frameStep) {
		sprintf(s, "%d", i);
		cvMocapFrame(s);
		redrawVNow();

		if (save) {
			sprintf(s, "anim/mc%04d.tga", frame);
			uiScreenshot(s);
			uiWait();
			frame++;
		}
	}
}

void cvSaveAni(const char *params) {
	char fname[256];
	params = extractString(params, fname, 256);
	ofstream out;
	if (!openOFStream(&out, fname, "ani file"))
		return;

	int minFrame = 0;
	int maxFrame = mcData.frames;
	int frameStep = 1;

	const int NUM_DOFS = 75;
	//const int NUM_DOFS = 15;
	const char dofNames[NUM_DOFS][20] = {
		"position_x", "position_y", "position_z",
		"base_x", "base_y", "base_z", "base_w",
		"waistQ_x", "waistQ_y", "waistQ_z", "waistQ_w",
		"abdomen_x", "abdomen_y", "abdomen_z", "abdomen_w",
		"lClavicle_x", "lClavicle_y", "lClavicle_z", "lClavicle_w",
		"lShoulder_x", "lShoulder_y", "lShoulder_z", "lShoulder_w",
		"lElbow_x", "lElbow_y", "lElbow_z", "lElbow_w",
		"lWrist_x", "lWrist_y", "lWrist_z", "lWrist_w",
		"rClavicle_x", "rClavicle_y", "rClavicle_z", "rClavicle_w",
		"rShoulder_x", "rShoulder_y", "rShoulder_z", "rShoulder_w",
		"rElbow_x", "rElbow_y", "rElbow_z", "rElbow_w",
		"rWrist_x", "rWrist_y", "rWrist_z", "rWrist_w",
		"neck_x", "neck_y", "neck_z", "neck_w",
		"lHip_x", "lHip_y", "lHip_z", "lHip_w",
		"lKnee_x", "lKnee_y", "lKnee_z", "lKnee_w",
		"lAnkle_x", "lAnkle_y", "lAnkle_z", "lAnkle_w",
		"rHip_x", "rHip_y", "rHip_z", "lHip_w",
		"rKnee_x", "rKnee_y", "rKnee_z", "rKnee_w",
		"rAnkle_x", "rAnkle_y", "rAnkle_z", "rAnkle_w"};

	out << maxFrame << " " << NUM_DOFS << endl << endl;

	double xFactor = -1;
	double yFactor = -1;
	double zFactor = -1;
	double wMult = -1;
	int i, j, ind;
	int index = 0;
	QuatNorm q;
	double *redofs = new double[maxFrame * NUM_DOFS];
	for (i = minFrame; i < maxFrame; i += frameStep) {
		int ofs = i * mcData.frameSize;

		// base
		redofs[index++] = 0;//mcData.data[ofs + mcData.mappings[0] + 0];
		redofs[index++] = 0;//mcData.data[ofs + mcData.mappings[0] + 1];
		redofs[index++] = 0;//mcData.data[ofs + mcData.mappings[0] + 2];
		ind = 0;
		q = QuatNorm(
			xFactor * mcData.data[ofs+mcData.mappings[ind]+0]*DEG_TO_RAD, 
			yFactor * mcData.data[ofs+mcData.mappings[ind]+1]*DEG_TO_RAD, 
			zFactor * mcData.data[ofs+mcData.mappings[ind]+2]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;

		// waist and abdomen
		ind = 2;
		q = QuatNorm(
			xFactor * mcData.data[ofs+mcData.mappings[ind]+0]*DEG_TO_RAD, 
			yFactor * mcData.data[ofs+mcData.mappings[ind]+1]*DEG_TO_RAD, 
			zFactor * mcData.data[ofs+mcData.mappings[ind]+2]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;
		ind = 3;
		q = QuatNorm(
			xFactor * mcData.data[ofs+mcData.mappings[ind]+0]*DEG_TO_RAD, 
			yFactor * mcData.data[ofs+mcData.mappings[ind]+1]*DEG_TO_RAD, 
			zFactor * mcData.data[ofs+mcData.mappings[ind]+2]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;

		// left arm
		ind = 4;
		q = QuatNorm(
			xFactor * mcData.data[ofs+mcData.mappings[ind]+0]*DEG_TO_RAD, 
			yFactor * mcData.data[ofs+mcData.mappings[ind]+1]*DEG_TO_RAD, 
			zFactor * mcData.data[ofs+mcData.mappings[ind]+2]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;
		ind = 5;
		q = QuatNorm(
			xFactor * mcData.data[ofs+mcData.mappings[ind]+0]*DEG_TO_RAD, 
			yFactor * mcData.data[ofs+mcData.mappings[ind]+1]*DEG_TO_RAD, 
			zFactor * mcData.data[ofs+mcData.mappings[ind]+2]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;
		q = QuatNorm(0, mcData.data[ofs + mcData.mappings[6] + 0]*DEG_TO_RAD, mcData.data[ofs + mcData.mappings[7] + 0]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;
		q = QuatNorm(0, mcData.data[ofs + mcData.mappings[8] + 0]*DEG_TO_RAD, 0);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;

		// right arm
		ind = 9;
		q = QuatNorm(
			xFactor * mcData.data[ofs+mcData.mappings[ind]+0]*DEG_TO_RAD, 
			yFactor * mcData.data[ofs+mcData.mappings[ind]+1]*DEG_TO_RAD, 
			zFactor * mcData.data[ofs+mcData.mappings[ind]+2]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;	
		ind = 10;
		q = QuatNorm(
			xFactor * mcData.data[ofs+mcData.mappings[ind]+0]*DEG_TO_RAD, 
			yFactor * mcData.data[ofs+mcData.mappings[ind]+1]*DEG_TO_RAD, 
			zFactor * mcData.data[ofs+mcData.mappings[ind]+2]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;
		q = QuatNorm(0, mcData.data[ofs + mcData.mappings[11] + 0]*DEG_TO_RAD, mcData.data[ofs + mcData.mappings[12] + 0]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;
		q = QuatNorm(0, mcData.data[ofs + mcData.mappings[13] + 0]*DEG_TO_RAD, 0);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;

		// head
		ind = 24;
		q = QuatNorm(
			xFactor * mcData.data[ofs+mcData.mappings[ind]+0]*DEG_TO_RAD, 
			yFactor * mcData.data[ofs+mcData.mappings[ind]+1]*DEG_TO_RAD, 
			zFactor * mcData.data[ofs+mcData.mappings[ind]+2]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;

		// left leg
		ind = 14;
		q = QuatNorm(
			xFactor * mcData.data[ofs+mcData.mappings[ind]+0]*DEG_TO_RAD, 
			yFactor * mcData.data[ofs+mcData.mappings[ind]+1]*DEG_TO_RAD, 
			zFactor * mcData.data[ofs+mcData.mappings[ind]+2]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;
		q = QuatNorm(0, mcData.data[ofs + mcData.mappings[15] + 0]*DEG_TO_RAD, mcData.data[ofs + mcData.mappings[16] + 0]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;
		q = QuatNorm(0, mcData.data[ofs + mcData.mappings[17] + 0]*DEG_TO_RAD, 0);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;

		// right leg
		ind = 19;
		q = QuatNorm(
			xFactor * mcData.data[ofs+mcData.mappings[ind]+0]*DEG_TO_RAD, 
			yFactor * mcData.data[ofs+mcData.mappings[ind]+1]*DEG_TO_RAD, 
			zFactor * mcData.data[ofs+mcData.mappings[ind]+2]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;
		q = QuatNorm(0, mcData.data[ofs + mcData.mappings[20] + 0]*DEG_TO_RAD, mcData.data[ofs + mcData.mappings[21] + 0]*DEG_TO_RAD);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;
		q = QuatNorm(0, mcData.data[ofs + mcData.mappings[22] + 0]*DEG_TO_RAD, 0);
		q.normalize();
		redofs[index++] = q.x;
		redofs[index++] = q.y;
		redofs[index++] = q.z;
		redofs[index++] = wMult*q.w;

		if (i == minFrame && index != NUM_DOFS) {
			cout << "uh-oh... index is " << index << "; should be " << NUM_DOFS << endl;
		}
	}

	for (i=0; i < NUM_DOFS; i++) {
		out << dofNames[i] << " ";
		for (j = minFrame; j < maxFrame; j += frameStep) {
			out << redofs[j*NUM_DOFS + i] << " ";
		}
		out << endl;
	}
}

void cvLoadMu(const char *params) {
	char fname[256];
	params = extractString(params, fname, 256);
	ifstream in(fname);
	int i, j;
	for (i=0; i < cvGF->numExamples; i++) {
		for (j=0; j < charVec.numComponents; j++) {
			in >> cvGF->curMu[i][j];
		}
	}
	in.close();
}

void cvAutoskinInit(const char *params) {
	if (!cvSGF)
		cvSGF = new CVSkinningGF();
	cvSGF->init(cvGF, &scSkin, matchPoses);
}

void cvAutoskinSolve(const char *params) {
	solver = new LBFGSSolver(cvSGF);
//	runTest(cvSGF, cvSGF->vars);
//	exit(0);
	solver->solve(1e-7, 1e-5, cvSGF->vars, 200);
	cvSGF->varsToSkin(cvSGF->vars);
}

void cvSkinMatchInit(const char *params) {
/*	int i;
	char markerName[255], mRefsName[255];

	params = extractString(params, markerName, 255);
	params = extractString(params, mRefsName, 255);

	// load lock shape
//	cvLoadTexture("data/lock-hi.disp");
//	int i;
//	for (i=0; i < cvGF->cv->numPts; i++) {
//		cvGF->cv->tm->getPtColor(i) = Vec3d(1, 1, 1);
//	}

	// initialize goal function
	if (!skinMatchGF) {
		skinMatchGF = new SkinMatchGF(cvGF->numExamples, cvGF->cv->tm, &scSkin, matchPoses);	

		// load marker refs
		ifstream is;
		int n;
		if (!openIFStream(&is, mRefsName, "marker refs"))
			return;
		is >> n;
		skinMatchGF->markerRefs.resize(n);
		for (i=0; i < n; i++) {
			is >> skinMatchGF->markerRefs[i];
		}
		is.close();
		cout << "loaded " << n << " marker refs" << endl;
	}
	else {
		int i;
		for (i=0; i < cvGF->cv->tm->numPts(); i++) {
			skinMatchGF->origVerts[i] = cvGF->cv->tm->getPt(i);
		}
	}

	// load markers
	skinMatchGF->markers = new MarkerSet;
	if (markerName[strlen(markerName)-3] == 'm')
		skinMatchGF->markers->loadFromMkr(markerName);
	else {
		cout << "loading from landmarks..." << endl;
		skinMatchGF->markers->loadFromLandmarks(markerName);
		// wipe out extraneous markers
		int i;
		for (i=73; i < 85; i++) {
			skinMatchGF->markers->markers[i].pos = Vec3d();
		}
	}
	
	// ignore patella
//	matchGF->markers->markers[73].pos = Vec3d();
//	matchGF->markers->markers[74].pos = Vec3d();

	cout << "initializing meshes..." << endl;
	for (i=0; i < cvGF->numExamples; i++) {
		matchMeshes[i].calcNormals();
		matchMeshes[i].calcHBB(16);
	}
	cout << "initializing goal function..." << endl;
	skinMatchGF->prepareTriMesh(matchMeshes);
	skinMatchGF->zeroDeformation();

	skinMatchGF->edgeMatchTM = new TriMesh();
	skinMatchGF->edgeMatchTM->loadPly("data/2a.ply");*/
}

void cvSkinMatchSolve(const char *params) {
/*	params = extractDouble(params, &(skinMatchGF->surfaceMatchWeight));
	params = extractDouble(params, &(skinMatchGF->smoothnessWeight));
	params = extractDouble(params, &(skinMatchGF->markerMatchWeight));
	int maxIter = 1000;
	params = extractInt(params, &maxIter);

	solver = new LBFGSSolver(skinMatchGF);
//	runTest(cvSGF, cvSGF->vars);
//	exit(0);
	solver->solve(1e-7, 1e-5, skinMatchGF->vars, maxIter);
//	skinMatchGF->varsToSkin(skinMatchGF->vars);*/
}

void cvEdgeMatchInit(const char *params) {
	int i;
	int example = 0;
	char markerName[255], mRefsName[255];

	params = extractInt(params, &example);
	params = extractString(params, markerName, 255);
	params = extractString(params, mRefsName, 255);

	TriMesh *m0 = new TriMesh(), *m1 = new TriMesh();
//	m0->loadPly("data/simple0.ply");
//	m0->loadPly("data/half-james-5k-m.ply");
//	m0->loadPly("data/pig1850.ply");
//	m0->loadPly("data/cube20.ply");
//	m0->loadPly("data/head.ply");
	m0->loadPly("y:ba-small.ply");
//	m0->loadPly("data/strip2.ply");
//	m1->loadPly("data/simple0.ply");
	m1->loadPly("data/strip.ply");
	m0->calcNormals();
//	cvGF->cv->tm = m0;

	// set global coordinates
//	cvGF->cv->updateTM(cvGF->examples[example].trans);

	// initialize goal function
	if (!edgeMatchGF) {
//		TriMesh *matchTM = new TriMesh();
//		matchTM->loadPly("data/2a.ply");
		edgeMatchGF = new EdgeMatchGF(cvGF->cv->tm, cvGF->cv->neighbors, cvGF->cv->tm);

		// load marker refs
		ifstream is;
		int i, n;
		if (!openIFStream(&is, mRefsName, "marker refs"))
			return;
		is >> n;
		edgeMatchGF->markerRefs.resize(n);
		for (i=0; i < n; i++) {
			is >> edgeMatchGF->markerRefs[i];
		}
		is.close();
		cout << "loaded " << n << " marker refs" << endl;

//		for (i=0; i < cvGF->cv->tm->numPts(); i++) {
//			edgeMatchGF->origVerts[i] = matchTM->getPt(i);
//		}
	}
	else {
		int i;
		for (i=0; i < cvGF->cv->tm->numPts(); i++) {
			edgeMatchGF->origVerts[i] = cvGF->cv->tm->getPt(i);
		}
	}

	// load markers
	edgeMatchGF->markers = new MarkerSet;
	if (markerName[strlen(markerName)-3] == 'm')
		edgeMatchGF->markers->loadFromMkr(markerName);
	else {
		cout << "loading from landmarks..." << endl;
		edgeMatchGF->markers->loadFromLandmarks(markerName);
		// wipe out extraneous markers
		int i;
		for (i=73; i < 85; i++) {
			edgeMatchGF->markers->markers[i].pos = Vec3d();
		}
	}
	
	// ignore patella
//	matchGF->markers->markers[73].pos = Vec3d();
//	matchGF->markers->markers[74].pos = Vec3d();

	matchMeshes[example].calcNormals();
	matchMeshes[example].calcHBB(16);
	edgeMatchGF->prepareTriMesh(&(matchMeshes[example]));
	edgeMatchGF->zeroDeformation();

	/*
	// color based on # of edges
	int *counts = new int[edgeMatchGF->curMesh->numPts()];
	memset(counts, 0, sizeof(int)*edgeMatchGF->curMesh->numPts());
	for (i=0; i < edgeMatchGF->numEdges; i++) {
		counts[edgeMatchGF->templateEdges[i].v0]++;
		counts[edgeMatchGF->templateEdges[i].v1]++;
	}
	for (i=0; i < edgeMatchGF->curMesh->numPts(); i++) {
		if (counts[i] < 2)
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d();
		else if (counts[i] == 3)
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d(1,0,0);
		else if (counts[i] == 4)
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d(0,1,0);
		else if (counts[i] == 5)
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d(0,0,1);
		else if (counts[i] == 6)
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d(1,1,0);
		else if (counts[i] == 7)
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d(0,1,1);
		else 
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d(1,1,1);
	}*/

	// initialization
/*	cvGF->cv->updateTM(cvGF->examples[example].trans);
	for (i=0; i < edgeMatchGF->curMesh->numPts(); i++) {
		edgeMatchGF->vars[i*3 + 0] = cvGF->cv->tm->getPt(i)[0];
		edgeMatchGF->vars[i*3 + 1] = cvGF->cv->tm->getPt(i)[1];
		edgeMatchGF->vars[i*3 + 2] = cvGF->cv->tm->getPt(i)[2];
	}*/
}

void cvEdgeMatchSolve(const char *params) {
	if (!edgeMatchGF) {
		cout << "warning: goal function not initialized in cvEdgeMatchSolve" << endl;
		return;
	}
	if (solver) {
		cout << "warning: solver already running in cvEdgeMatchSolve" << endl;
		return;
	}

	params = extractDouble(params, &(edgeMatchGF->surfaceMatchWeight));
	params = extractDouble(params, &(edgeMatchGF->smoothnessWeight));
	params = extractDouble(params, &(edgeMatchGF->bendWeight));
	params = extractDouble(params, &(edgeMatchGF->markerMatchWeight));
	int maxIter = 1000;
	params = extractInt(params, &maxIter);
	int minV = 0, maxV = 10;
	params = extractInt(params, &minV);
	params = extractInt(params, &maxV);

	if (maxIter == -1) {
		runTest(edgeMatchGF, edgeMatchGF->vars, minV, maxV);
		return;
	}
/*
	double *vel = new double[edgeMatchGF->vars.size()];
	memset(vel, 0, sizeof(double)*edgeMatchGF->vars.size());
	int i, j;
	double timeStep = 0.0001;
	for (i=0; i < maxIter; i++) {
		edgeMatchGF->evaluateFunction(edgeMatchGF->vars);
		for (j=0; j < edgeMatchGF->vars.size(); j++) {
			// damping
			edgeMatchGF->grad[j] = -edgeMatchGF->grad[j] - 0.1 * vel[j];
			// clamping
			if (edgeMatchGF->grad[j] > 10)
				edgeMatchGF->grad[j] = 10;
			if (edgeMatchGF->grad[j] < -10)
				edgeMatchGF->grad[j] = -10;

			if (j < 10) {
				cout << j << ": " << edgeMatchGF->vars[j] << " " << vel[j] << " " << edgeMatchGF->grad[j] << endl;
			}
			edgeMatchGF->vars[j] += timeStep * vel[j] + 0.5 * timeStep * sqr(edgeMatchGF->grad[j]);
			vel[j] += timeStep + edgeMatchGF->grad[j];
			if (j < 10) {
				cout << j << ": " << edgeMatchGF->vars[j] << " " << vel[j] << " " << edgeMatchGF->grad[j] << endl;
			}
		}
	}
*/

/*	NonlinearConjugateGradientSolver *s = new NonlinearConjugateGradientSolver(edgeMatchGF);
	s->solve(maxIter, 0.0001, edgeMatchGF->vars);
	delete s;*/
	solver = new LBFGSSolver(edgeMatchGF);
	solver->solve(1e+3, 1e-7, edgeMatchGF->vars, maxIter);
	edgeMatchGF->applyDef(edgeMatchGF->vars);
	redrawV();
//	solver->solve(1e+3, 1e-5, edgeMatchGF->vars, maxIter);

	delete solver;
	solver = NULL;
}

void cvEdgeMatchMangle(const char *params) {
	int i;
	int mode = 0;
	params = extractInt(params, &mode);

	if (mode == 0) {
		Vec3d scale(1, 1, 1);
		params = extractDouble(params, &scale[0]);
		params = extractDouble(params, &scale[1]);
		params = extractDouble(params, &scale[2]);

		for (i=0; i < edgeMatchGF->vars.size(); i++) {
			edgeMatchGF->vars[i] *= scale[i%3];
		}
	}
	else if (mode == 2) {
		for (i=0; i < edgeMatchGF->vars.size(); i+=3) {
			edgeMatchGF->vars[i] += 0.5*edgeMatchGF->vars[i+1];
		}
	}
	else {
		double noise = 0.01;
		params = extractDouble(params, &noise);

		for (i=0; i < edgeMatchGF->vars.size(); i++) {
			edgeMatchGF->vars[i] += boundedRand(-noise, noise);
		}
	}

	edgeMatchGF->applyDef(edgeMatchGF->vars);

}

void cvEdgeMatchDiagnose(const char *params) {
	int i;
	int vert = 0;
	params = extractInt(params, &vert);

	for (i=0; i < edgeMatchGF->numEdges; i++) {
		if (edgeMatchGF->templateEdges[i].v0 == vert || 
			edgeMatchGF->templateEdges[i].v1 == vert ||
			edgeMatchGF->templateEdges[i].opp0 == vert || 
			edgeMatchGF->templateEdges[i].opp1 == vert) {
			cout << "Edge (" << edgeMatchGF->templateEdges[i].v0 << ", " << edgeMatchGF->templateEdges[i].v1 << "): current angle = " <<
				edgeMatchGF->templateEdges[i].eCurAngle << "; rest angle = " << edgeMatchGF->templateEdges[i].eRestAngle << endl;
		}
	}
}

void cvEdgeMatchLoadPly(const char *params) {
	int i, j;
	char fname[80], mapName[80];
	params = extractString(params, fname, 80);
	params = extractString(params, mapName, 80);

	TriMesh *tm = new TriMesh();
	tm->loadFile(fname);

	// open map file
	FILE *mapF = NULL;
	int mapSmall, mapBig;
	if (!openFile(&mapF, mapName, "rb", "mapping"))
		return;
	fread(&mapSmall, sizeof(int), 1, mapF);
	fread(&mapBig, sizeof(int), 1, mapF);
	if (mapBig != charVec.numPts) {
		cout << "invalid map file target size; expected " << charVec.numPts << "; found " << mapBig << endl;
		fclose(mapF);
		return;
	}
	if (mapSmall != tm->numPts()) {
		cout << "invalid map file; expected " << tm->numPts() << "; found " << mapSmall << endl;
		fclose(mapF);
		return;
	}

	for (i=0; i < mapBig; i++) {
		int numSources;
		fread(&numSources, sizeof(int), 1, mapF);

		edgeMatchGF->vars[i*3+0] = 0;
		edgeMatchGF->vars[i*3+1] = 0;
		edgeMatchGF->vars[i*3+2] = 0;

		for (j=0; j < numSources; j++) {
			int curV;
			float curW;
			fread(&curV, sizeof(int), 1, mapF);
			fread(&curW, sizeof(float), 1, mapF);

			edgeMatchGF->vars[i*3+0] += tm->getPt(curV)[0] * curW;
			edgeMatchGF->vars[i*3+1] += tm->getPt(curV)[1] * curW;
			edgeMatchGF->vars[i*3+2] += tm->getPt(curV)[2] * curW;
		}
	}

	edgeMatchGF->applyDef(edgeMatchGF->vars);

	fclose(mapF);
	delete tm;
}

void cvEdgeMatchLoadDisp(const char *params) {
	static TriMesh *smallMesh = NULL;
	double *smallDisp;
	int i, j;
	char dispName[80], mapName[80], meshName[80];
	mapName[0] = 0;
	meshName[0] = 0;
	params = extractString(params, dispName, 80);
	params = extractString(params, mapName, 80);
	params = extractString(params, meshName, 80);

/*
	vector<Vec3d> orig, trans;
	Mat4d m = Mat4d(0,2,0,1,2,0,0,3,0,0,2,4,0,0,0,1);
	for (i=0; i < 6; i++) {
		orig.push_back(Vec3d(boundedRand(-1,1),boundedRand(-1,1),boundedRand(-1,1)));
		trans.push_back(m*orig[i]);
	}
	cout << m << endl;
	cout << transformFromPts(orig, trans) << endl;
	return;*/

	if (meshName[0]) {
		if (smallMesh) {
			delete smallMesh;
			smallMesh = NULL;
		}
		smallMesh = new TriMesh();
		smallMesh->loadFile(meshName);
	}

	if (mapName[0]) {
		int ofs = 0;
		int expectedSize;

		FILE *mapF = NULL;
		int mapSmall, mapBig;
		if (!openFile(&mapF, mapName, "rb", "mapping"))
			return;
		fread(&mapSmall, sizeof(int), 1, mapF);
		fread(&mapBig, sizeof(int), 1, mapF);
		expectedSize = mapSmall * edgeMatchGF->varsPerVert;

		if (mapBig != charVec.numPts) {
			cout << "invalid map file; expected " << charVec.numPts << "; found " << mapBig << endl;
			fclose(mapF);
			return;
		}

		FILE *f;
		if (!openFile(&f, dispName, "rb", "displacements"))
			return;

		char version;
		fread(&version, sizeof(char), 1, f);
		if (version != '0') {
			cout << "unsupported version" << endl;
			return;
		}

		int fileSize;
		fread(&fileSize, sizeof(int), 1, f);

		double *fileDisps = new double[fileSize];
		baAssert(fileDisps != NULL, "can't allocate memory for displacements");
		fread(fileDisps, sizeof(double), fileSize, f);

		if (fileSize != expectedSize) {
			cout << "file has " << fileSize << " points; expected " << expectedSize << endl;
			fclose(f);
			return;
		}

		// first, we'll make the mapping
		Mat4d *mapping = new Mat4d[smallMesh->numPts()];
		vector<Vec3d> *origPts, *targetPts;
		origPts = new vector<Vec3d>[smallMesh->numPts()];
		targetPts = new vector<Vec3d>[smallMesh->numPts()];
		for (i=0; i < smallMesh->numPts(); i++) {
			origPts[i].push_back(smallMesh->getPt(i));
			targetPts[i].push_back(Vec3d(fileDisps[i*3+0], fileDisps[i*3+1], fileDisps[i*3+2]));
		}
		for (i=0; i < smallMesh->numTris(); i++) {
			for (j=0; j < 3; j++) {
				int v0 = smallMesh->getTri(i, j);
				int v1 = smallMesh->getTri(i, (j+1)%3);
				origPts[v0].push_back(smallMesh->getPt(v1));
				targetPts[v0].push_back(Vec3d(fileDisps[v1*3+0], fileDisps[v1*3+1], fileDisps[v1*3+2]));
			}
		}
		for (i=0; i < smallMesh->numPts(); i++) {
			mapping[i] = transformFromPts(origPts[i], targetPts[i]);
		}
		delete []origPts;
		delete []targetPts;

		// implement me!
		cout << "applying mapping from " << mapSmall << " to " << mapBig << endl;
		for (i=0; i < mapBig; i++) {
			int numSources;
			fread(&numSources, sizeof(int), 1, mapF);

			int curOfs = ofs + i * edgeMatchGF->varsPerVert;
			int k;
//			for (k=0; k < edgeMatchGF->varsPerVert; k++)
//				edgeMatchGF->vars[curOfs + k] = 0;
			Mat4d curM(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

			for (j=0; j < numSources; j++) {
				int curV;
				float curW;
				fread(&curV, sizeof(int), 1, mapF);
				fread(&curW, sizeof(float), 1, mapF);
				curM += curW * mapping[curV];
//				for (k=0; k < edgeMatchGF->varsPerVert; k++)
//					edgeMatchGF->vars[curOfs + k] += curW * fileDisps[curV * edgeMatchGF->varsPerVert + k];
			}

			Vec3d newV = curM * edgeMatchGF->origVerts[i];
			edgeMatchGF->vars[curOfs+0] = newV[0];
			edgeMatchGF->vars[curOfs+1] = newV[1];
			edgeMatchGF->vars[curOfs+2] = newV[2];
		}

		delete []fileDisps;
		delete []mapping;
		fclose(f);

		if (mapF)
			fclose(mapF);

		edgeMatchGF->applyDef(edgeMatchGF->vars);
	}
	else {
		cout << "unmapped load not implemented!" << endl;
	}
}

void cvSetConfidence(const char *params) {
	if (edgeMatchGF) {
		int i;
		for (i=0; i < edgeMatchGF->targetTM->numPts(); i++)
			edgeMatchGF->targetTM->getPtConf(i) = 1;
	}
	if (matchGF) {
		int i;
		for (i=0; i < matchGF->targetTM->numPts(); i++)
			matchGF->targetTM->getPtConf(i) = 1;
	}
}

void cvUnfold(const char *params) {
	int i;
	int reps = 1;
	params = extractInt(params, &reps);

	baAssert(edgeMatchGF != NULL, "edgeMatchGF not initialized!", false);
	edgeMatchGF->unfold(reps);
	redrawV();
}

void cvColorVert(const char *params) {
	int vert = 0;
	params = extractInt(params, &vert);
	double r = 0, g = 1, b = 0;
	params = extractDouble(params, &r);
	params = extractDouble(params, &g);
	params = extractDouble(params, &b);

	cvGF->cv->tm->getPtColor(vert) = Vec3d(r, g, b);
	redrawV();
}

void cvThing(const char *params) {
/*	char fname[80];
	params = extractString(params, fname, 80);

	skelUI->savePose(fname);
	*/
}

Vec3d *ePts = NULL;
int *eHeads, *eTails, eNumEdges, eNumPts;

void loadEdges(const char *params) {
	int i;
	ifstream in("edges.txt");

	in >> eNumPts >> eNumEdges;
	ePts = new Vec3d[eNumPts];
	for (i=0; i < eNumPts; i++)
		ePts[i] = Vec3d(boundedRand(-1, 1), boundedRand(-1, 1), 0); //boundedRand(-1, 1));
	eHeads = new int[eNumEdges];
	eTails = new int[eNumEdges];
	
	for (i=0; i < eNumEdges; i++)
		in >> eHeads[i] >> eTails[i];
	in.close();
}

void cvRender() {
	glDisable(GL_LIGHTING);
	int i;
	for (i=0; i < eNumEdges; i++) {
		glBegin(GL_LINES);
		glColor3d(0, 0, 1);
		ePts[eHeads[i]].glVertex();
		glColor3d(0, 1, 0);
		ePts[eTails[i]].glVertex();
		glEnd();
	}
}

void optEdges(const char *params) {
	int numIter, i, j, k;
	params = extractInt(params, &numIter);

	Vec3d *deltas = new Vec3d[eNumPts];

	for (i = 0; i < numIter; i++) {
		memset(deltas, 0, sizeof(Vec3d)*eNumPts);

//		for (j=0; j < eNumPts; j++) {
//			deltas[j] += -ePts[j] * 0.001;
//		}

		for (j=0; j < eNumPts; j++) {
			for (k=j+1; k < eNumPts; k++) {
				Vec3d delta = ePts[j] - ePts[k];
				double dist = (delta).length();
				delta = delta / dist;

				if (dist < 0.10) {
					deltas[j] += delta * 0.25*(0.10 - dist);
					deltas[k] -= delta * 0.25*(0.10 - dist);
				}
//				deltas[j] += -ePts[j] * 0.001;
			}
		}


		for (j=0; j < eNumEdges; j++) {
			Vec3d delta = ePts[eHeads[j]] - ePts[eTails[j]];
			double dist = (delta).length();
			delta = delta / dist;

			if (dist < 0.20) {
				deltas[eHeads[j]] += delta * sqr(0.20 - dist);
				deltas[eTails[j]] -= delta * sqr(0.20 - dist);
			}
			else if (dist > 0.25) {
				deltas[eHeads[j]] -= delta * 0.01 * (dist-0.25);
				deltas[eTails[j]] += delta * 0.01 * (dist-0.25);
			}
		}

		for (j=0; j < eNumPts; j++) {
			ePts[j] += deltas[j];
		}

		redrawVNow();
	}
}

void buildAnalogy(const char *params) {
	int i, j, k;
	char origFN[80], newFN[80], avgFN[80], outFN[80];
	char pcaFN[80], pcaNewFN[80];
	pcaFN[0] = 0; pcaNewFN[0] = 0;
	params = extractString(params, origFN, 80);
	params = extractString(params, newFN, 80);
	params = extractString(params, avgFN, 80);
	params = extractString(params, outFN, 80);
	params = extractString(params, pcaFN, 80);
	params = extractString(params, pcaNewFN, 80);

	TriMesh origMesh, newMesh, avgMesh;
	if (!origMesh.loadFile(origFN) || !newMesh.loadFile(newFN) || !avgMesh.loadFile(avgFN)) {
		cout << "can't open meshes; aborted" << endl;
	}

	int *mapBaryV = new int[newMesh.numPts() * 3];
	double *mapBary = new double[newMesh.numPts() * 3];
	Vec3d *mapOrigPts = new Vec3d[newMesh.numPts()];

	origMesh.calcHBB(16);
	origMesh.closestRestrictNormal = false;
	for (i=0; i < newMesh.numPts(); i++) {
		if (!origMesh.calcClosestPoint(newMesh.getPt(i), 1.10)) {
			cout << "no closest point for vert " << i << endl;
			mapBaryV[i*3+0] = 0; mapBaryV[i*3+1] = 0; mapBaryV[i*3+2] = 0;
			mapBary[i*3+0] = 0; mapBary[i*3+1] = 0; mapBary[i*3+2] = 0;
		}
		else {
			mapBaryV[i*3+0] = origMesh.closestTri[0];
			mapBaryV[i*3+1] = origMesh.closestTri[1];
			mapBaryV[i*3+2] = origMesh.closestTri[2];
			mapBary[i*3+0] = origMesh.closestBary[0];
			mapBary[i*3+1] = origMesh.closestBary[1];
			mapBary[i*3+2] = origMesh.closestBary[2];
			mapOrigPts[i] = origMesh.closestPt;
		}
	}

	// map average onto new mesh
	for (i=0; i < newMesh.numPts(); i++) {
		Vec3d newPt;
		for (j=0; j < 3; j++) {
			if (mapBaryV[i*3+j] > -1)
				newPt += avgMesh.getPt(mapBaryV[i*3+j]) * mapBary[i*3+j];
		}
		newMesh.getPt(i) += newPt - mapOrigPts[i];
	}
	newMesh.saveFile(outFN);

	FILE *fIn, *fOut;
	if (openFile(&fIn, pcaFN, "rb", "PCA file") && openFile(&fOut, pcaNewFN, "wb", "PCA file")) {
		int numPCA;
		fread(&numPCA, sizeof(int), 1, fIn);
		fread(&i, sizeof(int), 1, fIn);
		if (i != origMesh.numPts()*3) {
			cout << "number of points doesn't match; expected " << (origMesh.numPts()*3) << "; read in " << i << endl;
			return;
		}
		Vec3d **pcaData = new Vec3d*[numPCA];

		for (i=0; i < numPCA; i++) {
			pcaData[i] = new Vec3d[origMesh.numPts()];
			for (j=0; j < origMesh.numPts() * 3; j++) {
				fread(&pcaData[i][j/3][j%3], sizeof(double), 1, fIn);
			}
		}
		fclose(fIn);

		fwrite(&numPCA, sizeof(int), 1, fOut);
		i = newMesh.numPts() * 3;
		fwrite(&i, sizeof(int), 1, fOut);

		for (i=0; i < numPCA; i++) {
			for (j=0; j < newMesh.numPts() * 3; j++) {
				double v = 0;
				for (k=0; k < 3; k++) {
					if (mapBaryV[(j/3)*3+k] > -1)
						v += pcaData[i][mapBaryV[(j/3)*3+k]][j%3] * mapBary[(j/3)*3+k];
				}
				fwrite(&v, sizeof(double), 1, fOut);
			}
		}
		fclose(fOut);
	}

	delete []mapBaryV;
	delete []mapBary;
	delete []mapOrigPts;
}


/*
void cvTestt(const char *params) {
	ifstream in;
	int i;
	Skeleton skel;
	TriMesh mesh;

	skel.copyFrom(charVec.skel);
	if (!openIFStream(&in, "../posedata/ply/2_A.po.txt", "pose"))
		return;
	skel.loadPose(in);
	in.close();

	mesh.loadPly("../posedata/ply/2_A-matched.ply");
	mesh.calcNormals();

	scSkin.skel->copyVals(&skel);
	scSkin.skel->updateCoords();
	scSkin.updateLocalFrames();
	cvGF->examples[i].buildExample(&matchMeshes[i], &scSkin);

	setIn >> cvGF->examples[i].fname;
//		cvGF->examples[i].save(temp);

	cvGF->init(cvGF->numExamples, cvGF->numChars);

	cvGF->cv->updateCurComponents(cvGF->examples[cvGF->numExamples-1].dofs);
	cvGF->cv->updateLocalPos(NULL);
	cvGF->cv->updateTM();
}*/


void initCVMaster() {
	registerFunction(cvInit, "cvInit");
	registerFunction(cvInitFromSkin, "cvInitFromSkin");
	registerFunction(cvAutoUpdate, "cvAutoUpdate");
	registerFunction(cvLoadExampleSet, "cvLoadExampleSet");
	registerFunction(cvBuildExampleSet, "cvBuildExampleSet");
	registerFunction(cvRunEM, "cvRunEM");
	registerFunction(cvStopEM, "cvStopEM");
	registerFunction(cvElbowAnim, "cvElbowAnim");
	registerFunction(cvCalcRMS, "cvCalcRMS");
	registerFunction(cvSaveSet, "cvSaveSet");
	registerFunction(cvLoadSet, "cvLoadSet");
	registerFunction(cvMakeMovie, "cvMakeMovie");
	registerFunction(cvSavePly, "cvSavePly");

	registerFunction(cvSave, "cvSave");
	registerFunction(cvLoad, "cvLoad");
	registerFunction(cvLoadWeighted, "cvLoadWeighted");
	registerFunction(cvDumpPt, "cvDumpPt");
	registerFunction(cvDumpDOFs, "cvDumpDOFs");
	registerFunction(cvShowInf, "cvShowInf");

	registerFunction(cvLoadTexture, "cvLoadTexture");
	registerFunction(cvSetLocalFromTM, "cvSetLocalFromTM");
	registerFunction(cvSaveLocal, "cvSaveLocal");
	registerFunction(cvLoadLocal, "cvLoadLocal");
	registerFunction(cvSaveExample, "cvSaveExample");
	registerFunction(cvSetCVFromExample, "cvSetCVFromExample");
	registerFunction(cvSaveExMatrix, "cvSaveExMatrix");
	registerFunction(cvLoadPose, "cvLoadPose");

	registerFunction(cvInitSurfaceMatch, "cvInitSurfaceMatch");
	registerFunction(cvRunSurfaceMatch, "cvRunSurfaceMatch");
	registerFunction(cvStop, "cvStop");
	registerFunction(cvDumpVariance, "cvDumpVariance");

	registerFunction(cvAnimPoses, "cvAnimPoses");
	registerFunction(cvCreateMapping, "cvCreateMapping");
	registerFunction(cvSaveDisp, "cvSaveDisp");
	registerFunction(cvLoadDisp, "cvLoadDisp");
	registerFunction(cvSaveMirrorMap, "cvSaveMirrorMap");

	registerFunction(cvLoadMocap, "cvLoadMocap");
	registerFunction(cvMocapFrame, "cvMocapFrame");
	registerFunction(cvMocapAnim, "cvMocapAnim");
	registerFunction(cvSaveAni, "cvSaveAni");

	registerFunction(cvLoadMu, "cvLoadMu");

	registerFunction(cvAutoskinInit, "cvAutoskinInit");
	registerFunction(cvAutoskinSolve, "cvAutoskinSolve");
	registerFunction(cvSkinMatchInit, "cvSkinMatchInit");
	registerFunction(cvSkinMatchSolve, "cvSkinMatchSolve");

	registerFunction(cvEdgeMatchInit, "cvEdgeMatchInit");
	registerFunction(cvEdgeMatchSolve, "cvEdgeMatchSolve");
	registerFunction(cvEdgeMatchMangle, "cvEdgeMatchMangle");
	registerFunction(cvEdgeMatchDiagnose, "cvEdgeMatchDiagnose");
	registerFunction(cvEdgeMatchLoadPly, "cvEdgeMatchLoadPly");
	registerFunction(cvEdgeMatchLoadDisp, "cvEdgeMatchLoadDisp");
	registerFunction(cvSetConfidence, "cvSetConfidence");
	registerFunction(cvUnfold, "cvUnfold");
	registerFunction(cvColorVert, "cvColorVert");
	registerFunction(cvBrush, "cvBrush");
	registerFunction(cvSaveColor, "cvSaveColor");
	registerFunction(cvLoadColor, "cvLoadColor");
	registerFunction(cvLoadNeighWeights, "cvLoadNeighWeights");
	registerFunction(cvLoadSurfWeights, "cvLoadSurfWeights");

	registerFunction(loadEdges, "loadEdges");
	registerFunction(optEdges, "optEdges");
	registerFunction(buildAnalogy, "buildAnalogy");
}