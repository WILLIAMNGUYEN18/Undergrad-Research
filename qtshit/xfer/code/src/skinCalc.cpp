// skinCalc.cpp

#include "ba.h"
#include "skinCalc.h"
#include "skeleton.h"
#include "trimesh.h"
#include "trimesh_util.h"
#include <float.h>
#include <list>

//char* skelName = "data/james.sk.txt";
char* skelName = "data/uni.sk.txt";

Vec3d *examplePoints;
Skeleton *matchSkel = NULL;
double *exampleDofs;
int numFrames, numDofs, numHandles;
vector<double> *infWeights;
vector<int> *infTrans;

TriMesh *scMesh;
Skeleton *scSkel;
Skin scSkin;

SkinCalcPt::SkinCalcPt() {
	free(false);
}

SkinCalcPt::~SkinCalcPt() {
	free();
}

double SkinCalcPt::getTransWeight(int t) {
	double ret = 0;
	int i;
	for (i=0; i < numTInf; i++) {
		if (tTransforms[i] == t) {
			ret = tWeights[i];
			break;
		}
	}
	return ret;
}


// Linear blend skinning ==================================

void LBSSkinCalcPt::free(bool delFirst) {
	if (delFirst) {
		if (tTransforms)
			delete []tTransforms;
		if (tWeights)
			delete []tWeights;
		if (tLocalPos)
			delete []tLocalPos;
	}
	numTInf = 0;
	tTransforms = NULL;
	tWeights = NULL;
	tLocalPos = NULL;
}

void LBSSkinCalcPt::initTrans(int iNumTInf) {
	numTInf = iNumTInf;
	tTransforms = new int[numTInf];
	tWeights = new double[numTInf];
	tLocalPos = new Vec3d[numTInf];

	numVars = numTInf * 4;
}

void LBSSkinCalcPt::initLocal(Skeleton *skel, Vec3d &v) {
	int j;
	for (j=0; j < numTInf; j++) {
		Mat4d m = skel->transforms.getT(tTransforms[j])->globalCoord.mat;
		m = m.inverse();
		tLocalPos[j] = m * v;
	}
}

void LBSSkinCalcPt::copyFromVars(double *vars) {
	int inf;
	int index = 0;

	for (inf=0; inf < numTInf; inf++) {
		tWeights[inf] = vars[index++];
		tLocalPos[inf][0] = vars[index++];
		tLocalPos[inf][1] = vars[index++];
		tLocalPos[inf][2] = vars[index++];
	}
}

void LBSSkinCalcPt::copyToVars(double *vars) {
	int inf;
	int index = 0;

	for (inf=0; inf < numTInf; inf++) {
		vars[index++] = tWeights[inf];
		vars[index++] = tLocalPos[inf][0];
		vars[index++] = tLocalPos[inf][1];
		vars[index++] = tLocalPos[inf][2];
	}
}

double LBSSkinCalcPt::calcGrad(double *grad, Vec3d &target, Skeleton *skel, double weight) {
	int inf;
	int index = 0;

	Vec3d delta = globalPos - target;
	double curErr = delta.length2();

	for (inf=0; inf < numTInf; inf++) {
		if (tWeights[inf] != 0) {
			double newErr;

			grad[index] += weight * 2.0 * delta *
				(skel->transforms.getT(tTransforms[inf])->globalCoord.mat * tLocalPos[inf]);
/*			tWeights[inf] += 0.0001;
			updateGlobal(skel);
			newErr = (globalPos - target).length2();
			grad[index] += (newErr - curErr) / 0.0001;
			tWeights[inf] -= 0.0001;
*/			index++;

			grad[index + 0] += weight * 2.0 * tWeights[inf] * 
				(delta * vec4to3(skel->transforms.getT(tTransforms[inf])->globalCoord.mat * Vec4d(1, 0, 0, 0)));
			grad[index + 1] += weight * 2.0 * tWeights[inf] * 
				(delta * vec4to3(skel->transforms.getT(tTransforms[inf])->globalCoord.mat * Vec4d(0, 1, 0, 0)));
			grad[index + 2] += weight * 2.0 * tWeights[inf] * 
				(delta * vec4to3(skel->transforms.getT(tTransforms[inf])->globalCoord.mat * Vec4d(0, 0, 1, 0)));
			index += 3;
		}
		else
			index += 4;
	}

	return weight * curErr;
}

void LBSSkinCalcPt::write(FILE *f) {
	// trim out unnecessary transforms
	int i, n;
	n = 0;
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] != 0)
			n++;
	}
	fwrite(&n, sizeof(int), 1, f);
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] != 0) {
			fwrite(&tTransforms[i], sizeof(int), 1, f);
		}
	}
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] != 0) {
			fwrite(&tWeights[i], sizeof(double), 1, f);
		}
	}
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] != 0) {
			fwrite(&tLocalPos[i].n, sizeof(double), 3, f);
		}
	}
}

void LBSSkinCalcPt::read(FILE *f) {
	fread(&numTInf, sizeof(int), 1, f);
	initTrans(numTInf);
	fread(tTransforms, sizeof(int), numTInf, f);
	fread(tWeights, sizeof(double), numTInf, f);
	fread(tLocalPos, sizeof(double)*3, numTInf, f);
}

void LBSSkinCalcPt::updateGlobal(Skeleton *skel) {
	globalPos = Vec3d();

	int inf;
	for (inf=0; inf < numTInf; inf++) {
		globalPos += tWeights[inf] * 
			((skel->transforms.getT(tTransforms[inf]))->globalCoord.mat * tLocalPos[inf]);
	}
}

void LBSSkinCalcPt::getPos(Vec3d &v) {
	v = globalPos;
}


// Multi-weight embedding skinning ========================

void MWESkinCalcPt::free(bool delFirst) {
	if (delFirst) {
		if (tTransforms)
			delete []tTransforms;
		if (tWeights)
			delete []tWeights;
		if (tMat)
			delete []tMat;
	}
	numTInf = 0;
	tTransforms = NULL;
	tWeights = NULL;
	tMat = NULL;
}

void MWESkinCalcPt::initTrans(int iNumTInf) {
	numTInf = iNumTInf;
	tTransforms = new int[numTInf];
	tWeights = new double[numTInf];
	tMat = new double12[numTInf];

	numVars = 3 + numTInf * 12;
}

void MWESkinCalcPt::initLocal(Skeleton *skel, Vec3d &v) {
	int i, j;
	for (j=0; j < numTInf; j++) {
		for (i=0; i < 12; i++)
			tMat[j][i] = tWeights[j];
	}

	updateGlobal(skel);
	localPos = localMat.inverse() * v;
}

void MWESkinCalcPt::copyFromVars(double *vars) {
	int inf, el;
	int index = 0;

	localPos[0] = vars[index++];
	localPos[1] = vars[index++];
	localPos[2] = vars[index++];

	for (inf=0; inf < numTInf; inf++) {
		for (el = 0; el < 12; el++) {
			tMat[inf][el] = vars[index++];
		}
	}
}

void MWESkinCalcPt::copyToVars(double *vars) {
	int inf, el;
	int index = 0;

	vars[index++] = localPos[0];
	vars[index++] = localPos[1];
	vars[index++] = localPos[2];

	for (inf=0; inf < numTInf; inf++) {
		for (el = 0; el < 12; el++) {
			vars[index++] = tMat[inf][el];
		}
	}
}

double MWESkinCalcPt::calcGrad(double *grad, Vec3d &target, Skeleton *skel, double weight) {
	int inf, el;
	int index = 0;

	Vec3d delta = localMat * localPos - target;
	double curErr = delta.length2();

	grad[index++] += weight * 2.0 * (delta * vec4to3(localMat * Vec4d(1, 0, 0, 0)));
	grad[index++] += weight * 2.0 * (delta * vec4to3(localMat * Vec4d(0, 1, 0, 0)));
	grad[index++] += weight * 2.0 * (delta * vec4to3(localMat * Vec4d(0, 0, 1, 0)));

	for (inf=0; inf < numTInf; inf++) {
		Mat4d &m = (skel->transforms.getT(tTransforms[inf]))->globalCoord.mat;
		if (tWeights[inf] != 0) {
			for (el=0; el < 12; el++) {
				grad[index++] += weight * 2.0 * delta[el/4] * m.n[el] * localPos[el%4];
			}
		}
		else
			index += 12;
	}

	return weight * curErr;
}

void MWESkinCalcPt::write(FILE *f) {
	fwrite(&localPos, sizeof(Vec3d), 1, f);

	// trim out unnecessary transforms
	int i, n;
	n = 0;
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] != 0)
			n++;
	}
	fwrite(&n, sizeof(int), 1, f);
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] != 0) {
			fwrite(&tTransforms[i], sizeof(int), 1, f);
		}
	}
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] != 0) {
			fwrite(&tWeights[i], sizeof(double), 1, f);
		}
	}
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] != 0) {
			fwrite(&tMat[i], sizeof(double), 12, f);
		}
	}
}

void MWESkinCalcPt::read(FILE *f) {
	fread(&localPos, sizeof(Vec3d), 1, f);

	fread(&numTInf, sizeof(int), 1, f);
	initTrans(numTInf);
	fread(tTransforms, sizeof(int), numTInf, f);
	fread(tWeights, sizeof(double), numTInf, f);
	fread(tMat, sizeof(double)*12, numTInf, f);
}

void MWESkinCalcPt::updateGlobal(Skeleton *skel) {
	localMat = Mat4d(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1);

	int inf, el;
	for (inf=0; inf < numTInf; inf++) {
		Mat4d &m = (skel->transforms.getT(tTransforms[inf]))->globalCoord.mat;
		for (el=0; el < 12; el++) {
			localMat.n[el] += tMat[inf][el] * m.n[el];
		}
	}
}

void MWESkinCalcPt::getPos(Vec3d &v) {
	v = localMat * localPos;
}


// Rot-Trans skinning =====================================

void RTSkinCalcPt::free(bool delFirst) {
	if (delFirst) {
		if (tTransforms)
			delete []tTransforms;
		if (tWeights)
			delete []tWeights;
		if (tPositions)
			delete []tPositions;
	}
	numTInf = 0;
	tTransforms = NULL;
	tWeights = NULL;
	tPositions = NULL;
}

void RTSkinCalcPt::initTrans(int iNumTInf) {
	numTInf = iNumTInf;
	tTransforms = new int[numTInf];
	tWeights = new double[numTInf];
	tPositions = new double[numTInf];

	numVars = 3 + numTInf * 2;
}

void RTSkinCalcPt::initLocal(Skeleton *skel, Vec3d &v) {
	updateGlobal(skel);
	localPos = localMat.inverse() * v;
}

void RTSkinCalcPt::copyFromVars(double *vars) {
	int inf;
	int index = 0;

	localPos[0] = vars[index++];
	localPos[1] = vars[index++];
	localPos[2] = vars[index++];

	for (inf=0; inf < numTInf; inf++) {
		tWeights[inf] = vars[index++];
		tPositions[inf] = vars[index++];
	}
}

void RTSkinCalcPt::copyToVars(double *vars) {
	int inf;
	int index = 0;

	vars[index++] = localPos[0];
	vars[index++] = localPos[1];
	vars[index++] = localPos[2];

	for (inf=0; inf < numTInf; inf++) {
		vars[index++] = tWeights[inf];
		vars[index++] = tPositions[inf];
	}
}

double RTSkinCalcPt::calcGrad(double *grad, Vec3d &target, Skeleton *skel, double weight) {
	double curErr = 0;
	/*
	int inf, el;
	int index = 0;

	Vec3d delta = localMat * localPos - target;
	double curErr = delta.length2();

	grad[index++] += 2.0 * (delta * vec4to3(localMat * Vec4d(1, 0, 0, 0)));
	grad[index++] += 2.0 * (delta * vec4to3(localMat * Vec4d(0, 1, 0, 0)));
	grad[index++] += 2.0 * (delta * vec4to3(localMat * Vec4d(0, 0, 1, 0)));

	for (inf=0; inf < numTInf; inf++) {
		Mat4d &m = (skel->transforms.getT(tTransforms[inf]))->globalCoord.mat;
		if (tWeights[inf] != 0) {
			for (el=0; el < 12; el++) {
				grad[index++] += 2.0 * delta[el/4] * m.n[el] * localPos[el%4];
			}
		}
		else
			index += 12;
	}
*/
	return weight * curErr;
}

void RTSkinCalcPt::write(FILE *f) {
	fwrite(&localPos, sizeof(Vec3d), 1, f);

	// trim out unnecessary transforms
	int i, n;
	n = 0;
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] != 0)
			n++;
	}
	fwrite(&n, sizeof(int), 1, f);
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] != 0) {
			fwrite(&tTransforms[i], sizeof(int), 1, f);
		}
	}
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] != 0) {
			fwrite(&tWeights[i], sizeof(double), 1, f);
		}
	}
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] != 0) {
			fwrite(&tPositions[i], sizeof(double), 1, f);
		}
	}
}

void RTSkinCalcPt::read(FILE *f) {
	fread(&localPos, sizeof(Vec3d), 1, f);

	fread(&numTInf, sizeof(int), 1, f);
	initTrans(numTInf);
	fread(tTransforms, sizeof(int), numTInf, f);
	fread(tWeights, sizeof(double), numTInf, f);
	fread(tPositions, sizeof(double), numTInf, f);
}

void RTSkinCalcPt::updateGlobal(Skeleton *skel) {
	int i, j;
	Vec3d trans;
	QuatNorm q;
	double weightSoFar = 0;
	for (i=0; i < numTInf; i++) {
		if (tWeights[i] == 0)
			continue;

		SkelTransform *t = skel->transforms.getT(tTransforms[i]);
		Vec3d pos = t->globalCoord.mat * (-t->curCoord.v * (1.0 - tPositions[i]));
		trans += pos * tWeights[i];
//		trans = pos;

		if (weightSoFar == 0)
			q = t->globalCoord.q;
		else {
			q = slerp(q, t->globalCoord.q, 1.0 - weightSoFar / (weightSoFar + tWeights[i]));
		}
		weightSoFar += tWeights[i];
	}
	localMat = Mat4d();
	localMat = q.toMatrixD();
	localMat[0][3] = trans[0];
	localMat[1][3] = trans[1];
	localMat[2][3] = trans[2];
}

void RTSkinCalcPt::getPos(Vec3d &v) {
	v = localMat * localPos;
}


// Skin ===================================================

void Skin::init(int iNumPts) {
	numPts = iNumPts;
	points = new POINT_TYPE[numPts];
	mirrorMap = NULL;
	baseVerts = new Vec3d[numPts];
}

void Skin::updatePoints() {
	int i;
	for (i=0; i < numPts; i++) {
		points[i].updateGlobal(skel);
	}
}

void Skin::updateMesh(TriMesh *tm) {
	int i;
	if (!baAssert(tm->numPts() == numPts, "Skin::updateLocalPos mismatched number of points", false))
		return;
	for (i=0; i < numPts; i++) {
		points[i].getPos(tm->getPt(i));
	}
	tm->calcNormals();
}

void Skin::renderPoints() {
/*	int i;
	for (i=0; i < numPts; i++) {
		points[i].updateLocalFrame(skel);
		Vec3d pt = points[i].localFrame * Vec3d(1,1,1);
		glColor3f(1, 1, 0);
		glbSphere(pt, 0.01, 10);
	}*/
}

void Skin::load(char *fname) {
	FILE *f;

	if (!openFile(&f, fname, "rb", "Skin"))
		return;

	fread(&numPts, sizeof(int), 1, f);
	init(numPts);
	fread(&numPreMirror, sizeof(int), 1, f);
	if (numPreMirror != numPts) {
		mirrorMap = new int[numPreMirror*2];
		fread(mirrorMap, sizeof(int), numPreMirror*2, f);
	}

	int i;
	for (i=0; i < numPts; i++)
		points[i].read(f);

	fclose(f);
}

void Skin::save(char *fname) {
	FILE *f;

	if (!openFile(&f, fname, "wb", "Skin"))
		return;

	fwrite(&numPts, sizeof(int), 1, f);
	fwrite(&numPreMirror, sizeof(int), 1, f);
	if (numPreMirror != numPts)
		fwrite(&mirrorMap, sizeof(int), numPreMirror*2, f);

	int i;
	for (i=0; i < numPts; i++)
		points[i].write(f);

	fclose(f);
}

void loadMarkerAnalysis(char *path, char *maFName) {
	int i, j, k, totalFrames, stride, index = 0;
	ifstream info;
	char str[1024], skelFN[80], markerFN[80], infFN[80], dofFN[80], *s;
	int numMarkers, firstFrame, lastFrame, *mhMapping;
	Vec3d v;

	// load info file
	sprintf(str, "%s%s", path, maFName);
	if (!openIFStream(&info, str, "Marker analysis")) {
		return;
	}
	info >> skelFN >> markerFN >> infFN >> dofFN;
	info >> numMarkers >> numHandles >> stride;
	info >> firstFrame >> lastFrame;
	mhMapping = new int[numMarkers];
	for (i=0; i < numMarkers; i++)
		mhMapping[i] = -1;
	while (info.good()) {
        info.getline(str, 1024);
		if (strlen(str) < 2)
			continue;
		i = atoi(strtok(str, ","));
		mhMapping[i] = atoi(strtok(NULL, ","));
	}
	info.close();

	// load skeleton
	sprintf(str, "%s%s", path, skelFN);
	matchSkel = Skeleton::load(str);

	// load marker data file
	sprintf(str, "%s%s", path, markerFN);
	if (!openIFStream(&info, str, "Marker file"))
		return;
	examplePoints = new Vec3d[numHandles * (lastFrame-firstFrame)];
	memset(examplePoints, 0, sizeof(Vec3d) * numHandles * (lastFrame-firstFrame));
	info.getline(str, 1024);
	strtok(str, "=");
	i = atoi(strtok(NULL, "="));
	if (i != numMarkers) {
		cout << "mismatch in number of markers in marker file; expected " << numMarkers << "; found " << i << endl;
		return;
	}
	info.getline(str, 1024);
	strtok(str, "=");
	totalFrames = atoi(strtok(NULL, "="));
	for (i=0; i < totalFrames; i++) {
		info.getline(str, 1024);
		info.getline(str, 1024);
		for (j=0; j < numMarkers; j++) {
			info >> k;
			if (k != j+1) {
				cout << "glitch reading marker file frame " << i << ": expected marker index " << (j+1) << "; got " << k << endl;
				return;
			}
			info >> v[0] >> v[1] >> v[2];
			if ((i % stride) == 0 && (i/stride) >= firstFrame && (i/stride) < lastFrame && mhMapping[j] > -1) {
				v *= 0.001;
				swap(v.n[1], v.n[2]);
				swap(v.n[0], v.n[2]);
				examplePoints[(i/stride - firstFrame) * numHandles + mhMapping[j]] = v;
			}
		}
		info.getline(str, 1024);
	}
	info.close();

	// load influence map
	sprintf(str, "%s%s", path, infFN);
	if (!openIFStream(&info, str, "DOF file"))
		return;
	infWeights = new vector<double>[numMarkers];
	infTrans = new vector<int>[numMarkers];
	while (1) {
		info >> i;
		if (i < 0 || !info.good())
			break;
		info >> j;
		for (k = 0; k < j; k++) {
			double weight;
			info >> weight;
			info >> str;
			int trans = matchSkel->transforms.lookupName(str);
			if (trans < 0) {
				cout << "unknown transform: " << str << endl;
				trans = 0;
			}
			infWeights[i].push_back(weight);
			infTrans[i].push_back(trans);
		}
	}
	info.close();

	// load dofs
	sprintf(str, "%s%s", path, dofFN);
	if (!openIFStream(&info, str, "DOF file"))
		return;
    info.getline(str, 1024);
	strtok(str, " ");
	strtok(NULL, " ");
	totalFrames = atoi(strtok(NULL, " "));
	strtok(NULL, " ");
	strtok(NULL, " ");
	numDofs = atoi(strtok(NULL, " "));
	cout << totalFrames << " frames; " << numDofs << " dofs" << endl;
	numFrames = lastFrame - firstFrame;
	exampleDofs = new double[numDofs * numFrames];
	info.getline(str, 1024); // dof names...
	for (i=0; i < totalFrames; i++) {
		for (j=0; j < numDofs; j++) {
			double d;
			info >> d;
			if (i >= firstFrame && i < lastFrame)
				exampleDofs[index++] = d;
		}
	}

	delete []mhMapping;

	setMatchSkelFrame(10);
}

void setMatchSkelFrame(int f) {
	if (!matchSkel)
		return;

	int ofs = numDofs * f;
	matchSkel->transforms.getT("pelvis_trans")->loadDofs(exampleDofs + ofs + 0);
	matchSkel->transforms.getT("pelvis_quat")->loadDofs(exampleDofs + ofs + 3);

	matchSkel->transforms.getT("l_thigh_quat")->loadDofs(exampleDofs + ofs + 7);
	matchSkel->transforms.getT("l_knee_euler_z")->loadDofs(exampleDofs + ofs + 11);
	matchSkel->transforms.getT("l_ankle_euler_z")->loadDofs(exampleDofs + ofs + 12);

	matchSkel->transforms.getT("r_thigh_quat")->loadDofs(exampleDofs + ofs + 13);
	matchSkel->transforms.getT("r_knee_euler_z")->loadDofs(exampleDofs + ofs + 17);
	matchSkel->transforms.getT("r_ankle_euler_z")->loadDofs(exampleDofs + ofs + 18);

	matchSkel->transforms.getT("abdomen_euler_x")->loadDofs(exampleDofs + ofs + 19);
	matchSkel->transforms.getT("abdomen_euler_x")->loadDofs(exampleDofs + ofs + 20);

	matchSkel->transforms.getT("spine_euler_y")->loadDofs(exampleDofs + ofs + 21);
	matchSkel->transforms.getT("head_quat")->loadDofs(exampleDofs + ofs + 22);

	double angle = exampleDofs[ofs + 26] - 0.9423;
	matchSkel->transforms.getT("l_scapula_euler_x")->loadDofs(&angle);
	matchSkel->transforms.getT("l_bicep_quat")->loadDofs(exampleDofs + ofs + 27);
	matchSkel->transforms.getT("l_elbow_euler_z")->loadDofs(exampleDofs + ofs + 31);
	matchSkel->transforms.getT("l_wrist_euler_x")->loadDofs(exampleDofs + ofs + 32);

	angle = exampleDofs[ofs + 33] + 0.9423;
	matchSkel->transforms.getT("r_scapula_euler_x")->loadDofs(&angle);
	matchSkel->transforms.getT("r_bicep_quat")->loadDofs(exampleDofs + ofs + 34);
	matchSkel->transforms.getT("r_elbow_euler_z")->loadDofs(exampleDofs + ofs + 38);
	matchSkel->transforms.getT("r_wrist_euler_x")->loadDofs(exampleDofs + ofs + 39);

//	Vec3d v(0.0000, -0.4018 * ((f+1) / 10.0), 0.0000);
//	matchSkel->transforms.getT("l_thigh_len")->loadDofs(v.n);

	matchSkel->updateCoords();
}

void renderMarkers(int frame) {
	int i, ofs = numHandles * frame;
	glColor3f(1, 1, 1);
	for (i=0; i < numHandles; i++) {
		glbSphere(examplePoints[ofs+i], 0.01);
	}
}

void saveSkelMatrices() {
	int hIndex, xyz, frame, joint;
	char fName[80];
	ofstream out, pOut;

	for (hIndex = 0; hIndex < numHandles; hIndex++) {
		if (infTrans[hIndex].size() < 2)
			continue;

		for (xyz = 0; xyz < 3; xyz++) {
			sprintf(fName, "data/skinmats/mat%04d-%d.txt", hIndex, xyz);
			if (!openOFStream(&out, fName, "matrix file"))
				return;
			sprintf(fName, "data/skinmats/vec%04d-%d.txt", hIndex, xyz);
			if (!openOFStream(&pOut, fName, "vector file"))
				return;

			for (frame=0; frame < numFrames; frame++) {
				setMatchSkelFrame(frame);

				for (joint = 0; joint < infTrans[hIndex].size(); joint++) {
					Mat4d m = matchSkel->transforms.getT(infTrans[hIndex][joint])->globalCoord.mat;
					out << m[xyz][0] << " " << m[xyz][1] << " " << m[xyz][2] << " " << m[xyz][3] << " ";
				}
				out << endl;

				pOut << examplePoints[frame * numHandles + hIndex][xyz] << endl;
			}
	
			out.close();
			pOut.close();
		}
	}
}

void meshKNNInterp(TriMesh *tm, int numPts, int *ptIndices) {
	double *dists = new double[tm->numPts() * numPts];
	
	vector<int> *neighbors = findTMNeighbors(tm);

	int pt, tPt, i;
	int ofs;

	for (pt=0; pt < numPts; pt++) {
		ofs = pt*tm->numPts();

		for (tPt=0; tPt < tm->numPts(); tPt++)
			dists[ofs + tPt] = 1e6;

		list <int> toSearch;
		list <int> toSearchD;
		toSearch.push_back(0);
		toSearchD.push_back(0);

		while (toSearch.size() > 0) {
			int curPt = toSearch.front();
			double curDist = toSearchD.front();
			toSearch.pop_front();
			toSearchD.pop_front();

			if (curDist < dists[ofs + curPt]) {
				dists[ofs + curPt] = curDist;

				Vec3d curV = tm->getPt(curPt);
				int neig;
				for (i = 0; i < neighbors[curPt].size(); i++) {
					toSearch.push_back(neighbors[curPt][i]);
					Vec3d neighV = tm->getPt(neighbors[curPt][i]);
					toSearchD.push_back(curDist + (curV - neighV).length());
				}
			}
		}
	}

	delete []dists;
}

#include "markers.h" 

void initSCMesh(char *fname, char *initPose, char *initMesh) {
	int i;
	scMesh = new TriMesh();

	/*
	MarkerSet *markers = new MarkerSet();
	markers->loadText("../doppel/james-data2/james-markers.txt");
	int i;
	ofstream f("../doppel/james-data2/james-new-markers.txt");
	f << markers->numMarkers << endl;
	Mat4d m(0.70488, 0.705371, 0.0748014, -0.086281, 
			-0.705544, 0.708088, -0.0286237, 0.0099844,
			-0.0731563, -0.0325995, 0.996788, 0.110465,
			0, 0, 0, 1);
//	m = m.inverse();
	for (i=0; i < markers->numMarkers; i++) {
		f << (vec4to3(m * vec3to4(markers->markers[i].pos)))	<< endl;
	}
	f.close();
	exit(0);*/

	if (strncmp(fname, "data/half", 9) == 0) {
		// mirror loaded surface
		scMesh->loadFile(fname);
		int n = scMesh->numPts();
		scSkin.numPreMirror = n;
		scSkin.mirrorMap = new int[n*2];
		int i;
		for (i=0; i < n; i++) {
			Vec3d &v = scMesh->getPt(i); 
/*			swap(v[1], v[2]);
			swap(v[0], v[1]);
			v[2] -= 1.0;
			v[0] = -v[0];
			v[1] = -v[1];*/
			if (v[1] < 0.001) {
				v[1] = 0;
				scSkin.mirrorMap[i] = i;
			}
			else {
				scSkin.mirrorMap[i] = scMesh->numPts();
				scSkin.mirrorMap[scMesh->numPts()] = i;
				scMesh->addPoint(Vec3d(v[0], -v[1], v[2]));
			}
		}
		n = scMesh->numTris();
		for (i=0; i < n; i++) {
			scMesh->addTri(
				scSkin.mirrorMap[scMesh->getTri(i, 0)], 
				scSkin.mirrorMap[scMesh->getTri(i, 2)], 
				scSkin.mirrorMap[scMesh->getTri(i, 1)]);
		}

		scMesh->savePly("test.ply");
//		scMesh->alpha = 0.5;

//		matchSkel = Skeleton::load(skelName);
		scSkel = Skeleton::load(skelName);

		// hack: use subject 2 for initialization
		if (initPose[0] != 0) {
			ifstream in(initPose);
			if (in.good()) {
				scSkel->loadPose(in);
				in.close();
				scSkel->updateCoords();
				TriMesh mesh;
				mesh.loadPly(initMesh);
				for (i=0; i < scMesh->numPts(); i++) {
					scMesh->getPt(i) = mesh.getPt(i);
				}
			}
			else 
				cout << "WARNING: can't open initialization pose" << endl;
		}
	}
	else {
		scSkin.numPreMirror = scSkin.numPts;
		scMesh->loadFile(fname);
		scSkel = Skeleton::load(skelName);
//		ifstream in;
//		if (!openIFStream(&in, "data/james.po.txt", "James pose"))
//			return;
//		scSkel->loadPose(in);
		//in.close();

		if (initPose[0] != 0) {
			ifstream in(initPose);
			if (in.good()) {
				scSkel->loadPose(in);
				in.close();
				scSkel->updateCoords();
				TriMesh mesh;
				mesh.loadPly(initMesh);
				int i;
				for (i=0; i < scMesh->numPts(); i++) {
					scMesh->getPt(i) = mesh.getPt(i);
				}
			}
			else 
				cout << "WARNING: can't open initialization pose" << endl;
		}
	}
	scMesh->calcNormals();

	scSkin.init(scMesh->numPts());
	scSkin.skel = scSkel;
	for (i=0; i < scMesh->numPts(); i++) {
		scSkin.baseVerts[i] = scMesh->getPt(i);
	}
}

static double *siTemp;
static int siSkelSize;
static vector<int> *siNeighbors;

void smoothSurfRecurse(double depth, double totalDepth, int index, int frame) {
	int tInd = index * siSkelSize + frame;
	double v = sqr(depth / totalDepth);
	if (siTemp[tInd] >= v)
		return;

	siTemp[tInd] = v;

	int i;
	double newDepth;
	for (i=0; i < siNeighbors[index].size(); i++) {
		int n = siNeighbors[index][i];

		if (n < 0)
			return;

		newDepth = depth - (scMesh->getPt(index) - scMesh->getPt(n)).length();

		if (newDepth > 0)
			smoothSurfRecurse(newDepth, totalDepth, n, frame);
	}
}

void smoothSurfInterp(int numFrames, int *minFrames, int *boneIndices) {
	siSkelSize = numFrames;
	siTemp = new double[scMesh->numPts() * siSkelSize];

	memset(siTemp, 0, sizeof(double) * scMesh->numPts() * siSkelSize);
	siNeighbors = findTMNeighbors(scMesh);

	int i, j;
	for (i=0; i < scMesh->numPts(); i++) {
		double recurseDepth = 0.20; // 0.03
		char *name = scSkel->transforms.getT(boneIndices[minFrames[i]])->name;
		if (//(strncmp(name, "chestT", 6) == 0) ||
			(strncmp(name+1, "Clav", 4) == 0) ||
			(strncmp(name+1, "UpperArm", 8) == 0) ||
			(strncmp(name, "abdomenT", 8) == 0))
			recurseDepth = 0.25;
		if ((strncmp(name+1, "UpperArm", 8) == 0))
			recurseDepth = 0.30;
/*		if ((strncmp(name+1, "Clav", 4) == 0) ||
			(strncmp(name, "chestT", 6) == 0) ||
			(strncmp(name+1, "Pelvis", 6) == 0) ||
			(strncmp(name, "abdomenT", 8) == 0))
			recurseDepth = 0.25;*/

//		recurseDepth = 0.05;
		
		smoothSurfRecurse(recurseDepth, recurseDepth, i, minFrames[i]);
	}

	cout << "done recursive assignment" << endl;

	// normalize
	for (i=0; i < scMesh->numPts(); i++) {
		int transInd = 0;
		double sum = 0;

		for (j=0; j < siSkelSize; j++) {
			sum += siTemp[i * siSkelSize + j];
		}

		for (j=0; j < siSkelSize; j++)
			scSkin.points[i].tWeights[j] = siTemp[i * siSkelSize + j] / sum;
	}

	delete []siTemp;
	delete []siNeighbors;
}

void calcTransInfluences() {
	int i, j, k;
	int numBones;
	int boneIndices[80];
	int *minFrames;

	// first off, determine the total number of bones
	numBones = 0;
	for (j=2; j < scSkel->transforms.size(); j++) {
		SkelTransform *curTrans = scSkel->transforms.getT(j);
		if (strcmp(curTrans->className, "SkelTranslation") != 0 
			&& strcmp(curTrans->className, "SkelSymmetricTranslation") != 0
			&& strcmp(curTrans->className, "SkelPartialTransform") != 0)
			continue;
		boneIndices[numBones] = j;
		numBones++;
	}

	// now, initialize each skin point to have a maximal # of influences...
	for (i=0; i < scSkin.numPts; i++) {
		scSkin.points[i].initTrans(numBones);

		for (j=0; j < numBones; j++) {
			scSkin.points[i].tWeights[j] = 0;
			scSkin.points[i].tTransforms[j] = boneIndices[j];
		}
	}

	// simple segmentation -- based on nearest bone
	minFrames = new int[scMesh->numPts()];
	for (i=0; i < scMesh->numPts(); i++) {
		Vec3d curPt = scMesh->getPt(i);
		Vec3d curNorm = scMesh->getPtNormal(i);
		curNorm.normalize();

//			searchMinSegment(0, 0, scSkel->transforms.getT(0)->globalCoord.v, scSkel);
		double minDist;
		int minFrame;
		minDist = DBL_MAX; 
		minFrame = 0;
		int boneIndex = 0;
		for (j=2; j < scSkel->transforms.size(); j++) {
			SkelTransform *curTrans = scSkel->transforms.getT(j);
			if (strcmp(curTrans->className, "SkelTranslation") != 0 
				&& strcmp(curTrans->className, "SkelSymmetricTranslation") != 0
				&& strcmp(curTrans->className, "SkelPartialTransform") != 0)
				continue;
			if ((strcmp(curTrans->name, "leftT") == 0) ||
				(strcmp(curTrans->name, "rightT") == 0) ||
				strcmp(curTrans->className, "SkelPartialTransform") == 0) {
				boneIndex++;
				continue;
			}

			double adjustFactor = 0.15;
			if (strcmp(curTrans->name, "neckT") == 0)
				adjustFactor = 0;
			if (strncmp(curTrans->name+1, "Hand", 4) == 0)
				adjustFactor = 0;
//			if (strncmp(curTrans->name+1, "Shin", 4) == 0)
//				adjustFactor = 0.05;
//			if (strncmp(curTrans->name+1, "Clav", 4) == 0)
//				adjustFactor = 0.2;
			if (strncmp(curTrans->name, "chestT", 6) == 0)
				adjustFactor = 0.05;
//			if (strncmp(curTrans->name+1, "UpperArm", 8) == 0)
//				adjustFactor = 0.5;
//			if (strncmp(curTrans->name+1, "Forearm", 7) == 0)
//				adjustFactor = 0.5;
			if (strncmp(curTrans->name+1, "Thigh", 5) == 0)
				adjustFactor = 0.25;
			if (strncmp(curTrans->name+1, "Shin", 4) == 0)
				adjustFactor = 0.0;
			if (strncmp(curTrans->name+1, "Foot", 4) == 0)
				adjustFactor = 0.0;

			Vec3d boneVec = curTrans->globalCoord.q.toMatrixD() * curTrans->curCoord.v;
			double boneLen = boneVec.length();
			boneVec /= boneLen;
			Vec3d distV = (curTrans->globalCoord.v - curPt);
			double linePos = (boneVec * distV);

#ifdef USE_POSITIONS
			scSkin.points[i].tPositions[boneIndex] = 1.0 - (linePos / boneLen);
			if (scSkin.points[i].tPositions[boneIndex] < 0)
				scSkin.points[i].tPositions[boneIndex] = 0;
			if (scSkin.points[i].tPositions[boneIndex] > 1)
				scSkin.points[i].tPositions[boneIndex] = 1;
#endif

			if (linePos < 0) {
				// distV already calculated
			}
			else if (linePos > boneLen) {
				distV = curPt - (curTrans->globalCoord.v - boneVec * boneLen);
			}
			else {
				distV = curPt - (curTrans->globalCoord.v - linePos * boneVec);
			}

			double d = distV.length();

			// only allow CLOSE matches on the arm
			if (((strncmp(curTrans->name+1, "UpperArm", 8) == 0) ||
				(strncmp(curTrans->name+1, "Forearm", 7) == 0)) &&
				(d > 0.1)) {
				d += 0.5;
			}
			// hell, let's be reasonable on the leg, too
			if (((strncmp(curTrans->name+1, "Thigh", 5) == 0) ||
				(strncmp(curTrans->name+1, "Shin", 4) == 0)) &&
				(d > 0.15)) {
				d += 0.5;
			}
			if (((strncmp(curTrans->name+1, "Foot", 4) == 0)) &&
				(d > 0.10)) {
				d += 0.5;
			}

			distV.normalize();
			double angleAdjust = distV * curNorm;
			angleAdjust += 0.75;
			if (angleAdjust > 0) {
				angleAdjust = adjustFactor * angleAdjust;
				d += angleAdjust;
			}

			if (d < minDist) {
				minDist = d;
				minFrame = boneIndex;
				minFrames[i] = boneIndex;
			}

			boneIndex++;
		}

//			double totalWeight = 0;
//			for (j=0; j < 3; j++)
//				if (minDists[j] < 1e6)
//					totalWeight += 1.0 / minDists[j];
		if (minDist < 1e6) {
			scSkin.points[i].tWeights[minFrame] = 1.0; //minDist;
//			scMesh->evalPts[i].transInterp[j] = minDists[j]; //(1.0 / minDists[j]) / totalWeight;
//			scMesh->evalPts[i].transFrame[j] = minFrames[j];
		}
	}

	smoothSurfInterp(numBones, minFrames, boneIndices);

	// copy the bone weights into the partial transforms
	for (i=0; i < scSkin.numPts; i++) {
		for (j=0; j < numBones; j++) {
			SkelTransform *curTrans = scSkel->transforms.getT(boneIndices[j]);
			SkelTransform *nextTrans;

			if (strcmp(curTrans->className, "SkelPartialTransform") != 0)
				continue;

			// find next transform
			for (k=j; k < numBones; k++) {
				nextTrans = scSkel->transforms.getT(boneIndices[k]);
				if (strcmp(nextTrans->className, "SkelPartialTransform") != 0)
					break;
			}

			if (scSkin.points[i].tWeights[k] > 0) {
				scSkin.points[i].tWeights[j] = 0.5 * scSkin.points[i].tWeights[k];
				scSkin.points[i].tWeights[k] = 0.5 * scSkin.points[i].tWeights[k];
			}
		}
	}

	// assign local positions
	for (i=0; i < scSkin.numPts; i++) {
		scSkin.points[i].initLocal(scSkin.skel, scMesh->getPt(i));
	}

	cout << "assigning colors..." << endl;
	scMesh->gsPtColors(1);
	scMesh->showColor = true;

	for (i=0; i < scSkin.numPts; i++) {
		Vec3d color;
		for (j=0; j < numBones; j++) 
			if (scSkin.points[i].tWeights[j] > 0)
				color += scSkin.points[i].tWeights[j] *
					scSkel->transforms.getT(scSkin.points[i].tTransforms[j])->color;
		scMesh->getPtColor(i) = color;
	}

	delete []minFrames;
}