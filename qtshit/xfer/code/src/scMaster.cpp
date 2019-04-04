#include <fstream>
using namespace std;

#include "doppel2.h"
#include "scMaster.h"
#include "skinCalc.h"
#include "cli.h"
#include "dof.h"
#include "skeleton.h"
#include "trimesh.h"

extern TriMesh *scMesh;
extern Skin scSkin;
Skin curSkin;

void scLoadMarkerAnalysis(const char *params) {
	char path[80], fName[80];
	params = extractString(params, path, 80);
	params = extractString(params, fName, 80);
	loadMarkerAnalysis(path, fName);
}

void scSaveMatrices(const char *) {
	saveSkelMatrices();
}

void scLoadSkin(const char *params) {
	char fName[80];
	params = extractString(params, fName, 80);
	scSkin.load(fName);
}

void scSaveSkin(const char *params) {
	char fName[80];
	params = extractString(params, fName, 80);
	scSkin.save(fName);
}

void scInitSCMesh(const char *params) {
	char fName[256], initPose[256], initMesh[256];
	initPose[0] = 0; initMesh[0] = 0;
	params = extractString(params, fName, 255);
	params = extractString(params, initPose, 255);
	params = extractString(params, initMesh, 255);

	initSCMesh(fName, initPose, initMesh);
}

void scCalcTransInfluences(const char *params) {
	calcTransInfluences();
}

void scUpdateLocal(const char *params) {
	scSkin.updatePoints();
	scSkin.updateMesh(scMesh);
}

void scUpdate(const char *params) {
	scSkel->updateCoords();
	scSkin.updatePoints();
	scSkin.updateMesh(scMesh);
	redrawV();
}

static Vec3d **ssdData = NULL;

void scInitSSD(const char *params) {
	int i, j;

	ssdData = new Vec3d*[scSkin.numPts];
	for (i=0; i < scSkin.numPts; i++) {
		ssdData[i] = new Vec3d[scSkin.points[i].numTInf];

		for (j=0; j < scSkin.points[i].numTInf; j++) {
			Mat4d trans = scSkel->transforms.getT(scSkin.points[i].tTransforms[j])->globalCoord.mat;
			trans = trans.inverse();
			ssdData[i][j] = trans * scMesh->getPt(i);
		}
	}
}

void scUpdateSSD(const char *params) {
	scSkel->updateCoords();
	int i, j;

	for (i=0; i < scSkin.numPts; i++) {
		scMesh->getPt(i) = Vec3d();

		for (j=0; j < scSkin.points[i].numTInf; j++) {
			Mat4d trans = scSkel->transforms.getT(scSkin.points[i].tTransforms[j])->globalCoord.mat;
			scMesh->getPt(i) += scSkin.points[i].tWeights[j] * (trans * ssdData[i][j]);
		}
	}
	scMesh->calcNormals();
	redrawV();
}

void scColor(const char *params) {
	int i, j;
	int inf = -1;
	params = extractInt(params, &inf);

	if (inf < 0) {
		for (i=0; i < scMesh->numPts(); i++)
			scMesh->getPtColor(i) = Vec3d(1, 1, 1);
	}
	else if (inf == 0) {
		for (i=0; i < scMesh->numPts(); i++) {
			Vec3d color;
			for (j=0; j < scSkin.points[i].numTInf; j++) 
				if (scSkin.points[i].tWeights[j] > 0)
					color += scSkin.points[i].tWeights[j] *
						scSkin.skel->transforms.getT(scSkin.points[i].tTransforms[j])->color;
			scMesh->getPtColor(i) = color;
		}
	}
	else {
		for (i=0; i < scMesh->numPts(); i++) {
			Vec3d color = Vec3d(0.5, 0.5, 0.5);
			for (j=0; j < scSkin.points[i].numTInf; j++) 
				if (scSkin.points[i].tTransforms[j] == inf && scSkin.points[i].tWeights[j] > 0) {
					color = Vec3d(1, 1.0 - scSkin.points[i].tWeights[j], 1.0 - scSkin.points[i].tWeights[j]);
				}
			scMesh->getPtColor(i) = color;
		}
	}
	redrawV();
}

void scClipPly(const char *params) {
	char fName[80], fName2[80];
	params = extractString(params, fName, 80);
	params = extractString(params, fName2, 80);

	TriMesh *tm = new TriMesh();
	if (!tm->loadFile(fName))
		return;

	int i, index = 0;
	int *vertMap = new int[tm->numPts()];
	TriMesh *tm2 = new TriMesh();
	for (i=0; i < tm->numPts(); i++) {
		Vec3d v = tm->getPt(i);
		if (v[1] < -0.01)
			vertMap[i] = -1;
		else {
			tm2->addPoint(v); // Vec3d(-v[2], -v[0], v[1]-1));
			vertMap[i] = index++;
		}
	}

	for (i=0; i < tm->numTris(); i++) {
		int v[3];
		v[0] = vertMap[tm->getTri(i, 0)];
		v[1] = vertMap[tm->getTri(i, 1)];
		v[2] = vertMap[tm->getTri(i, 2)];
		if (v[0] < 0 || v[1] < 0 || v[2] < 0)
			continue;
		tm2->addTri(v[0], v[1], v[2]);
	}
	tm2->saveFile(fName2);

	delete []vertMap;
	delete tm;
	delete tm2;
}

void scSaveWeights(const char *params) {
	int i, j;
	char fname[80];
	FILE *f;
	params = extractString(params, fname, 80);
	if (!openFile(&f, fname, "wb", "weight file"))
		return;
	int numT = scSkin.skel->transforms.size();
	double *weights = new double[numT];
	for (i=0; i < scSkin.numPts; i++) {
		for (j=0; j < numT; j++)
			weights[j] = -1;
		for (j=0; j < scSkin.points[i].numTInf; j++) {
			if (scSkin.points[i].tWeights[j] > 0) {
				weights[scSkin.points[i].tTransforms[j]] = scSkin.points[i].tWeights[j];
			}
		}
		fwrite(weights, sizeof(double), numT, f);
	}
	delete []weights;
	fclose(f);
}

void initSCMaster() {
	registerFunction(scLoadMarkerAnalysis, "scLoadMarkerAnalysis");
	registerFunction(scSaveMatrices, "scSaveMatrices");
	registerFunction(scSaveSkin, "scSaveSkin");
	registerFunction(scLoadSkin, "scLoadSkin");
	registerFunction(scInitSCMesh, "scInitSCMesh");
	registerFunction(scCalcTransInfluences, "scCalcTransInfluences");
	registerFunction(scUpdateLocal, "scUpdateLocal");
	registerFunction(scUpdate, "scUpdate");
	registerFunction(scColor, "scColor");

	registerFunction(scInitSSD, "scInitSSD");
	registerFunction(scUpdateSSD, "scUpdateSSD");

	registerFunction(scClipPly, "scClipPly");
	registerFunction(scSaveWeights, "scSaveWeights");
}