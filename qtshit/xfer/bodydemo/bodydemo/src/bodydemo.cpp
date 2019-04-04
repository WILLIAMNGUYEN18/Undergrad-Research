// bodydemo.cpp
// by Brett Allen (allen@cs.washington.edu)

// This tool demonstrates our editable body shape framework.
//
// There are two editing modes:
//   1. A slider-based mode where the user can manipulate certain
//      feature values (such as height, weight, etc.) of the current
//      body shape.
//   2. A point-dragging interface where the user can move any 
//      point on the body to any location, and the most-likely body
//      shape that meets those constraints is found.
//
// The two modes are independent (for now).  That is, changing the 
// body shape in one mode will not affect the body shape in the 
// other mode.
//
// The body shape is stored as a vector of weights, called curPCA.  
// The "average" person's vector is all zeros.  For the slider mode,
// the current vector is determined relative to a base vector, stored
// in curBase.

#include "main_win.h"
#include "bodydemo.h"
#include "viewer.h"
#include "ba.h"
#include "trimesh.h"
#include "pcaFitGF.h"
#include <GL/glu.h>

MainWin *mainWin;

FastTriMesh curMesh;

// PCA data
char80 setName;
const int maxPCA = 120;
int numReducedPCA;
float *pcaData, *avg, *variance;

// mode and selection
int curMode = SLIDER_MODE;
int selectedPt = -1;
int selectedPtInd = -1;
bool dragging = false;
bool coneDrag = false;
double conePosDelta = 0;
Vec3d coneDragPt, coneDragNorm;
bool showConstraints = true, showNormals = true, showLines = true;

// feature stuff
int curFeature;
int numFeatureMats;
char80 *fNames;
FeatureMat *fCirc;
int maxFeatures;
float *featureBase;
float *curPCA, *curBase;

float *markerPCA;
PCAFitGF *pcaGF;

// rendered size of indicators
const double POINT_DOT_RADIUS = 0.015;
const double CYLINDER_RADIUS = 0.005;
const double CYLINDER_LENGTH = 0.075;
const double CONE_RADIUS = 0.01;
const double CONE_LENGTH = 0.025;


// FeatureMat class =================================================
// This class encapsulates the data for editing a particular set of
// features (sliders).

// load a feature matrix (and its inverse) from text files
void FeatureMat::load(char *fname, char *invFName) {
	int i, j;
	ifstream in;
	if (!openIFStream(&in, fname, "feature matrix"))
		return;

	in >> numPCA >> numFeatures;
	numFeatures--;
	data = new float[numPCA*(numFeatures+1)];
	for (i=0; i < numPCA; i++) {
		for (j=0; j < numFeatures+1; j++) {
			in >> data[i * (numFeatures+1) + j];
		}
	}
	in.close();

	if (!openIFStream(&in, invFName, "inverse feature matrix"))
		return;
	invData = new float[numPCA*numFeatures];
	for (i=0; i < numFeatures; i++) {
		for (j=0; j < numPCA; j++) {
			in >> invData[i * numPCA + j];
		}
	}
	in.close();

	name = new char80[numFeatures];
	iName = new char80[numFeatures];
	minV = new float[numFeatures];
	maxV = new float[numFeatures];
	mult = new float[numFeatures];
	iMult = new float[numFeatures];
	power = new float[numFeatures];
	for (i=0; i < numFeatures; i++) power[i] = 1;
}

// look up the i'th entry for feature f in the feature matrix
float FeatureMat::v(int i, int f) {
	if (f < 0)
		return data[i * (numFeatures+1) + numFeatures];

	if (f >= numFeatures)
		return 0;

	return data[i * (numFeatures+1) + f];
}


// ==================================================================

inline Vec3d getPt(FastTriMesh *mesh, int pt) {
	int ofs = pt*3;
	return Vec3d(mesh->verts[ofs+0], mesh->verts[ofs+1], mesh->verts[ofs+2]);
}

inline Vec3d getPtNormal(FastTriMesh *mesh, int pt) {
	int ofs = pt*3;
	return Vec3d(mesh->norms[ofs+0], mesh->norms[ofs+1], mesh->norms[ofs+2]);
}


// calculate the likelihood of the curPCA values
void updateLikelihood() {
	double score = 0;
	int i;

	if (curMode == SLIDER_MODE) {
		for (i=0; i < numReducedPCA; i++) {
			score += sqr(curPCA[i]) / variance[i];
		}
		mainWin->scoreVO->value(max(0,10-(0.1*score)));
	}
	else {
		for (i=0; i < numReducedPCA; i++) {
			score += sqr(pcaGF->vars[i]) / variance[i];
		}
		mainWin->pScoreVO->value(max(0,10-(0.1*score)));
	}
}

void meshFromCurPCA() {
	// reconstruct the mesh from curPCA
	int i, j;

	memcpy(curMesh.verts, avg, curMesh.numPts * 3 * sizeof(float));
	float *pcaPtr = pcaData;
	for (j=0; j < numReducedPCA; j++) {
		double w = curPCA[j];

		if (fabs(w) < 1e-4) {
			// skip very small weights
			pcaPtr += curMesh.numPts*3;
			continue;
		}

		i = curMesh.numPts*3;
		int ind;
		for (ind = 0; ind < i; ind++) {
			curMesh.verts[ind] += w * *(pcaPtr++);
		}
	}
	curMesh.dirtyVerts = true;
	updateLikelihood();
	redrawV();
}

// given a set of feature (slider) values, update the mesh relative to curBase
void setFeatures(int mode, float *fVals) {
	int i, j;

	FeatureMat *f = &fCirc[mode];

	// update curPCA
	for (i=0; i < numReducedPCA; i++) {
		curPCA[i] = curBase[i];

		for (j=0; j < f->numFeatures; j++) {
			double val = fVals[j] - featureBase[j];
			curPCA[i] += f->v(i, j) * val;
		}
	}

	meshFromCurPCA();
}

// update the current feature (slider) values based on curBase
void updateFromBase() {
	int i, j;
	memset(featureBase, 0, sizeof(float)*maxFeatures);

	int numF = fCirc[curFeature].numFeatures;
	for (i=0; i < numF; i++) {
		for (j=0; j < numReducedPCA; j++) {
			featureBase[i] += fCirc[curFeature].invData[i*fCirc[curFeature].numPCA + j] * 
				(curBase[j] - fCirc[curFeature].data[j*(numF+1) + numF]);
		}
	}
	featureBase[i] = 1;
}

// make the current PCA values the base PCA values
void resetFromCur() {
	memcpy(curBase, curPCA, sizeof(float)*numReducedPCA);
	updateFromBase();
}

// set the current and base PCA values to be zero
void resetToAverage() {
	memset(curPCA, 0, sizeof(float)*numReducedPCA);
	resetFromCur();
	setFeatures(curFeature, featureBase);
}

// set the current and base PCA values randomly (according to a Gaussian distribution)
void resetToRandom() {
	int i;
	for (i=0; i < numReducedPCA; i++) {
		curPCA[i] = sqrt(-2 * log(boundedRand(0, 1) + 0.0000001)) * cos(boundedRand(0, PI*2)) * 0.8 * sqrt(variance[i]);
	}
	resetFromCur();
	setFeatures(curFeature, featureBase);
}

/* experimental code
void resetFromPly(char *fname) {
	int i, j;
	TriMesh mesh;
	if (!mesh.loadFile(fname))
		return;

	if (mesh.numPts() != curMesh.numPts) {
		cout << "point size mismatch when resetting from ply" << endl;
		return;
	}

	cout << "resetting from " << fname << endl;

	float *vec = new float[curMesh.numPts * 3];
	for (i=0; i < curMesh.numPts * 3; i++)
		vec[i] = mesh.getPt(i/3)[i%3];
	for (i=0; i < curMesh.numPts * 3; i++)
		vec[i] -= avg[i];
	for (j=0; j < numReducedPCA; j++) {
		curPCA[j] = 0;
		for (i=0; i < curMesh.numPts * 3; i++)
			curPCA[j] += vec[i] * pcaData[j*(curMesh.numPts*3) + i];
	}
	delete []vec;

	for (j=0; j < numReducedPCA; j++) {
		for (i=0; i < curMesh.numPts * 3; i++) {
			mesh.getPt(i/3)[i%3] -= curPCA[j] * pcaData[j*(curMesh.numPts*3) + i];
		}
	}
	mesh.savePly("reavg.ply");

	resetFromCur();
	setFeatures(curFeature, featureBase);
}

void applyDelta() {
	TriMesh tmM, tmF;
	if (!tmM.loadFile("data/ea-avg.ply") || !tmF.loadFile("data/ea-avg-f.ply"))
		return;

	int i;
	for (i=0; i < curMesh.numPts*3; i++) {
		avg[i] += tmF.getPt(i/3)[i%3] - tmM.getPt(i/3)[i%3];
	}
	for (i=0; i < curMesh.numPts; i++) {
		tmF.getPt(i) = Vec3d(avg[i*3+0], avg[i*3+1], avg[i*3+2]);
	}
	tmF.saveFile("data/boxer-f.ply");
}*/

// load a binary "points" file (used to initialize a dataset)
bool loadPoints(char *fname, float *&pts) {
	FILE *f;
	if (!openFile(&f, fname, "rb", "points"))
		return false;

	char version;
	fread(&version, sizeof(char), 1, f);

	if (version != '0') {
		cout << "error loading points " << fname << ": wrong version!" << endl;
		fclose(f);
		return false;
	}

	int fileSize;
	fread(&fileSize, sizeof(int), 1, f);
	if (fileSize != curMesh.numPts * 3) {
		cout << "error loading points " << fname << ": wrong number of points; expected " << (curMesh.numPts*3) << ", found " << fileSize << endl;
		fclose(f);
		return false;
	}

	pts = new float[fileSize];

	int i;
	double d;
	for (i=0; i < fileSize; i++) {
		fread(&d, sizeof(double), 1, f);
		pts[i] = d;
	}

	fclose(f);
	return true;
}

// load a binary "displacements" file (used to initialize a dataset)
double *loadDisp(char *fname) {
	FILE *f;
	if (!openFile(&f, fname, "rb", "displacements"))
		return NULL;

	char version;
	fread(&version, sizeof(char), 1, f);

	if (version != '0') {
		cout << "displacement version not supported in " << fname << endl;
		return NULL;
	}

	int fileSize;
	fread(&fileSize, sizeof(int), 1, f);

	double *fileDisps = new double[fileSize];
	baAssert(fileDisps != NULL, "can't allocate memory for displacements");
	fread(fileDisps, sizeof(double), fileSize, f);
	fclose(f);

	return fileDisps;
}

// load a data set in non-encapsulated form
void loadDataSet(char *fname) {
	ifstream dataSet, in;
	char s[256], s2[256];
	int i, j, index, pca;
	double d;
	FILE *f;
	TriMesh *tempMesh = new TriMesh();

	if (!openIFStream(&dataSet, fname, "data set")) {
		exit(0);
	}

	// load camera
	dataSet >> s;
	if (!openIFStream(&in, s, "camera"))
		return;
	in >> s;
	in >> mainWin->viewer->camera.rot.x >> mainWin->viewer->camera.rot.y >>
		mainWin->viewer->camera.rot.z >> mainWin->viewer->camera.rot.w;
	in >> s;
	in >> mainWin->viewer->camera.lightRot.x >> mainWin->viewer->camera.lightRot.y >>
		mainWin->viewer->camera.lightRot.z >> mainWin->viewer->camera.lightRot.w;
	in >> s;
	in >> mainWin->viewer->camera.trans[0] >> mainWin->viewer->camera.trans[1] >>
		mainWin->viewer->camera.trans[2];
	in.close();

	// load mesh
	dataSet >> s;
	if (!tempMesh->loadPly(s)) {
		cout << "can't load mesh " << s << endl;
		return;
	}
	curMesh.copyFromTriMesh(tempMesh);
	if (curMesh.colors) {
		for (i=0; i < curMesh.numPts * 3; i++) {
			curMesh.colors[i] = 0.8f;
		}
	}
	delete tempMesh;

	// load average
	dataSet >> s;
	if (s[0] == '-') {
		avg = new float[curMesh.numPts * 3];
		memcpy(avg, curMesh.verts, curMesh.numPts * 3 * sizeof(float));
	}
	else
		loadPoints(s, avg);

	// load reduced dimensions
	dataSet >> numReducedPCA;

	// load PCA data
	dataSet >> s;
	if (!openFile(&f, s, "rb", "pca data"))
		return;
	fread(&i, sizeof(int), 1, f); // num PCA vectors (I'll assume >= numReducedPCA)
	fread(&i, sizeof(int), 1, f);
	if (i != curMesh.numPts*3) {
		cout << "number of PCA points doesn't match; expected " << (curMesh.numPts * 3) << "; read in " << i << endl;
		return;
	}
	pcaData = new float[numReducedPCA * curMesh.numPts * 3];
	index = 0;
	for (i=0; i < numReducedPCA; i++) {
		for (j=0; j < curMesh.numPts * 3; j++) {
			fread(&d, sizeof(double), 1, f);
			pcaData[index++] = d;
		}
	}
	fclose(f);

	// load PCA variances
	dataSet >> s;
	if (!openIFStream(&in, s, "variances"))
		return;
	variance = new float[numReducedPCA];
	for (i=0; i < numReducedPCA; i++)
		in >> variance[i];
	in.close();

	// load feature matrices
	dataSet >> numFeatureMats;
	fCirc = new FeatureMat[numFeatureMats];
	fNames = new char80[numFeatureMats];
	maxFeatures = 0;
	for (i=0; i < numFeatureMats; i++) {
		dataSet >> fNames[i];
		dataSet >> s >> s2;
		fCirc[i].load(s, s2);

		if (fCirc[i].numFeatures > maxFeatures)
			maxFeatures = fCirc[i].numFeatures;

		for (j=0; j < fCirc[i].numFeatures; j++) {
			dataSet >> fCirc[i].name[j] >> fCirc[i].iName[j] >> 
				fCirc[i].mult[j] >> fCirc[i].iMult[j] >> 
				fCirc[i].power[j] >> fCirc[i].minV[j] >> fCirc[i].maxV[j];
		}
	}
	maxFeatures++;
}

// save the current dataset in encapsulated form
void saveDataSetBin(char *fname) {
	FILE *f;
	char s[256], s2[256];
	int i, j, index, pca;
	double d;

	if (!openFile(&f, fname, "wb", "data set"))
		return;

	// load camera
	fwrite(&mainWin->viewer->camera.rot, sizeof(QuatNorm), 1, f);
	fwrite(&mainWin->viewer->camera.lightRot, sizeof(QuatNorm), 1, f);
	fwrite(&mainWin->viewer->camera.trans, sizeof(Vec3d), 1, f);

	// load mesh
	i = curMesh.numPts;
	fwrite(&i, sizeof(int), 1, f);
	i = curMesh.numTris;
	fwrite(&i, sizeof(int), 1, f);
	fwrite(curMesh.tris, sizeof(int), i*3, f);

	// load average
	fwrite(avg, sizeof(float), curMesh.numPts*3, f);

	// load reduced dimensions
	fwrite(&numReducedPCA, sizeof(int), 1, f);

	// load PCA data
	fwrite(pcaData, sizeof(float), numReducedPCA * curMesh.numPts * 3, f);

	// load PCA variances
	fwrite(variance, sizeof(float), numReducedPCA, f);

	// load feature matrices
	fwrite(&numFeatureMats, sizeof(int), 1, f);
	for (i=0; i < numFeatureMats; i++) {
		j = (int)strlen(fNames[i]) + 1;
		fwrite(&j, sizeof(int), 1, f);
		fwrite(fNames[i], sizeof(char), j, f);
		
		fwrite(&fCirc[i].numFeatures, sizeof(int), 1 ,f);
		fwrite(fCirc[i].data, sizeof(float), numReducedPCA * (fCirc[i].numFeatures+1), f);
		for (j=0; j < fCirc[i].numFeatures; j++)
			fwrite(fCirc[i].invData + j*fCirc[i].numPCA, sizeof(float), numReducedPCA, f);

		for (j=0; j < fCirc[i].numFeatures; j++) {
			index = (int)strlen(fCirc[i].name[j])+1;
			fwrite(&index, sizeof(int), 1, f);
			fwrite(fCirc[i].name[j], sizeof(char), index, f);
			index = (int)strlen(fCirc[i].iName[j])+1;
			fwrite(&index, sizeof(int), 1, f);
			fwrite(fCirc[i].iName[j], sizeof(char), index, f);
			fwrite(&fCirc[i].mult[j], sizeof(float), 1, f);
			fwrite(&fCirc[i].iMult[j], sizeof(float), 1, f);
			fwrite(&fCirc[i].power[j], sizeof(float), 1, f);
			fwrite(&fCirc[i].minV[j], sizeof(float), 1, f);
			fwrite(&fCirc[i].maxV[j], sizeof(float), 1, f);
		}
	}
}

// load a dataset in encapsulated form
void loadDataSetBin(char *fname) {
	FILE *f;
	char s[256], s2[256];
	int i, j, index, pca;
	double d;

	if (!openFile(&f, fname, "rb", "data set"))
		return;

	// load camera
	fread(&mainWin->viewer->camera.rot, sizeof(QuatNorm), 1, f);
	fread(&mainWin->viewer->camera.lightRot, sizeof(QuatNorm), 1, f);
	fread(&mainWin->viewer->camera.trans, sizeof(Vec3d), 1, f);

	// load mesh
	fread(&i, sizeof(int), 1, f);
	fread(&j, sizeof(int), 1, f);
	curMesh.init(i, j, false);
	fread(curMesh.tris, sizeof(int), j*3, f);

	// load average
	avg = new float[curMesh.numPts*3];
	fread(avg, sizeof(float), curMesh.numPts*3, f);

	// load reduced dimensions
	fread(&numReducedPCA, sizeof(int), 1, f);

	// load PCA data
	pcaData = new float[numReducedPCA * curMesh.numPts * 3];
	fread(pcaData, sizeof(float), numReducedPCA * curMesh.numPts * 3, f);

	// load PCA variances
	variance = new float[numReducedPCA];
	fread(variance, sizeof(float), numReducedPCA, f);

	// load feature matrices
	fread(&numFeatureMats, sizeof(int), 1, f);
	fCirc = new FeatureMat[numFeatureMats];
	fNames = new char80[numFeatureMats];
	maxFeatures = 0;
	for (i=0; i < numFeatureMats; i++) {
		fread(&j, sizeof(int), 1, f);
		fread(fNames[i], sizeof(char), j, f);
		
		fCirc[i].numPCA = numReducedPCA;
		fread(&fCirc[i].numFeatures, sizeof(int), 1 ,f);
		if (fCirc[i].numFeatures > maxFeatures)
			maxFeatures = fCirc[i].numFeatures;
		fCirc[i].data = new float[numReducedPCA * (fCirc[i].numFeatures+1)];
		fread(fCirc[i].data, sizeof(float), numReducedPCA * (fCirc[i].numFeatures+1), f);
		fCirc[i].invData = new float[numReducedPCA * fCirc[i].numFeatures];
		fread(fCirc[i].invData, sizeof(float), numReducedPCA * fCirc[i].numFeatures, f);

		fCirc[i].name = new char80[fCirc[i].numFeatures];
		fCirc[i].iName = new char80[fCirc[i].numFeatures];
		fCirc[i].mult = new float[fCirc[i].numFeatures];
		fCirc[i].iMult = new float[fCirc[i].numFeatures];
		fCirc[i].power = new float[fCirc[i].numFeatures];
		fCirc[i].minV = new float[fCirc[i].numFeatures];
		fCirc[i].maxV = new float[fCirc[i].numFeatures];
		for (j=0; j < fCirc[i].numFeatures; j++) {
			fread(&index, sizeof(int), 1, f);
			fread(fCirc[i].name[j], sizeof(char), index, f);
			fread(&index, sizeof(int), 1, f);
			fread(fCirc[i].iName[j], sizeof(char), index, f);
			fread(&fCirc[i].mult[j], sizeof(float), 1, f);
			fread(&fCirc[i].iMult[j], sizeof(float), 1, f);
			fread(&fCirc[i].power[j], sizeof(float), 1, f);
			fread(&fCirc[i].minV[j], sizeof(float), 1, f);
			fread(&fCirc[i].maxV[j], sizeof(float), 1, f);
		}
	}
	maxFeatures++;
}

// switch between slider mode and point-dragging mode
void setCurMode(int m) {
	curMode = m;

	if (curMode == SLIDER_MODE) {
		meshFromCurPCA();
	}
	else if (curMode == MARKER_MODE) {
		pcaGF->applyDef(pcaGF->vars);
		curMesh.calcNormals();
	}
	else {
		meshFromCurPCA();
		mainWin->updatePCASliders();
	}
	
	updateLikelihood();
	redrawV();
}

// render the markers, lines, and normals for point-dragging mode
void drawMarkerStuff(float alpha) {
	int i;

	for (i=0; i < pcaGF->markerIndices.size(); i++) {
		Vec3d meshPt = getPt(&curMesh, pcaGF->markerIndices[i]);
		Vec3d normal = getPtNormal(&curMesh, pcaGF->markerIndices[i]);
		normal.normalize();

		if (showLines) {
			glDisable(GL_LIGHTING);
			glColor4f(1, 1, 1, alpha);
			glBegin(GL_LINES);
			pcaGF->markerPositions[i].glVertex();
			meshPt.glVertex();
			glEnd();
			glEnable(GL_LIGHTING);
		}

		if (showConstraints) {
			if (i == selectedPtInd)
				glColor4f(1, 0.75f, 0, alpha);
			else
				glColor4f(1, 0, 0, alpha);
			glbSphere(pcaGF->markerPositions[i], POINT_DOT_RADIUS, 10);
		}

		if (showNormals) {
			if (i == selectedPtInd && coneDrag) {
				glPushMatrix();
					glTranslated(coneDragPt[0], coneDragPt[1], coneDragPt[2]);
					glColor4f(0.0f, 0.0, 1.0f, alpha);
					glbDirectedCyl(coneDragNorm, conePosDelta, CYLINDER_RADIUS, CYLINDER_RADIUS);
					glTranslated(coneDragNorm[0] * conePosDelta, coneDragNorm[1] * conePosDelta, coneDragNorm[2] * conePosDelta);
					glColor4f(0.0f, 1.0f, 0.0f, alpha);
					glbDirectedCyl(coneDragNorm, CONE_LENGTH, CONE_RADIUS, 0);
				glPopMatrix();
			}
			else {
				glPushMatrix();
					glTranslated(meshPt[0], meshPt[1], meshPt[2]);
					glColor4f(0.0f, 0.0, 1.0f, alpha);
					glbDirectedCyl(normal, CYLINDER_LENGTH, CYLINDER_RADIUS, CYLINDER_RADIUS);
					glTranslated(normal[0] * CYLINDER_LENGTH, normal[1] * CYLINDER_LENGTH, normal[2] * CYLINDER_LENGTH);
					glColor4f(0.0f, 1.0f, 0.0f, alpha);
					glbDirectedCyl(normal, CONE_LENGTH, CONE_RADIUS, 0);
				glPopMatrix();
			}
		}
	}
}

// render everything
void bodydemoDraw() {
	// draw mesh
	glColor3f(0.8f, 0.8f, 0.8f);
	if (dispMode == VM_SMOOTH)
		curMesh.renderSmooth();
	else if (dispMode == VM_WIREFRAME)
		curMesh.renderWF(bkgColor[0], bkgColor[1], bkgColor[2]);
	else if (dispMode == VM_HIDDEN_LINE)
		curMesh.renderSmoothWF(true, bkgColor[0], bkgColor[1], bkgColor[2]);
	else
		curMesh.renderSmoothWF(false, bkgColor[0], bkgColor[1], bkgColor[2]);

	if (curMode == MARKER_MODE) {
		glDisable(GL_DEPTH_TEST);
		drawMarkerStuff(0.25);
		glEnable(GL_DEPTH_TEST);

		drawMarkerStuff(1);
	}
}


// ray intersection code ============================================
double tPos, tNeg;
bool hitPos, hitNeg;
Vec3d rayOrig, rayDir;

bool calcRayIntersection(FastTriMesh *mesh, Vec3d orig, Vec3d dir) {
	tPos = 1e6;
	tNeg = -1e6;
	hitPos = false;
	hitNeg = false;
	rayOrig = orig;
	rayDir = dir;

	int j;
	double u, v;
	double curDist;
	bool hit = false;
	bool curAway;

	for (j=0; j < mesh->numTris*3; j += 3) {
		curDist = 1e6;
		if (intersect_triangle(rayOrig, rayDir,
				getPt(mesh, mesh->tris[j + 2]),
				getPt(mesh, mesh->tris[j + 1]),
				getPt(mesh, mesh->tris[j + 0]),
				&curDist, &u, &v, &curAway)) {

			if (curDist < 0) {
				if (curDist > tNeg) {
					tNeg = curDist;
					hitNeg = curAway; //true;
				}
			}
			else {
				if (curDist < tPos) {
					tPos = curDist;
					hitPos = curAway; //true;
				}
			}
		}
	}

	return hitPos || hitNeg;
}

bool sphereIntersection(Vec3d pnt, Vec3d dir, Vec3d center, double radius) {
	double t;
	//Calculate the coefficients for the quadratic that defines where the line intersects the sphere
	double a = dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2];
	double b = 2*pnt[0]*dir[0] - 2*dir[0]*center[0] + 2*pnt[1]*dir[1] - 2*dir[1]*center[1] +
			   2*pnt[2]*dir[2] - 2*dir[2]*center[2];
	double c = pnt[0]*pnt[0] - 2*pnt[0]*center[0] + center[0]*center[0] +
			   pnt[1]*pnt[1] - 2*pnt[1]*center[1] + center[1]*center[1] +
			   pnt[2]*pnt[2] - 2*pnt[2]*center[2] + center[2]*center[2] - 
			   radius*radius;
	//if the determinant of the quadratic is greater than or equal to 0, there is an intersection, else not
	return b*b - 4*a*c >= 0;
}

// ==================================================================

// handle a mouse click
bool bodydemoClick(int x, int y, bool alt) {
	int i;

	dragging = false;
	coneDrag = false;
	conePosDelta = 0;

	if (curMode != MARKER_MODE)
		return false;

	Vec3d orig, direction;
	y = mainWin->viewer->viewPort[3] - y - 1;
	gluUnProject(GLdouble(x), GLdouble(y), GLdouble(0), mainWin->viewer->modelMatrix, 
		mainWin->viewer->projMatrix, mainWin->viewer->viewPort, &orig[0], &orig[1], &orig[2]);
	gluUnProject(GLdouble(x), GLdouble(y), GLdouble(1), mainWin->viewer->modelMatrix, 
		mainWin->viewer->projMatrix, mainWin->viewer->viewPort, &direction[0], &direction[1], &direction[2]);
	direction = -(direction - orig);
	direction.normalize();

	if (alt) {
		// check if a ray from the eye through the mouseclick intersects the object
		if (calcRayIntersection(&curMesh, orig, direction)){
			if (hitNeg) {
				double time = tNeg;
				Vec3d intersectPt = direction * time + orig;
				double nearestDist = -1;
				Vec3d nearestPt;

				// Check each point in the mesh to see if it is closer to the intersection than any others seen so far.
				for(i = 0; i < curMesh.numPts; i++){
					double dist = (getPt(&curMesh, i) - intersectPt).length();
					//calculate absolute value of dist.
					if(dist < 0){
						dist = -dist;
					}
					if(dist < nearestDist || nearestDist < 0){
						nearestDist = dist;
						nearestPt = getPt(&curMesh, i);
						selectedPt = i;
					}
				}

				selectedPtInd = -1;
				for (i=0; i < pcaGF->markerIndices.size(); i++) {
					if (pcaGF->markerIndices[i] == selectedPt) {
						selectedPtInd = i;
						break;
					}
				}
				if (selectedPtInd < 0) {
					selectedPtInd = (int)pcaGF->markerIndices.size();
					pcaGF->markerIndices.push_back(selectedPt);
					pcaGF->markerPositions.push_back(getPt(&curMesh, selectedPt));
				}
				dragging = true;

				bodydemoFit();
				return true;
			}
		}
	}
	else {
		// check for sphere intersections
		for (i = 0; i < pcaGF->markerPositions.size(); i++){
			if (sphereIntersection(orig, direction, pcaGF->markerPositions[i], POINT_DOT_RADIUS)) {
				selectedPt = pcaGF->markerIndices[i];
				selectedPtInd = i;
				dragging = true;
				return true;
			}
		}

		if (showNormals) {
			// check for cone intersections
			for (i = 0; i < pcaGF->markerPositions.size(); i++){
				Vec3d normal = getPtNormal(&curMesh, pcaGF->markerIndices[i]);
				normal.normalize();
				Vec3d center = getPt(&curMesh, pcaGF->markerIndices[i]) + (normal * (CYLINDER_LENGTH + (CONE_LENGTH / 2)));
				if (sphereIntersection(orig, direction, center, CONE_LENGTH / 2)){
					selectedPt = pcaGF->markerIndices[i];
					selectedPtInd = i;
					dragging = true;
					coneDrag = true;
					coneDragPt = getPt(&curMesh, pcaGF->markerIndices[i]);
					coneDragNorm = normal;
					return true;
				}
			}
		}
	}

	return false;
}

// handle a mouse drag
bool bodydemoDrag(int x, int y) {
	if (!dragging || curMode != MARKER_MODE)
		return false;

	GLdouble xDir;
	GLdouble yDir;
	GLdouble zDir;
	GLdouble selX;
	GLdouble selY;
	GLdouble selZ;
	y = mainWin->viewer->viewPort[3] - y - 1;
	Vec3d selPt = pcaGF->markerPositions[selectedPtInd];
	gluProject(selPt[0], selPt[1], selPt[2], mainWin->viewer->modelMatrix, mainWin->viewer->projMatrix, mainWin->viewer->viewPort, 
		&selX, &selY, &selZ);
	gluUnProject(GLdouble(x), GLdouble(y), selZ, mainWin->viewer->modelMatrix, 
		mainWin->viewer->projMatrix, mainWin->viewer->viewPort, &xDir, &yDir, &zDir);
	if (coneDrag) {
		Vec3d normal = coneDragNorm;
		conePosDelta = Vec3d(xDir - coneDragPt[0], yDir - coneDragPt[1], zDir - coneDragPt[2]) * normal;
		normal = (Vec3d(xDir - coneDragPt[0], yDir - coneDragPt[1], zDir - coneDragPt[2])* normal - CYLINDER_LENGTH) * normal;
		pcaGF->markerPositions[selectedPtInd] = coneDragPt + normal;
	}
	else {
		pcaGF->markerPositions[selectedPtInd] = Vec3d(xDir, yDir, zDir);
	}

	bodydemoFit();
	return true;
}

// handle a mouse-up event
void bodydemoRelease() {
	dragging = false;
	coneDrag = false;
	redrawV();
}

// fit the body to the current markers (for point-dragging mode)
void bodydemoFit() {
//	cout << pcaGF->markerIndices[0] << " " << pcaGF->markerPositions[0] << endl;
	pcaGF->zeroDeformation(false);
	pcaGF->solve(30);
	curMesh.calcNormals();
	updateLikelihood();
	redrawV();
}

// remove all constraints for point-dragging mode
void bodydemoClearPts() {
	pcaGF->markerPositions.clear();
	pcaGF->markerIndices.clear();
	selectedPt = -1;
	selectedPtInd = -1;
	bodydemoFit();
}

// set conformation value for point-dragging mode
void bodydemoSetConform(double f) {
	pcaGF->markerVariance = sqr(f);
	bodydemoFit();
}

// save the current mesh
void bodydemoSaveMesh(char *fname) {
	int i;
	ofstream out(fname);
	if (!out.good()) {
		cout << "can't open " << fname << endl;
		return;
	}

	out << "g default" << endl;
	for (i=0; i < curMesh.numPts; i++) {
		out << "v " << 
			curMesh.verts[i*3+0] << " " << 
			curMesh.verts[i*3+1] << " " <<
			curMesh.verts[i*3+2] << endl;
	}
	for (i=0; i < curMesh.numTris; i++) {
		out << "f " << 
			(curMesh.tris[i*3+0]+1) << " " << 
			(curMesh.tris[i*3+2]+1) << " " << 
			(curMesh.tris[i*3+1]+1) << endl;
	}
	out.close();
}

void updatePCAVal(int ind, double val) {
	ind--;
	if (ind >= 0 && ind < numReducedPCA) {
		curPCA[ind] = val;
		meshFromCurPCA();
	}
}

int main(int argc, char **argv) {
	// default dataset
	strcpy(setName, "male.dat");
//	strcpy(setName, "data/dataset-female.txt");
//	strcpy(setName, "male-ea.dat");
	if (argc > 1)
		strcpy(setName, argv[1]);

	Fl::visual(FL_DOUBLE|FL_RGB);
	dispMode = VM_SMOOTH;

	mainWin = new MainWin();
	mainWin->window->show();

	// load dataset (check extension to determine format)
	if (setName[strlen(setName)-2] == 'x' || setName[strlen(setName)-2] == 'X')
        loadDataSet(setName);
	else
		loadDataSetBin(setName);

	// initialize variables
	featureBase = new float[maxFeatures];
	curPCA = new float[numReducedPCA];
	curBase = new float[numReducedPCA];
	curFeature = 0;
	resetToAverage();
	PCAData data;
	data.activeComponents = numReducedPCA;
	data.numComponents = numReducedPCA;
	data.average = avg;
	data.components = pcaData;
	data.sigma2 = variance;
	pcaGF = new PCAFitGF(&curMesh, &data);

	if (argc > 2)
		saveDataSetBin(argv[2]);
	mainWin->updateSets();
	mainWin->updateSliders();

	// add help text for point-dragging mode
	Fl_Text_Buffer *buff = new Fl_Text_Buffer();
	mainWin->infoText->buffer(buff);
	buff->text("Alt-click a point to add a \nconstraint.  Click and drag a \nconstraint to change the \nbody shape.");

	return Fl::run();
}
