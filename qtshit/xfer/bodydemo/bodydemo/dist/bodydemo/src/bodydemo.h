#ifndef BODY_DEMO_H
#define BODY_DEMO_H

#include "ba.h"
#include <math.h>
#include <iostream>
#include <vector>
using namespace std;
#include "vec.h"
#include "mat.h"
#include "quatnorm.h"
#include "trimesh_render.h"
#include "fast_trimesh.h"

class MainWin;
extern MainWin *mainWin;

const int SLIDER_MODE = 0;
const int MARKER_MODE = 1;

void redrawV();
void setCurMode(int m);
void bodydemoDraw();
bool bodydemoClick(int x, int y, bool alt);
bool bodydemoDrag(int x, int y);
void bodydemoRelease();
void bodydemoFit();
void bodydemoClearPts();
void bodydemoSetConform(double v);
void bodydemoSaveMesh(char *fname);


extern int dispMode;
extern Vec3d bkgColor;
extern bool showConstraints, showNormals, showLines;

extern FastTriMesh curMesh;

typedef char char80[80];

class FeatureMat {
public:
	int numPCA, numFeatures;
	float *data;
	float *invData;
	char80 *name, *iName;
	float *minV, *maxV, *mult, *iMult, *power;

	void load(char *fname, char *invFName);
	float v(int i, int f);
};


extern char80 setName;
extern int numFeatureMats;
extern char80 *fNames;

extern int curFeature;
extern int maxFeatures;
extern FeatureMat *fCirc;
extern float *featureBase;

void resetFromCur();
void resetToAverage();
void resetToRandom();
void setFeatures(int mode, float *fVals);

#endif