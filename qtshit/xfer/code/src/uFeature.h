#ifndef U_FEATURE_H
#define U_FEATURE_H

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

void initFeatures(int iNumPCA);
void loadFeatures(char *fname);
void setFeatures(float *fVals, double *pca);
void setFeatureBase(float *fVals, double *pca);

#endif