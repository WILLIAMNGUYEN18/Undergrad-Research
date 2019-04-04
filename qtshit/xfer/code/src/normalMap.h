#ifndef NORMAL_MAP_HDR
#define NORMAL_MAP_HDR

#include "trimesh.h"
#include "quatnorm.h"

const int NORMAL_MAP_W = 256;
const int NORMAL_MAP_H = 256;
const int NORMAL_MAP_SIZE = NORMAL_MAP_W * NORMAL_MAP_H;

class TexBary {
public:
	bool exists;
	int tri;
	double bary[3];
};

void calcTexBary(TriMesh *srcMesh, TexBary *texBary);
void calcTangentSpace(TriMesh *mesh, Mat3d **tangentSpace);

void buildNormalMap(TriMesh *srcMesh, TexBary *texBary, TriMesh *targetMesh, Vec3d *normalMap, 
				   Mat3d *tangentSpace);
void diffuseNormalMap(Vec3d *normalMap);
void symmetrifyNormalMap(Vec3d *normalMap);
void renderNormalMap(Vec3d *normalMap, unsigned char *tex, 
					 TexBary *texBary, Mat3d *tangentSpace, 
					 QuatNorm lightRot, Mat4d modelMat);

#endif