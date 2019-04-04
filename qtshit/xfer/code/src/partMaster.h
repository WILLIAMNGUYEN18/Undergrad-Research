#ifndef PART_MASTER
#define PART_MASTER

#include "tool.h"
#include "vec.h"
#include "mat.h"
#include "vl/VLd.h"
#include <vector>
using namespace std;

extern bool showPart, showTarget;

class TriMesh;

class PCAPart {
public:
	int numPts, numComponents;
	TriMesh *mesh;
	double *weights;
	Vecd pcaW;
	Mat4d trans;

	TVec mean, variance;
	vector<TVec> components;

	bool save(char *fname);
	bool load(char *fname);
	void updateMesh();
};

void initPartMaster();
void updatePartWeights(double *w, int numW);
void loadTargetMesh(const char *fname);

#endif
