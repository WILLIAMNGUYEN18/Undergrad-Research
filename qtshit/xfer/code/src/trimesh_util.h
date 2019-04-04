#ifndef TRIMESH_UTIL_H
#define TRIMESH_UTIL_H

class TriMesh;

#include <vector>
using namespace std;

class TMNeigh {
public:
	int vert, vert2;
	double dist, weight;
	TMNeigh() { vert=-1; vert2 = -1; dist = 0; weight = 1.0; }
	TMNeigh(int iVert, int iVert2, double iDist, double iWeight = 1.0) { vert=iVert; vert2=iVert2; dist=iDist; weight = iWeight; }
};


vector<int> *findTMNeighbors(TriMesh *tm);
vector<TMNeigh> *findTMGeodesics(TriMesh *tm, double maxDist, vector<int> *neighbors = NULL);
void findPtGeodesics(TriMesh *tm, int numSeeds, int *pts, double *dists, double maxDist, vector<TMNeigh> &geo, vector<int> *neighbors = NULL);

bool saveGeodesics(vector<TMNeigh> *geodesics, int numPts, const char *fname);
bool loadGeodesics(vector<TMNeigh> *&geodesics, const char *fname);

void butterflySubdiv(TriMesh *mesh, int **mappingVerts = NULL, double **mappingWeights = NULL);

#endif