#include "trimesh.h"
#include "trimesh_util.h"
#include <float.h>
#include <stdio.h>

vector<int> *findTMNeighbors(TriMesh *tm) {
	vector<int> *neighbors;

	int i;

	neighbors = new vector<int>[tm->numPts()];
	for (i=0; i < tm->numTris(); i++) {
		int v0 = tm->getTri(i, 0);
		int v1 = tm->getTri(i, 1);
		int v2 = tm->getTri(i, 2);
		neighbors[v0].push_back(v1);
		neighbors[v1].push_back(v2);
		neighbors[v2].push_back(v0);
	}

	return neighbors;
}

vector<TMNeigh> *findTMGeodesics(TriMesh *tm, double maxDist, vector<int> *neighbors) {
	double *dists;
	bool *front;
	int pt, i;

	// make a neighbor map if it's not supplied
	bool delNeighbors = false;
	if (!neighbors) {
		neighbors = findTMNeighbors(tm);
		delNeighbors = true;
	}

	vector<TMNeigh> *geodesics = new vector<TMNeigh>[tm->numPts()];
	dists = new double[tm->numPts()];
	front = new bool[tm->numPts()];
	for (pt=0; pt < tm->numPts(); pt++) {
		front[pt] = false;
	}
	
	for (pt=0; pt < tm->numPts(); pt++) {
		for (i=0; i < tm->numPts(); i++) {
			dists[i] = DBL_MAX;
		}

		// build a distance map for the current point
		int curExpand = pt;
		dists[curExpand] = 0;
		while (curExpand >= 0) {
			// expand the current vertex
			for (i=0; i < neighbors[curExpand].size(); i++) {
				int n = neighbors[curExpand][i];
				double dist  = dists[curExpand] + (tm->getPt(curExpand) - tm->getPt(n)).length();
				if (dist < dists[n] && dist < maxDist) {
					dists[n] = dist;
					front[n] = true;
				}
			}
			front[curExpand] = false;

			// find the next vertex to expand
			curExpand = -1;
			double minExpand = DBL_MAX;
			for (i=0; i < tm->numPts(); i++) {
				if (front[i] && dists[i] < minExpand) {
					minExpand = dists[i];
					curExpand = i;
				}
			}
		}

		// store geodesic info for this point
		for (i=0; i < tm->numPts(); i++) {
			if (dists[i] < maxDist && i != pt) {
				geodesics[pt].push_back(TMNeigh(i, pt, dists[i]));
			}
		}
	}

	// delete the neighbor map if we made one
	if (delNeighbors)
		delete neighbors;

	return geodesics;
}

bool saveGeodesics(vector<TMNeigh> *geodesics, int numPts, const char *fname) {
	int pt, i;

	FILE *f = fopen(fname, "wb");
	if (!f)
		return false;

	fwrite(&numPts, sizeof(int), 1, f);
	for (pt=0; pt < numPts; pt++) {
		i = (int)geodesics[pt].size();
		fwrite(&i, sizeof(int), 1, f);
		for (i=0; i < geodesics[pt].size(); i++) {
			fwrite(&geodesics[pt][i], sizeof(TMNeigh), 1, f);
		}
	}

	fclose(f);
	return true;
}

bool loadGeodesics(vector<TMNeigh> *&geodesics, const char *fname) {
	int numPts, numN;
	int pt, i;

	FILE *f = fopen(fname, "rb");
	if (!f)
		return false;

	fread(&numPts, sizeof(int), 1, f);
	geodesics = new vector<TMNeigh>[numPts];
	for (pt=0; pt < numPts; pt++) {
		fread(&numN, sizeof(int), 1, f);
		geodesics[pt].resize(numN);
		for (i=0; i < numN; i++) {
			fread(&geodesics[pt][i], sizeof(TMNeigh), 1, f);
		}
	}

	fclose(f);
	return true;
}

struct EdgeInfo {
	int v0, v1;
	int leftTri, rightTri;
	int count;

	EdgeInfo() {
		count = 0;
		leftTri = -1;
		rightTri = -1;
	}
};

// This class is just a list of every edge in a mesh.  We'll
// use it to find out if an edge appears twice.
class TMUEdgeList {
public:
	int numVerts;
	typedef vector<EdgeInfo*> EdgeVec;
	EdgeVec allEdges, *vertEdges;

	TMUEdgeList() {
		vertEdges = NULL;
	}

	void init(int nv) {
		numVerts = nv;
		if (vertEdges)
			delete []vertEdges;
		vertEdges = new EdgeVec[numVerts];
	}

	EdgeInfo *findEdge(int v0, int v1) {
		if (v1 < v0)
			swap(v1, v0);

		unsigned int i;
		for (i=0; i < vertEdges[v0].size(); i++)
			if (vertEdges[v0][i]->v1 == v1)
				return vertEdges[v0][i];
		return NULL;
	}

	EdgeInfo *addEdge(int v0, int v1, int tri) {
		bool isNew = false;
		EdgeInfo *ei = findEdge(v0, v1);
		if (ei == NULL) {
			ei = new EdgeInfo;
			isNew = true;
		}

		if (v1 < v0) {
			swap(v0, v1);
			ei->rightTri = tri;
		}
		else {
			ei->leftTri = tri;
		}
		ei->v0 = v0;
		ei->v1 = v1;
		ei->count++;

		if (isNew) {
			vertEdges[v0].push_back(ei);
			allEdges.push_back(ei);
		}

		return ei;
	}

	void buildFromTriMesh(TriMesh &tm) {
		init(tm.numPts());
		int i;

		for (i=0; i < tm.numTris(); i ++) {
			addEdge(tm.getTri(i, 0), tm.getTri(i, 1), i);
			addEdge(tm.getTri(i, 1), tm.getTri(i, 2), i);
			addEdge(tm.getTri(i, 2), tm.getTri(i, 0), i);
		}
	}
/*
	void markVerts(char *verts) {
		memset(verts, 0, sizeof(char) * numVerts);
		unsigned int i, j;
		for (i=0; i < numVerts; i++) {
			if (edges[i].size() > 0) {
				for (j=0; j < edges[i].size(); j++) {
					if (edges[i][j].count == 1) {
						verts[i] = 1;
						verts[edges[i][j].vert] = 1;
					}
				}
			}
		}
	}*/
};

int thirdVert(TriMesh *mesh, int tri, int v0, int v1) {
	if (mesh->getTri(tri, 0) != v0 && mesh->getTri(tri, 0) != v1)
		return mesh->getTri(tri, 0);
	if (mesh->getTri(tri, 1) != v0 && mesh->getTri(tri, 1) != v1)
		return mesh->getTri(tri, 1);
	return mesh->getTri(tri, 2);
}

void butterflySubdiv(TriMesh *mesh, int **mappingVerts, double **mappingWeights) {
	TMUEdgeList edgeList;
	EdgeInfo **triEdges;
	int *triVerts;
	int tri, edge;
	int numTris = mesh->numTris();

	triEdges = new EdgeInfo*[numTris * 3];
	triVerts = new int[numTris * 3];

	edgeList.init(mesh->numPts());
	for (tri=0; tri < numTris; tri ++) {
		triVerts[tri*3+0] = mesh->getTri(tri, 0);
		triVerts[tri*3+1] = mesh->getTri(tri, 1);
		triVerts[tri*3+2] = mesh->getTri(tri, 2);

		triEdges[tri*3+0] = edgeList.addEdge(mesh->getTri(tri, 0), mesh->getTri(tri, 1), tri);
		triEdges[tri*3+1] = edgeList.addEdge(mesh->getTri(tri, 1), mesh->getTri(tri, 2), tri);
		triEdges[tri*3+2] = edgeList.addEdge(mesh->getTri(tri, 2), mesh->getTri(tri, 0), tri);
	}

	if (mappingVerts)
		*mappingVerts = new int[edgeList.allEdges.size() * 8];
	if (mappingWeights)
		*mappingWeights = new double[edgeList.allEdges.size() * 8];

	for (edge = 0; edge < edgeList.allEdges.size(); edge++) {
		edgeList.allEdges[edge]->count = mesh->numPts();
		// simple subdivision
		//mesh->addPoint(0.5 * (
		//	mesh->getPt(edgeList.allEdges[edge]->v0) +
		//	mesh->getPt(edgeList.allEdges[edge]->v1)));

		int v1 = edgeList.allEdges[edge]->v0;
		int v2 = edgeList.allEdges[edge]->v1;
		int v3, v4, v5, v6, v7, v8;
		EdgeInfo *e;

		v3 = thirdVert(mesh, edgeList.allEdges[edge]->rightTri, v1, v2);
		e = edgeList.findEdge(v1, v3);
		if (e->rightTri == edgeList.allEdges[edge]->rightTri)
			v5 = thirdVert(mesh, e->leftTri, v1, v3);
		else
			v5 = thirdVert(mesh, e->rightTri, v1, v3);
		e = edgeList.findEdge(v2, v3);
		if (e->rightTri == edgeList.allEdges[edge]->rightTri)
			v6 = thirdVert(mesh, e->leftTri, v2, v3);
		else
			v6 = thirdVert(mesh, e->rightTri, v2, v3);

		v4 = thirdVert(mesh, edgeList.allEdges[edge]->leftTri, v1, v2);
		e = edgeList.findEdge(v1, v4);
		if (e->rightTri == edgeList.allEdges[edge]->leftTri)
			v7 = thirdVert(mesh, e->leftTri, v1, v4);
		else
			v7 = thirdVert(mesh, e->rightTri, v1, v4);
		e = edgeList.findEdge(v2, v4);
		if (e->rightTri == edgeList.allEdges[edge]->leftTri)
			v8 = thirdVert(mesh, e->leftTri, v2, v4);
		else
			v8 = thirdVert(mesh, e->rightTri, v2, v4);

		double w = 1.0/16.0;
		mesh->addPoint(
			0.5 * (mesh->getPt(v1) + mesh->getPt(v2)) + 
			2.0 * w * (mesh->getPt(v3) + mesh->getPt(v4)) -
			w * (mesh->getPt(v5) + mesh->getPt(v6) + mesh->getPt(v7) + mesh->getPt(v8)));

		if (mappingVerts) {
			(*mappingVerts)[edge*8 + 0] = v1;
			(*mappingVerts)[edge*8 + 1] = v2;
			(*mappingVerts)[edge*8 + 2] = v3;
			(*mappingVerts)[edge*8 + 3] = v4;
			(*mappingVerts)[edge*8 + 4] = v5;
			(*mappingVerts)[edge*8 + 5] = v6;
			(*mappingVerts)[edge*8 + 6] = v7;
			(*mappingVerts)[edge*8 + 7] = v8;
		}
		if (mappingWeights) {
			(*mappingWeights)[edge*8 + 0] = 0.5;
			(*mappingWeights)[edge*8 + 1] = 0.5;
			(*mappingWeights)[edge*8 + 2] = 2.0 * w;
			(*mappingWeights)[edge*8 + 3] = 2.0 * w;
			(*mappingWeights)[edge*8 + 4] = -w;
			(*mappingWeights)[edge*8 + 5] = -w;
			(*mappingWeights)[edge*8 + 6] = -w;
			(*mappingWeights)[edge*8 + 7] = -w;
		}
	}

	mesh->setNumTris(numTris * 4);
	for (tri=0; tri < numTris; tri++) {
		mesh->addTri(triVerts[tri*3+0], triEdges[tri*3+0]->count, triEdges[tri*3+2]->count);
		mesh->addTri(triVerts[tri*3+1], triEdges[tri*3+1]->count, triEdges[tri*3+0]->count);
		mesh->addTri(triVerts[tri*3+2], triEdges[tri*3+2]->count, triEdges[tri*3+1]->count);
		mesh->addTri(triEdges[tri*3+0]->count, triEdges[tri*3+1]->count, triEdges[tri*3+2]->count);
	}
	mesh->calcNormals();
}

void findPtGeodesics(TriMesh *tm, int numSeeds, int *seedPts, double *seedDists, double maxDist, vector<TMNeigh> &geo, vector<int> *neighbors) {
	double *dists;
	bool *front;
	int pt, i;

	// make a neighbor map if it's not supplied
	bool delNeighbors = false;
	if (!neighbors) {
		neighbors = findTMNeighbors(tm);
		delNeighbors = true;
	}

	dists = new double[tm->numPts()];
	front = new bool[tm->numPts()];
	for (pt=0; pt < tm->numPts(); pt++) {
		front[pt] = false;
		dists[pt] = DBL_MAX;
	}
	
	for (pt=0; pt < numSeeds; pt++) {
		dists[seedPts[pt]] = seedDists[pt];
		front[seedPts[pt]] = true;
	}

	// build a distance map for the current point
	int curExpand;
	while (1) {
		// find the next vertex to expand
		curExpand = -1;
		double minExpand = DBL_MAX;
		for (i=0; i < tm->numPts(); i++) {
			if (front[i] && dists[i] < minExpand) {
				minExpand = dists[i];
				curExpand = i;
			}
		}
		if (curExpand == -1)
			break;

		// expand the current vertex
		for (i=0; i < neighbors[curExpand].size(); i++) {
			int n = neighbors[curExpand][i];
			double dist  = dists[curExpand] + (tm->getPt(curExpand) - tm->getPt(n)).length();
			if (dist < dists[n] && dist < maxDist) {
				dists[n] = dist;
				front[n] = true;
			}
		}
		front[curExpand] = false;
	}

	// store geodesic info for this point
	for (i=0; i < tm->numPts(); i++) {
		if (dists[i] < maxDist) {
			geo.push_back(TMNeigh(i, 0, dists[i]));
		}
	}
}