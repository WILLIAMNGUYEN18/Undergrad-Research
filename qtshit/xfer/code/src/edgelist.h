#ifndef EDGELIST_H
#define EDGELIST_H

#include <vector>
using namespace std;
#include "trimesh.h"

struct EdgeInfo {
	int vert;
	int count;
	Vec3d normal;

	EdgeInfo() {
		count = 0;
		normal = Vec3d(0, 0, 0);
	}
};

// This class is just a list of every edge in a mesh.  We'll
// use it to find out if an edge appears twice.
class EdgeList {
public:
	int numVerts;
	typedef vector<EdgeInfo> edgeVec;
	edgeVec *edges;

	EdgeList() {
		edges = NULL;
	}

	void init(int nv) {
		numVerts = nv;
		if (edges)
			delete []edges;
		edges = new edgeVec[numVerts];
	}

	EdgeInfo *findEdge(int v0, int v1) {
		if (v1 < v0)
			swap(v1, v0);

		unsigned int i;
		for (i=0; i < edges[v0].size(); i++)
			if (edges[v0][i].vert == v1)
				return &edges[v0][i];
		return NULL;
	}

	void addEdge(int v0, int v1, Vec3d normal) {
		if (v1 < v0)
			swap(v1, v0);

		EdgeInfo newEI;
		bool isNew = false;
		EdgeInfo *ei = findEdge(v0, v1);
		if (ei == NULL) {
			ei = &newEI;
			isNew = true;
		}

		ei->vert = v1;
		ei->count++;
		ei->normal += normal;

		if (isNew)
			edges[v0].push_back(*ei);
	}

	void buildFromTriMesh(TriMesh &tm) {
		init(tm.numPts());
		int i;

		for (i=0; i < tm.numTris(); i ++) {
			addEdge(tm.getTri(i, 0), tm.getTri(i, 1), tm.getTriNormal(i));
			addEdge(tm.getTri(i, 1), tm.getTri(i, 2), tm.getTriNormal(i));
			addEdge(tm.getTri(i, 2), tm.getTri(i, 0), tm.getTriNormal(i));
		}
	}

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
	}
};

#endif