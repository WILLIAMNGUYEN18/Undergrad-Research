#ifndef FAST_TRIMESH
#define FAST_TRIMESH

typedef unsigned long indexType;

class TriMeshInterface;

class FastTriMesh {
public:
	indexType numPts, numTris;
	float *verts, *norms, *colors;
	indexType *tris;

	bool dirtyVerts, dirtyColors;

	FastTriMesh(indexType nPts = 0, indexType nTris = 0, bool hasColors = true);
	~FastTriMesh();

	void init(indexType nPts = 0, indexType nTris = 0, bool hasColors = true);
	void free();

	void copyFromTriMesh(TriMeshInterface *tm);
	void calcNormals();

	void uploadVerts();
	void uploadColors();
	void renderSmooth();
	void renderSmoothWF(bool wfOnly, float r, float g, float b);
	void renderWF(float r, float g, float b);
};

#endif
