#ifndef FAST_TRIMESH
#define FAST_TRIMESH

typedef unsigned long indexType;

class TriMeshInterface;
class CCEvalMesh;

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
	void copyFromEvalMesh(CCEvalMesh *tm);
	void calcNormals();

	void uploadVerts();
	void uploadColors();
	void renderSmooth();
};

#endif
