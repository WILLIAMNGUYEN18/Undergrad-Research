// trimesh.h
//
// A class for manipulating triangle meshes.
// by Brett Allen
//
// *** For license information, please see license.txt. ***

#ifndef _TRIMESH_H_
#define _TRIMESH_H_

#include <vector>
using namespace std;
#include "vec.h"
#include "mat.h"

typedef vector<int> intVec;
typedef vector<Vec3d> ptVec;

class HBBox;


int
intersect_triangle(Vec3d &orig, Vec3d &dir,
                   Vec3d &vert0, Vec3d &vert1, Vec3d &vert2,
                   double *t, double *u, double *v, bool *away);

// used for return-by-reference when no real reference is available:
extern Vec3d TEMP_VEC;
extern double TEMP_DOUBLE;


class TriMeshInterface {
public:
	virtual bool init(int pts = 0, int tris = 0);

	virtual int numPts() = 0;
	virtual bool setNumPts(int size) = 0;
	virtual void addPoint(const Vec3d &pt) = 0;

	virtual Vec3d &getPt(int index) = 0;
	virtual Vec3d &getPtColor(int index) { return TEMP_VEC; }
	virtual Vec3d &getPtNormal(int index) { return TEMP_VEC; }
	virtual double &getPtConf(int index) { return TEMP_DOUBLE; }

	virtual int numTris() = 0;
	virtual bool setNumTris(int size) = 0;
	virtual void addTri(int v0, int v1, int v2) = 0;

	virtual int &getTri(int face, int vert) = 0;
	virtual Vec3d &getTriColor(int face) { return TEMP_VEC; }
	virtual Vec3d &getTriNormal(int face) { return TEMP_VEC; }

	// The following functions get/set the existence of various
	// information.  If one or zero are passed in, then the flag 
	// will be set to true or false and the return value will 
	// indicate success or failure.  All other parameter values
	// return the current flag value without changing it.
	virtual bool gsPtColors(int action = -1) { return false; }
	virtual bool gsPtNormals(int action = -1) { return false; }
	virtual bool gsPtConfs(int action = -1) { return false; }
	virtual bool gsTriColors(int action = -1) { return false; }
	virtual bool gsTriNormals(int action = -1) { return false; }
};


class TriMesh : public TriMeshInterface {
protected:
	bool hasPtColors, hasPtConfs, hasTriColors;
	int addedPts, addedTris;

	ptVec pts;
	ptVec ptColors;
	ptVec ptNormals;
	vector<double> ptConfs;

	intVec tris;
	ptVec triColors;
	ptVec triNormals;

	int getPointIndex(Vec3d pt);
	HBBox *calcHBBRecursive(int maxSize, intVec &tris);

public:
/*	ptVec  m_points;
	ptVec  m_triNormals;
	ptVec  m_vertexNormals;
	ptVec  m_colors;
	intVec m_tris;
	vector<double> m_conf;*/
	double xClip;
	Vec3d solidColor;
	double alpha;
	bool showColor;

	HBBox  *m_hbb;
	
	TriMesh();
	TriMesh(TriMesh *t);
	~TriMesh();

	// TriMeshInterface functions:
	virtual bool init(int pts = 0, int tris = 0);
	virtual int numPts();
	virtual bool setNumPts(int size);
	virtual void addPoint(const Vec3d &pt);
	virtual Vec3d &getPt(int index);
	virtual Vec3d &getPtColor(int index);
	virtual Vec3d &getPtNormal(int index);
	virtual double &getPtConf(int index);
	virtual int numTris();
	virtual bool setNumTris(int size);
	virtual void addTri(int v0, int v1, int v2);
	virtual int &getTri(int face, int vert);
	virtual Vec3d &getTriColor(int face);
	virtual Vec3d &getTriNormal(int face);
	virtual bool gsPtColors(int action = -1);
	virtual bool gsPtNormals(int action = -1);
	virtual bool gsPtConfs(int action = -1);
	virtual bool gsTriColors(int action = -1);
	virtual bool gsTriNormals(int action = -1);

	void resetMesh();

	bool loadFile(char *fname);
	bool saveFile(char *fname);
	bool loadRaw(istream &in);
	bool loadTxt(istream &in);
	bool loadVrml(istream &in);
	bool loadObj(istream &in);
	bool saveObj(ostream &out);
	bool loadPly(const char *fn, bool color = false);
	bool savePly(const char *fn);
	void makeCylinder(double len, double rad, int lenSegs, int radSegs, double eRad = 0, double zMult = 1.0);
	void calcNormals();

	void copyPoints(TriMesh *source);

	void scale(double xs, double ys, double zs);
	void translate(double xt, double yt, double zt);
	void transform(Mat4d &m);

	void calcTriBounds(double *bb, intVec &tris);
	void calcPointBounds(double *bb, int *pts, int numPts);
	void calcPointBounds(double *bb);
	void dumpBounds();

	void calcHBB(int maxSize);

	Vec3d rayOrig, rayDir;
	double tPos, tNeg;
	bool hitPos, hitNeg;
	bool calcRayIntersection(Vec3d orig, Vec3d dir);
	void calcRayIntersection(intVec &tris);
	void calcRayIntersection(HBBox *hbb);
	double closestDist;
	Vec3d closestPt, closestBary;
	int closestTri[3];
	bool closestRestrictNormal;
	Vec3d closestNormalRestriction;
	bool calcClosestPoint(Vec3d orig, double closest, Vec3d normalMatch = Vec3d());
	bool calcClosestPoint(intVec &tris);
	bool calcClosestPoint(HBBox *hbb);

	void blur();
};

void calcFaceNormals(TriMeshInterface *tm);
void calcVertNormals(TriMeshInterface *tm);
void calcNormals(TriMeshInterface *tm);

#endif
