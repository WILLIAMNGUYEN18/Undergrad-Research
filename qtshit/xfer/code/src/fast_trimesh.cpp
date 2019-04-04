#ifdef WIN32
#include <windows.h>
#endif
#include <stdlib.h>
#include <GL/gl.h>
#include "fast_trimesh.h"
#include "trimesh.h"
//#include "evalCC.h"

static FastTriMesh *lastDrawn = NULL;

FastTriMesh::FastTriMesh(indexType nPts, indexType nTris, bool hasColors) {
	verts = NULL;
	norms = NULL;
	colors = NULL;
	tris = NULL;
	init(nPts, nTris, hasColors);
}

FastTriMesh::~FastTriMesh() {
	free();
	dirtyVerts = true;
	dirtyColors = true;
}

void FastTriMesh::init(indexType nPts, indexType nTris, bool hasColors) {
	free();
	numPts = nPts;
	numTris = nTris;
	if (numPts > 0) {
		verts = new float[numPts*3];
		norms = new float[numPts*3];
		if (hasColors)
			colors = new float[numPts*3];
	}
	tris =  new indexType[numTris*3];

	dirtyVerts = true;
	dirtyColors = true;
}

void FastTriMesh::free() {
	if (verts) {
		delete []verts;
		verts = NULL;
	}
	if (norms) {
		delete []norms;
		norms = NULL;
	}
	if (colors) {
		delete []colors;
		colors = NULL;
	}
	if (tris) {
		delete []tris;
		tris = NULL;
	}
}

void FastTriMesh::copyFromTriMesh(TriMeshInterface *tm) {
	bool hasColors = tm->gsPtColors();
	init(tm->numPts(), tm->numTris(), hasColors);

	int i;
	for (i=0; i < numPts; i++) {
		Vec3d v = tm->getPt(i);
		verts[i*3 + 0] = v[0];
		verts[i*3 + 1] = v[1];
		verts[i*3 + 2] = v[2];

		if (hasColors) {
			Vec3d v = tm->getPt(i);
			colors[i*3 + 0] = v[0];
			colors[i*3 + 1] = v[1];
			colors[i*3 + 2] = v[2];
		}
	}

	for (i=0; i < numTris * 3; i++) {
		tris[i] = tm->getTri(i/3, i%3);
	}

	dirtyVerts = true;
	dirtyColors = true;
}

void FastTriMesh::copyFromEvalMesh(CCEvalMesh *tm) {
/*	bool hasColors = true;
	init(tm->numEvalPts, tm->numTriangles, hasColors);

	int i;
	for (i=0; i < numPts; i++) {
		Vec3d v = tm->evalPts[i].dispPos;
		verts[i*3 + 0] = v[0];
		verts[i*3 + 1] = v[1];
		verts[i*3 + 2] = v[2];

		if (hasColors) {
			Vec3d v = tm->evalPts[i].color;
			colors[i*3 + 0] = v[0];
			colors[i*3 + 1] = v[1];
			colors[i*3 + 2] = v[2];
		}
	}

	for (i=0; i < numTris * 3; i++) {
		tris[i] = tm->triangles[i/3].verts[2-i%3];
	}

	dirtyVerts = true;
	dirtyColors = true;*/
}

void FastTriMesh::uploadVerts() {
	if (verts)
		glVertexPointer(3, GL_FLOAT, 0, verts);
	if (norms)
		glNormalPointer(GL_FLOAT, 0, norms);
	dirtyVerts = false;
}

void FastTriMesh::uploadColors() {
	if (colors)
		glColorPointer(3, GL_FLOAT, 0, colors);
	dirtyColors = false;
}

inline float sqr(float x) {
	return x*x;
}

void FastTriMesh::calcNormals() {
	int i;
	Vec3d pt0, pt1, pt2, norm;
	int t0, t1, t2;

	memset(norms, 0, sizeof(float)*numPts*3);

	double a0, a1, a2, b0, b1, b2;
	double n0, n1, n2, nl;

	for (i=0; i < numTris*3; i+=3) {
		t0 = tris[i+0] * 3;
		t1 = tris[i+1] * 3;
		t2 = tris[i+2] * 3;

//		pt0 = Vec3d(verts[t0+0], verts[t0+1], verts[t0+2]);
//		pt1 = Vec3d(verts[t1+0], verts[t1+1], verts[t1+2]);
//		pt2 = Vec3d(verts[t2+0], verts[t2+1], verts[t2+2]);
//		norm = (pt1-pt0)^(pt2-pt0);
//		norm.normalize();

		a0 = verts[t1+0] - verts[t0+0];
		a1 = verts[t1+1] - verts[t0+1];
		a2 = verts[t1+2] - verts[t0+2];
		b0 = verts[t2+0] - verts[t0+0];
		b1 = verts[t2+1] - verts[t0+1];
		b2 = verts[t2+2] - verts[t0+2];
		n0 = a1*b2 - a2*b1;
		n1 = a2*b0 - a0*b2;
		n2 = a0*b1 - a1*b0;
		nl = sqrt(n0*n0 + n1*n1 + n2*n2);
		n0 /= nl;
		n1 /= nl;
		n2 /= nl;

		norms[t0+0] += n0;
		norms[t0+1] += n1;
		norms[t0+2] += n2;
		norms[t1+0] += n0;
		norms[t1+1] += n1;
		norms[t1+2] += n2;
		norms[t2+0] += n0;
		norms[t2+1] += n1;
		norms[t2+2] += n2;
	}

	// normalization
/*	for (i=0; i < numPts*3; i += 3) {
		float d = sqrt(sqr(norms[i+0]) + sqr(norms[i+1]) + sqr(norms[i+2]));
		norms[i+0] = 0; //= d;
		norms[i+1] = 0; //= d;
		norms[i+2] = 1; //= d;
	}*/
}

void FastTriMesh::renderSmooth() {
	if (lastDrawn != this) {
		dirtyVerts = true;
		dirtyColors = true;
		lastDrawn = this;
	}

	if (dirtyVerts) {
		calcNormals();
		uploadVerts();
	}

	if (dirtyColors)
		uploadColors();

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	if (colors)
		glEnableClientState(GL_COLOR_ARRAY);
	else
		glDisableClientState(GL_COLOR_ARRAY);

	glDrawElements(GL_TRIANGLES, numTris*3, GL_UNSIGNED_INT, tris);

/*
	if (dirtyVerts) {
		calcNormals();
	}
	glBegin(GL_TRIANGLES);
		int i, pt0, pt1, pt2;
		for (i=0; i < numTris; i++) {
			pt0 = tris[i*3+0] * 3;
			pt1 = tris[i*3+1] * 3;
			pt2 = tris[i*3+2] * 3;

			glNormal3f(norms[pt0+0], norms[pt0+1], norms[pt0+2]);
			if (colors)
				glColor3f(colors[pt0+0], colors[pt0+1], colors[pt0+2]);
			glVertex3f(verts[pt0+0], verts[pt0+1], verts[pt0+2]);

			glNormal3f(norms[pt1+0], norms[pt1+1], norms[pt1+2]);
			if (colors)
				glColor3f(colors[pt1+0], colors[pt1+1], colors[pt1+2]);
			glVertex3f(verts[pt1+0], verts[pt1+1], verts[pt1+2]);

			glNormal3f(norms[pt2+0], norms[pt2+1], norms[pt2+2]);
			if (colors)
				glColor3f(colors[pt2+0], colors[pt2+1], colors[pt2+2]);
			glVertex3f(verts[pt2+0], verts[pt2+1], verts[pt2+2]);
		}
	glEnd();*/
}
