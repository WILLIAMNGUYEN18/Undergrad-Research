#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include "trimesh.h"
#include "trimesh_render.h"

using namespace std;

#define EPSILON 0.000000001

void renderTriMesh(TriMesh *tm, int viewMode, Vec3d bkg) {
	int i, j;
	int numTris = tm->numTris();
	int pt[3];
	Vec3d col;


	if (viewMode & VM_SURF_ON) {
		if (viewMode & VM_WF_ON) {
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(1, 1);
		}

		if (viewMode & VM_SURF_COLOR_S) {
			tm->solidColor.glColor(tm->alpha);
		}
		else if (!((viewMode & VM_SURF_COLOR_F) || (viewMode & VM_SURF_COLOR_V))) {
			bkg.glColor();
		}

		if (!((viewMode & VM_SURF_FLAT) || (viewMode & VM_SURF_SMOOTH))) {
			glDisable(GL_LIGHTING);
		}
		else {
			glEnable(GL_LIGHTING);
		}

		glBegin(GL_TRIANGLES);
			for (i=0; i < numTris; i++) {
				pt[0] = tm->getTri(i, 0);
				pt[1] = tm->getTri(i, 1);
				pt[2] = tm->getTri(i, 2);

				if (viewMode & VM_SURF_FLAT) {
					(-tm->getTriNormal(i)).glNormal();
				}
				if ((viewMode & VM_SURF_COLOR_F)) {
					tm->getTriColor(i).glColor(tm->alpha);
				}

				for (j = 0; j < 3; j++) {
					if (viewMode & VM_SURF_SMOOTH) {
						(-tm->getPtNormal(pt[j])).glNormal();
					}
					if ((viewMode & VM_SURF_COLOR_V)) {
						tm->getPtColor(pt[j]).glColor(tm->alpha);
					}
					tm->getPt(pt[j]).glVertex();
				}
			}
		glEnd();

		if (viewMode & VM_WF_ON) {
			glDisable(GL_POLYGON_OFFSET_FILL);
		}
	}

	if (viewMode & VM_WF_ON) {
		if (viewMode & VM_WF_COLOR_S) {
			tm->solidColor.glColor(tm->alpha);
		}
		else if (!((viewMode & VM_WF_COLOR_F) || (viewMode & VM_WF_COLOR_V))) {
			(Vec3d(1,1,1)-bkg).glColor();
		}

		if (!((viewMode & VM_WF_FLAT) || (viewMode & VM_WF_SMOOTH))) {
			glDisable(GL_LIGHTING);
		}
		else {
			glEnable(GL_LIGHTING);
		}

		for (i=0; i < numTris; i++) {
			glBegin(GL_LINE_STRIP);
				pt[0] = tm->getTri(i, 0);
				pt[1] = tm->getTri(i, 1);
				pt[2] = tm->getTri(i, 2);

				if (viewMode & VM_WF_FLAT) {
					(-tm->getTriNormal(i)).glNormal();
				}
				if ((viewMode & VM_WF_COLOR_F)) {
					tm->getTriColor(i).glColor(tm->alpha);
				}

				for (j = 0; j < 4; j++) {
					if (viewMode & VM_WF_SMOOTH) {
						(-tm->getPtNormal(pt[j%3])).glNormal();
					}
					if ((viewMode & VM_WF_COLOR_V)) {
						tm->getPtColor(pt[j%3]).glColor(tm->alpha);
					}
					tm->getPt(pt[j%3]).glVertex();
				}
			glEnd();
		}
	}
}
