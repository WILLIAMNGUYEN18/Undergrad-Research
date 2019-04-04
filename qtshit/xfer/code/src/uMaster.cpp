#include "doppel2.h"
#include "cli.h"
#include "uMaster.h"
#include "uSkin.h"
#include "uMocap.h"
#include "skeleton.h"
#include "markers.h"
#include "jactest.h"
#include "trimesh.h"
#include "trimesh_util.h"
#include "uniUI.h"
#include "propui.h"
#include "surfaceDef.h"
#include "cvgf.h"
#include "skinMatch.h"
#include "rotMatch.h"
#include "normalMap.h"
#include "vl/VLd.h"
#include "main_win.h"
#include "springmatch.h"
#include "uFeature.h"
#include <float.h>
#include <GL/glu.h>

extern UniUI *uniUI;

USkin skin;
UDataSet dataSet;
USolver uSolver;
UMocap mocap;
UMocapPoses mocapPoses;
TriMesh *uMesh, *dispMesh = NULL;
vector<int> *uMeshNeigh;
LBFGSSolver *lbfgs = NULL;
int visSkel = 0;
bool showData = true;
bool showDress = false;
bool showGeodesics = false;
bool showMocap = true;
bool showMarkers = false;
bool showFrames = false;
bool showTex = false; //true;
bool showTangentSpace = false;
bool showMarkerSprings = false;
bool showAxes = false;
bool softwareNM = false;
int showMesh = -1;
static SkinMatchGF *edgeMatchGF = NULL;
double curComps[NUM_CUR_COMPS];
char80 *markerNames;
vector<Vec3d> extraMarkers;
int pddWJoint = -1;

Vec3d *normalMap = NULL;
TexBary *texBary = NULL;
unsigned char *curTex = NULL;
Mat3d *tangentSpace = NULL;

// skin painting stuff
int spTransform = 0, spMode = 0;
double spIntensity = 1.0, spInnerRadius = 0.01, spOuterRadius = 0.01;
bool spDragging = false, spAutoUpdate = false, spAutoExtend = false, spGeodesic = true;

extern double spandex;

class KNNStruct {
public:
	int numPts, k;
	double *dists;
	int *ids;

	KNNStruct(int iNumPts, int iK) {
		numPts = iNumPts;
		k = iK;
		dists = new double[numPts * (k+1)];
		ids = new int[numPts * (k+1)];
		int i;
		for (i=0; i < numPts * (k+1); i++) {
			dists[i] = DBL_MAX;
			ids[i] = -1;
		}
	}

	~KNNStruct() {
		delete []dists;
		delete []ids;
	}

	bool insertDist(int pt, double dist, int id) {
		int i, j;
		int oldPos = -1;

		// first, determine if this id is already here
		for (i=0; i < k+1; i++) {
			if (ids[pt*(k+1) + i] == id) {
				oldPos = i;
				break;
			}
		}
		if (oldPos > -1) {
			if (dist >= dists[pt*(k+1) + oldPos])
				return false;
		}

		for (i=0; i < k+1; i++) {
			if (dist < dists[pt*(k+1) + i]) {
				// push others down
				if (oldPos > -1) {
					for (j=oldPos; j > i; j--) {
						dists[pt*(k+1) + j] = dists[pt*(k+1) + j-1];
						ids[pt*(k+1) + j] = ids[pt*(k+1) + j-1];
					}
				}
				else {
					for (j=k; j > i; j--) {
						dists[pt*(k+1) + j] = dists[pt*(k+1) + j-1];
						ids[pt*(k+1) + j] = ids[pt*(k+1) + j-1];
					}
				}
				// add this point
				dists[pt*(k+1) + i] = dist;
				ids[pt*(k+1) + i] = id;
				return true;
			}
		}
		return false;
	}

	void getWeights(int pt, int *wIds, double *wWeights) {
		int i;

		for (i=0; i < k; i++) {
			wIds[i] = -1;
			wWeights[i] = 0;
		}

		if (dists[pt*(k+1) + 0] < 1e-5) {
			wIds[0] = dists[pt*(k+1) + 0];
			wWeights[0] = 1;
		}
		else {
			double cur, sum = 0;
			int last;
			for (last=k; last > 0; last--)
				if (dists[pt*(k+1) + last] != DBL_MAX) break;
			double thresh = dists[pt*(k+1) + last];
			if (dists[pt*(k+1) + 0] - thresh <= 1e-5) {
				for (i=0; i < last; i++) {
					wWeights[i] = 1;
					sum++;
					wIds[i] = ids[pt*(k+1) + i];
				}
			}
			else {
				for (i=0; i < last; i++) {
					cur = -(1.0 - (dists[pt*(k+1) + i] / thresh)); // / (dist[i].dist * dist[i].dist);
					sum += cur;
					wWeights[i] = cur;
					wIds[i] = ids[pt*(k+1) + i];
				}
			}
			// normalize
			for (i=0; i < last; i++)
				wWeights[i] /= sum;
		}
	}
};

void uDumpPtInfo(int pt);

void calcInnerSprings(TriMesh *mesh) {
	/*
	if (!edgeMatchGF) {
		cout << "edgeMatchGF not initialized; aborting calcInnerSprings" << endl;
		return;
	}

	mesh->calcHBB(16);

	int pt;
	for (pt = 0; pt < mesh->numPts(); pt++) {
		mesh->calcRayIntersection(mesh->getPt(pt) + mesh->getPtNormal(pt) * 0.005, mesh->getPtNormal(pt));
		if (mesh->ptPos > 0) {
			double dist = (mesh->getPt(mesh->ptPos) - mesh->getPt(pt)).length();
			if (dist > 1e-6 && dist < 0.75)
				edgeMatchGF->geodesics[pt].push_back(TMNeigh(mesh->ptPos, pt, dist));
		}
	}*/
}


void UniUI::drawGL() {
	int i;

	if (dispMesh && showTemplate) {
		renderTriMesh(dispMesh, dispMode);
	}

	if (dataSet.numExamples > 0) {
		if (showMesh >= 0) {
			renderTriMesh(dataSet.examples[showMesh].mesh, dispMode);
		}
		else if (uMesh) {
			smDrawInnerSprings();

			// if we're in mesh mode, draw the current mesh
			if (showDef) {
				if (showTex && normalMap) {
					if (!softwareNM) {
						// normal mapping
						glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE_EXT);
						glTexEnvf(GL_TEXTURE_ENV, GL_COMBINE_RGB_EXT, GL_DOT3_RGB_EXT);
						glTexEnvf(GL_TEXTURE_ENV, GL_SOURCE0_RGB_EXT, GL_PRIMARY_COLOR_EXT);
						glTexEnvf(GL_TEXTURE_ENV, GL_OPERAND0_RGB_EXT, GL_SRC_COLOR);
						glTexEnvf(GL_TEXTURE_ENV, GL_SOURCE1_RGB_EXT, GL_TEXTURE);
						glTexEnvf(GL_TEXTURE_ENV, GL_OPERAND1_RGB_EXT, GL_SRC_COLOR);

						glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
						for (i=0; i < 256*256*3; i++) {
							if (i%3 == 0) {
								normalMap[i/3].normalize();
								if (normalMap[i/3].iszero() || normalMap[i/3][2] >= 0)
									normalMap[i/3] = Vec3d(0, 0, -1);
							}
							curTex[i] = max(0, min(255, 
								(int)(normalMap[i/3][i%3] * 128 + 128.5)));
							if (i % 3 == 2)
								curTex[i] = 255 - curTex[i];
						}
						gluBuild2DMipmaps(GL_TEXTURE_2D, 3, NORMAL_MAP_W, NORMAL_MAP_H, 
							GL_RGB, GL_UNSIGNED_BYTE, curTex);
						renderTriMesh(uMesh, dispMode | VM_TEX);
					}
					else {
						glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

						glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
						for (i=0; i < 256*256*3; i++) {
							curTex[i] = max(0, min(255, 
								(int)(normalMap[i/3][i%3] * 128 + 128.5)));
							if (i % 3 == 2)
								curTex[i] = 255 - curTex[i];
						}
						if (!tangentSpace)
							calcTangentSpace(uMesh, &tangentSpace);
						renderNormalMap(normalMap, curTex, texBary, tangentSpace, 
							mainWin->viewer->camera.lightRot, 
							mainWin->viewer->curModelViewMat);
						gluBuild2DMipmaps(GL_TEXTURE_2D, 3, NORMAL_MAP_W, NORMAL_MAP_H, 
							GL_RGB, GL_UNSIGNED_BYTE, curTex);
						renderTriMesh(uMesh, dispMode | VM_TEX);
					}
				}
				else
					renderTriMesh(uMesh, dispMode);
			}

			if (showFrames && edgeMatchGF) {
				for (i=0; i < uMesh->numPts(); i++) {
					glPushMatrix();
					uMesh->getPt(i).glTranslate();
					QuatNorm q(edgeMatchGF->vars[i*edgeMatchGF->varsPerVert+3],
						edgeMatchGF->vars[i*edgeMatchGF->varsPerVert+4],
						edgeMatchGF->vars[i*edgeMatchGF->varsPerVert+5],
						edgeMatchGF->vars[i*edgeMatchGF->varsPerVert+6]);
					Mat4d m = q.toMatrixD();
					double m2[16];
					m.getGLMatrix(m2);
					glMultMatrixd(m2);
					glScaled(0.02, 0.02, 0.02);
					glbAxes();
					glPopMatrix();
				}
			}

/*			if (showGeodesics && edgeMatchGF && edgeMatchGF->geodesics) {
				int j;
				glDisable(GL_LIGHTING);
				for (i=0; i < uMesh->numPts(); i++) {
					for (j=0; j < edgeMatchGF->geodesics[i].size();j++) {
						TMNeigh &n = edgeMatchGF->geodesics[i][j];
						if (n.weight == 1)
							glColor3f(0, 1, 0);
						else
							glColor3f(1, 0, 0);
						glBegin(GL_LINES);
						uMesh->getPt(n.vert).glVertex();
						uMesh->getPt(n.vert2).glVertex();
						glEnd();
					}
				}
				glEnable(GL_LIGHTING);
			}*/

			if (showMocap && visSkel >= 0 && (dataSet.examples[visSkel].lookup != NULL) || (uSolver.markerAssignments != NULL)) {
				// draw spheres at mocap markers
				for (i=0; i < dataSet.examples[visSkel].numPts; i++) {
					Vec3d v = dataSet.examples[visSkel].pts[i];
					glColor3f(1, 0, 0);
					if (!v.iszero())
						glbSphere(v, 0.005);
				}
				if (dataSet.examples[visSkel].lookup != NULL) {
					// draw error lines
					for (i=0; i < dataSet.examples[visSkel].luSize; i++) {
						glDisable(GL_LIGHTING);
						glColor3f(1,1,1);
						int ind = dataSet.examples[visSkel].lookup[i];
						if (ind >= 0 && dataSet.examples[visSkel].ptsConf[ind] > 0) {
							glBegin(GL_LINES);
							dataSet.examples[visSkel].pts[dataSet.examples[visSkel].lookup[i]].glVertex();
							uMesh->getPt(i).glVertex();
							glEnd();
						}
						glEnable(GL_LIGHTING);
					}
				}
			}

			if (showMarkerSprings && visSkel >= 0 && (dataSet.charIndex[dataSet.examples[visSkel].character] == visSkel)) {
				for (i = 0; i < uSolver.mSpringWeights.size(); i++) {
					double w = uSolver.mSpringWeights[i];
					if (w == 0)
						continue;
					char side = skin.skel->transforms.getT(uSolver.mSpringTrans[i])->name[0];
					if (side == 'l')
						glColor3f(0, 1, 0);
					else if (side == 'r')
						glColor3f(1, 0, 0);
					else
						glColor3f(1, 1, 0);
					Vec3d v = uSolver.mSpringVerts2[dataSet.examples[visSkel].character * uSolver. mSpringWeights.size() + i];
					glbSphere(v, 0.005);
				}
			}
		}
		else {
			// otherwise, show the markers
			for (i=0; i < dataSet.examples[visSkel].numPts; i++) {
				if (i < skin.numPts) {
					// reconstructed markers in green
					glColor3f(0, 1, 0);
					if (!skin.curPts[i].iszero())
						glbSphere(skin.curPts[i], 0.005);
				}
				double exConf;
				Vec3d exPt;
				dataSet.examples[visSkel].getPt(i, &exPt, &exConf);
				if (showData && (exConf > 0)) {
					// observed markers in red
					glColor3f(1, 0, 0);
					glbSphere(exPt, 0.005);

					// make lines between the observed & reconstructed markers
					if (i < skin.numPts) {
						glColor3f(1, 1, 1);
						glDisable(GL_LIGHTING);
						glBegin(GL_LINES);
						skin.curPts[i].glVertex();
						exPt.glVertex();
						glEnd();
						glEnable(GL_LIGHTING);
					}
				}
			}
		}
		// render the skeleton
		if (showSkel) {
			skin.skel->drawGL();

			if (showAxes) {
				int tr;
				for (tr=0; tr < skin.skel->transforms.size(); tr++) {
					SkelTransform *trans = skin.skel->transforms.getT(tr);
					if (strstr(trans->className, "Translation") == NULL)
						continue;

					glPushMatrix();
					float m[16];
					trans->globalCoord.mat.getGLMatrixF(m);
					glMultMatrixf(m);
					glScalef(0.2, 0.2, 0.2);
					glbAxes();
					glScalef(-1, -1, -1);
					glbAxes();
					glPopMatrix();
				}
			}

//			cout << skin.skel->transforms.getT("rShoulderCR")->globalCoord.q << endl;
//			cout << skin.skel->transforms.getT("rShoulderYR")->globalCoord.q << endl;
		}
	}


	if (showMarkers && edgeMatchGF) {
		for (i=0; i < edgeMatchGF->srcMarkers->numMarkers; i++) {
			if (!edgeMatchGF->srcMarkers->v(i).iszero()) {
				// reconstructed markers in green
				glColor3f(0, 1, 0);
				glbSphere(edgeMatchGF->srcMarkers->v(i), 0.005);
			}
			if (!edgeMatchGF->markers->markers[i].pos.iszero()) {
				// observed markers in red
				glColor3f(1, 0, 0);
				glbSphere(edgeMatchGF->markers->markers[i].pos, 0.005);

				if (!edgeMatchGF->srcMarkers->v(i).iszero()) {
					glColor3f(1,1,1);
					glDisable(GL_LIGHTING);
					glBegin(GL_LINES);
						edgeMatchGF->srcMarkers->v(i).glVertex();
						edgeMatchGF->markers->markers[i].pos.glVertex();
					glEnd();
					glEnable(GL_LIGHTING);
				}
			}
		}
	}

	for (i=0; i < extraMarkers.size(); i++) {
		glColor3f(0, 1, 1);
		glbSphere(extraMarkers[i], 0.005);
	}

	// tangent space visualization
	if (showTangentSpace && tangentSpace) {
		for (i=0; i < uMesh->numTris() * 3; i++) {
			Mat4d m;
			Vec3d v = uMesh->getPt(uMesh->getTri(i/3, i%3));
			m[0][0] = tangentSpace[i][0][0];
			m[0][1] = tangentSpace[i][0][1];
			m[0][2] = tangentSpace[i][0][2];
			m[1][0] = tangentSpace[i][1][0];
			m[1][1] = tangentSpace[i][1][1];
			m[1][2] = tangentSpace[i][1][2];
			m[2][0] = tangentSpace[i][2][0];
			m[2][1] = tangentSpace[i][2][1];
			m[2][2] = tangentSpace[i][2][2];
			m[0][3] = v[0];
			m[1][3] = v[1];
			m[2][3] = v[2];
			glPushMatrix();
			double m2[16];
			m.getGLMatrix(m2);
			glMultMatrixd(m2);
			glScaled(0.02, 0.02, 0.02);
			glbAxes();
			glPopMatrix();
		}
	}
}

void updateMesh() {
	if (uMesh) {
		int i;
		for (i=0; i < uMesh->numPts(); i++) {
			uMesh->getPt(i) = skin.curPts[i];
//			uMesh->getPtColor(i) = Vec3d(0.8, 0.8, 0.8);
		}
		uMesh->calcNormals();
		if (normalMap && tangentSpace && showTex && softwareNM) {
			calcTangentSpace(uMesh, &tangentSpace);
		}
	}
}

void setPtWeight(int curPt, int curTrans, double newWeight) {
	int inf, curInf = -1, curNumInf = 0;

	if (curPt >= uSolver.numOrigPts) {
		curPt = uSolver.mirrorMap[curPt];
		curTrans = uSolver.mirrorTrans[curTrans];
	}

	for (inf = 0; inf < skin.maxInf; inf++) {
		if (skin.infJoints[curPt * skin.maxInf + inf] == curTrans)
			curInf = inf;
		if (skin.infJoints[curPt * skin.maxInf + inf] >= 0)
			curNumInf++;
	}

	if (curInf < 0 && spAutoExtend && curNumInf < skin.maxInf) {
		skin.infJoints[curPt * skin.maxInf + curNumInf] = curTrans;
		skin.infJoints[uSolver.mirrorMap[curPt] * skin.maxInf + curNumInf] = uSolver.mirrorTrans[curTrans];
		curInf = curNumInf;
		curNumInf++;
	}

	double oldWeight = skin.infWeights[curPt * skin.maxInf + curInf];

	for (inf = 1; inf  < skin.maxInf; inf++) {
		if (inf != curInf) {
			if (oldWeight > 0.99999)
				uSolver.vWeight[curPt * (skin.maxInf - 1) + (inf - 1)] = (1.0 - newWeight) / curNumInf;
			else
				uSolver.vWeight[curPt * (skin.maxInf - 1) + (inf - 1)] *= (1.0 - newWeight) / (1.0 - oldWeight);
		}
		if (uSolver.vWeight[curPt * (skin.maxInf - 1) + (inf - 1)] < 0)
			uSolver.vWeight[curPt * (skin.maxInf - 1) + (inf - 1)] = 0;
		else if (uSolver.vWeight[curPt * (skin.maxInf - 1) + (inf - 1)] > 1)
			uSolver.vWeight[curPt * (skin.maxInf - 1) + (inf - 1)] = 1;
	}
	if (curInf > 0)
		uSolver.vWeight[curPt * (skin.maxInf - 1) + (curInf - 1)] = newWeight;
}

void paintSkin(Vec3d pt, int tri = -1) {
	int i, inf;
	bool hasChanged = false;

	vector<TMNeigh> affectedPoints;

	if (spGeodesic) {
		// find points that are within the radius geodesically (approximately)
		if (tri < 0) {
			// no triangle given; use closest vertex only
			double dist, minDist = 1e10;
			int seedPt = 0;
			for (i = 0; i < uMesh->numPts(); i++) {
				dist = (uMesh->getPt(i) - pt).length();
				if (dist < minDist) {
					minDist = dist;
					seedPt = i;
				}
			}
			findPtGeodesics(uMesh, 1, &seedPt, &minDist, max(spInnerRadius, spOuterRadius), 
				affectedPoints, uMeshNeigh);
		}
		else {
			// seed geodesic calculation with triangle vertices
			int seedPts[3];
			double seedDists[3];
			for (i=0; i < 3; i++) {
				seedPts[i] = uMesh->getTri(tri, i);
				seedDists[i] = (uMesh->getPt(seedPts[i]) - pt).length();
			}
			findPtGeodesics(uMesh, 3, seedPts, seedDists, max(spInnerRadius, spOuterRadius), 
				affectedPoints, uMeshNeigh);
		}
	}
	else {
		// use geometric distance
		for (i = 0; i < uMesh->numPts(); i++){
			// first, check each point to find out if it's within the current paint radius
			double dist = (uMesh->getPt(i) - pt).length();

			if (dist < max(spInnerRadius, spOuterRadius))
				affectedPoints.push_back(TMNeigh(i, 0, dist));
		}
	}

	double *newAffWeights = NULL;
	if (spMode == 4) {
		newAffWeights = new double[affectedPoints.size()];
	}

	for (i=0; i < affectedPoints.size(); i++) {
		int curPt = affectedPoints[i].vert;
		double scale;
		if (affectedPoints[i].dist < spInnerRadius) {
			scale = 1.0;
		}
		else if (affectedPoints[i].dist < spOuterRadius) {
			scale = 1.0 - (affectedPoints[i].dist-spInnerRadius) / (spOuterRadius-spInnerRadius);
		}

		// now figure out which weight to change
		double oldWeight, newWeight;
		int curInf = -1;
		
		for (inf = 0; inf < skin.maxInf; inf++) {
			if (skin.infJoints[curPt * skin.maxInf + inf] == spTransform)
				curInf = inf;
		}

		if (curInf >= 0)
			oldWeight = skin.infWeights[curPt * skin.maxInf + curInf];
		else {
			oldWeight = 0;
		}

		if (spMode == 4) {
			double totalW = uMeshNeigh[curPt].size();
			newWeight = oldWeight * totalW;
			int neigh;
			for (neigh = 0; neigh < uMeshNeigh[curPt].size(); neigh++) {
				for (inf = 0; inf < skin.maxInf; inf++) {
					if (skin.infJoints[uMeshNeigh[curPt][neigh] * skin.maxInf + inf] == spTransform) {
						newWeight += scale * skin.infWeights[uMeshNeigh[curPt][neigh] * skin.maxInf + inf];
						totalW += scale;
						break;
					}
				}
			}
			newAffWeights[i] = newWeight / totalW;
			hasChanged = true;
			continue;
		}

		if (spMode == 1)
			newWeight = scale * spIntensity + (1.0 - scale) * oldWeight;
		else if (spMode == 2)
			newWeight = spIntensity * scale;
		else
			newWeight = 1.0 - scale * (1.0 - spIntensity);

		if ((spMode == 1) || 
			(spMode == 2 && newWeight > oldWeight) || 
			(spMode == 3 && newWeight < oldWeight)) {
			setPtWeight(curPt, spTransform, newWeight);
			hasChanged = true;
		}
	}

	if (spMode == 4) {
		for (i=0; i < affectedPoints.size(); i++) {
			setPtWeight(affectedPoints[i].vert, spTransform, newAffWeights[i]);
		}
		delete []newAffWeights;
	}

	if (hasChanged) {
		char s[80];
		sprintf(s, "%d", spTransform);
		uSolver.updateWeights();
		if (spAutoUpdate)
			updateSkel(NULL);
		uColor(s);
	}
	redrawV();
}

bool UniUI::startDragGL(GLuint *nameBuf, int x, int y) {
//void UniUI::clickGL(GLuint *nameBuf, double x, double y, int button) {
	Vec3d orig, direction;
	y = mainWin->viewer->viewPort[3] - y - 1;
	gluUnProject(GLdouble(x), GLdouble(y), GLdouble(0), mainWin->viewer->modelMatrix, 
		mainWin->viewer->projMatrix, mainWin->viewer->viewPort, &orig[0], &orig[1], &orig[2]);
	gluUnProject(GLdouble(x), GLdouble(y), GLdouble(1), mainWin->viewer->modelMatrix, 
		mainWin->viewer->projMatrix, mainWin->viewer->viewPort, &direction[0], &direction[1], &direction[2]);
	direction = -(direction - orig);
	direction.normalize();

	spDragging = false;

	if (uMesh) {
		if (uMesh->calcRayIntersection(orig, direction)){
			//if it hits the positive side
			if(uMesh->hitNeg){
				double time = uMesh->tNeg;
				Vec3d intersectPt = direction * time + orig;
				double nearestDist = -1;
				Vec3d nearestPt;
				int selectedMeshIndex = -1;

				if (spMode > 0) {
					spDragging = true;
					paintSkin(intersectPt, uMesh->ptNeg);
					return true;
				}
				else {
					//Check each point in the mesh to see if it is closer to the intersection than any others seen so far.
					for(int i = 0; i < uMesh->numPts(); i++){
						double dist = (uMesh->getPt(i) - intersectPt).length();
						//calculate absolute value of dist.
						if(dist < 0){
							dist = -dist;
						}
						if(dist < nearestDist || nearestDist < 0){
							nearestDist = dist;
							nearestPt = uMesh->getPt(i);
							selectedMeshIndex = i;
						}
					}
					if (selectedMeshIndex > -1)
						uDumpPtInfo(selectedMeshIndex);
		//			cout << selectedMeshIndex << endl;
					if (edgeMatchGF) {
						int i;
						for (i=0; i < edgeMatchGF->varsPerVert; i++) {
							cout << edgeMatchGF->vars[selectedMeshIndex*edgeMatchGF->varsPerVert + i] << " ";
						}
						cout << endl;
						QuatNorm q(edgeMatchGF->vars[selectedMeshIndex*edgeMatchGF->varsPerVert + 3],
							edgeMatchGF->vars[selectedMeshIndex*edgeMatchGF->varsPerVert + 4],
							edgeMatchGF->vars[selectedMeshIndex*edgeMatchGF->varsPerVert + 5],
							edgeMatchGF->vars[selectedMeshIndex*edgeMatchGF->varsPerVert + 6]);
						cout << q.toMatrixD() << endl;
					}
				}
			}
		}
	}
	else {
		int i;
		for (i=0; i < skin.numPts; i++) {
			if (skin.curPts[i].iszero())
				continue;

			if (sphereIntersection(orig, direction, skin.curPts[i], 0.005)) {
				cout << "vertex " << i << ": " << skin.curPts[i] << " (" << markerNames[i] << ")" << endl;
			}
		}

		if (visSkel >= 0) {
			for (i=0; i < dataSet.examples[visSkel].numPts; i++) {
				Vec3d v;
				dataSet.examples[visSkel].getPt(i, &v);
				if (v.iszero())
					continue;
				
				if (sphereIntersection(orig, direction, v, 0.005)) {
					cout << "marker " << i << ": " << v << " (" << markerNames[i] << ")" << endl;
				}
			}
		}
	}
	return false;
}

void UniUI::dragGL(int x, int y) {
	if (spDragging) {
		Vec3d orig, direction;
		y = mainWin->viewer->viewPort[3] - y - 1;
		gluUnProject(GLdouble(x), GLdouble(y), GLdouble(0), mainWin->viewer->modelMatrix, 
			mainWin->viewer->projMatrix, mainWin->viewer->viewPort, &orig[0], &orig[1], &orig[2]);
		gluUnProject(GLdouble(x), GLdouble(y), GLdouble(1), mainWin->viewer->modelMatrix, 
			mainWin->viewer->projMatrix, mainWin->viewer->viewPort, &direction[0], &direction[1], &direction[2]);
		direction = -(direction - orig);
		direction.normalize();

		if (uMesh->calcRayIntersection(orig, direction)){
			//if it hits the positive side
			if(uMesh->hitNeg){
				double time = uMesh->tNeg;
				Vec3d intersectPt = direction * time + orig;
				double nearestDist = -1;
				Vec3d nearestPt;
				int selectedMeshIndex = -1;

				paintSkin(intersectPt);
			}
		}
	}
}

void updateSkel(PropUI *) {
	// update the mesh based on the current skeleton (callback for uniUI)
	int i;

	skin.skel->updateCoords();
	skin.updateMats();
//	cout << "shoulder: " << skin.skel->transforms.getT("lShoulderCR")->curCoord.q << endl;
//	cout << "lKneeA: " << skin.skel->transforms.getT("lKneeA")->curCoord.q << endl;
//	cout << "lClavicle: " << skin.skel->transforms.getT("lClavicleCQ")->curCoord.q << endl;
	if (showDress) {
		for (i=0; i < skin.skel->transforms.size(); i++) {
			skin.curMats[i] = Mat4d();
			skin.curMats[i][0][3] = skin.dressJoints[i][0];
			skin.curMats[i][1][3] = skin.dressJoints[i][1];
			skin.curMats[i][2][3] = skin.dressJoints[i][2];
		}
/*		double w[KNN_MAX_SAMPLES];
		Vec3d v[KNN_MAX_SAMPLES];
		knnQuatInterp(skin.skel->transforms.getT("lShoulderQ")->curCoord.q, v, w);
		for (i=0; i < KNN_MAX_SAMPLES; i++) {
			cout << w[i] << " ";
		}
		cout << endl;*/
	}
	skin.updatePts();

	// copy points into the TriMesh
	updateMesh();
	redrawV();

	updatePDDW();
}

bool loadMapping(char *fname, int &numSmallPts, int &numBigPts, int expectedPts, int *&mapPts, float *&mapWeights) {
	int i, j;
	FILE *f2;

	if (!openFile(&f2, fname, "rb", "mapping file"))
		return false;

	// load mapping header
	fread(&numSmallPts, sizeof(int), 1, f2);
	fread(&numBigPts, sizeof(int), 1, f2);
	if (numBigPts != expectedPts) {
		cout << "point size mismatch: mapping is for " << numBigPts << "; expected " << expectedPts << endl;
		fclose(f2);
		return false;
	}

	// load mapping data
	mapPts = new int[numBigPts * MAX_MAPPING_SIZE];
	mapWeights = new float[numBigPts * MAX_MAPPING_SIZE];
	for (i=0; i < numBigPts; i++) {
		int numSources;
		fread(&numSources, sizeof(int), 1, f2);

		for (j=0; j < numSources; j++) {
			fread(mapPts + MAX_MAPPING_SIZE*i + j, sizeof(int), 1, f2);
			fread(mapWeights + MAX_MAPPING_SIZE*i + j, sizeof(float), 1, f2);
		}
		for (; j < MAX_MAPPING_SIZE; j++) {
			mapPts[MAX_MAPPING_SIZE*i+j] = -1;
			mapWeights[MAX_MAPPING_SIZE*i+j] = 0;
		}
	}

	fclose(f2);
	return true;
}

bool loadPoints(const char *datFName) {
	int i;
	FILE *f;
	if (openFile(&f, datFName, "rb", "displacements")) {
		if (datFName[strlen(datFName)-1] == 'x') {
			// load example file (matched mesh)
			int numPoints;
			fread(&numPoints, sizeof(int), 1, f);
			if (numPoints != uMesh->numPts()) {
				cout << "point size mismatch: expected " << uMesh->numPts() << ", found " <<
					numPoints << endl;
				return false;
			}
			for (i=0; i < numPoints; i++) {
				Vec3d &v = uMesh->getPt(i);
				fread(v.n, sizeof(Vec3d), 1, f);
			}
		}
		else {
			char version;
			fread(&version, sizeof(char), 1, f);

			if (version != '0') {
				cout << "unsupported version" << endl;
				return false;
			}

			int fileSize;
			fread(&fileSize, sizeof(int), 1, f);
			if (fileSize != skin.numPts * 3) {
				cout << "wrong filesize: " << fileSize << endl;
				return false;
			}

			for (i=0; i < uMesh->numPts(); i++) {
				Vec3d v;
				fread(v.n, sizeof(double), 3, f);
				uMesh->getPt(i) = v;
			}
		}

		fclose(f);
		uMesh->calcNormals();
		return true;
	}
	else
		return false;
}

void loadMarkerSprings(char *fname) {
	ifstream msIn;
	int i;
	
	if (strlen(fname) > 1 && fname[0] != '-' && openIFStream(&msIn, fname, "marker springs")) {
		int numSprings;
		msIn >> numSprings;
		for (i=0; i < numSprings; i++) {
			char tName[80];
			int p;
			double d;
			msIn >> tName >> d;
			p = skin.skel->transforms.lookupName(tName);
			if (p < 0)
				p = atoi(tName);
			uSolver.mSpringTrans.push_back(p);
			uSolver.mSpringWeights.push_back(d);
		}
		msIn >> uSolver.mSpringNumEx;
		int ex;
		for (ex = 0; ex < uSolver.mSpringNumEx; ex++) {
			for (i=0; i < numSprings; i++) {
				Vec3d v;
				msIn >> v;
				uSolver.mSpringVerts2.push_back(v);
			}
		}
	}
}

void uLoadExampleSet(const char *params) {
	ifstream in;
	int i, j, numComp, numPDComp, ex;
	char fname[80], s[80], mirrorFN[80], mSpringFN[80];
	int charID, lastID;

	params = extractString(params, fname, 80);

	if (!openIFStream(&in, fname, "example set"))
		return;
	uMesh = NULL;
	markerNames = NULL;

	// load skin data
	in >> fname;
	Skeleton *skel = Skeleton::load(fname);
	in >> i >> j;	// numPts, numInf
	if (i == 0) {
		// if i is zero, we're in mesh mode
		uMesh = new TriMesh();
		in >> fname;
		uMesh->loadFile(fname);
		uMesh->calcNormals();
		i = uMesh->numPts();
		uMeshNeigh = findTMNeighbors(uMesh);

		in >> mirrorFN;
	}
	skin.init(skel, i, j);
	if (!uMesh) {
		// if we're in marker mode, then load the initialization weights
		cout << "loading initialization weights" << endl;
		for (i=0; i < skin.numPts; i++) {
			in >> skin.ptNames[i];
			for (j=0; j < skin.maxInf; j++) {
				in >> s;
				if (s[0] == '-') {
					skin.infJoints[i*skin.maxInf + j] = -1;
				}
				else {
					skin.infJoints[i*skin.maxInf + j] = skin.skel->transforms.lookupName(s);
					if (skin.infJoints[i*skin.maxInf + j] < 0) {
						cout << "WARNING -- unknown joint: " << s << endl;
					}
				}
				in >> skin.infWeights[i*skin.maxInf + j];
			}
		}

		in >> mirrorFN;
		in >> mSpringFN;
		loadMarkerSprings(mSpringFN);
	}
	else {
		in >> mSpringFN;
		loadMarkerSprings(mSpringFN);
	}

	// load translation initialization
	in >> skin.numTransInit;
	skin.tiFrames = new char80[skin.numTransInit];
	skin.tiMarkers = new int[skin.numTransInit*2];
	for (i=0; i < skin.numTransInit; i++) {
		in >> skin.tiFrames[i] >> skin.tiMarkers[i*2];
		if (skin.tiMarkers[i*2] < 0) {
			in >> skin.tiMarkers[i*2] >> skin.tiMarkers[i*2+1];
		}
		else {
			skin.tiMarkers[i*2+1] = -1;
		}
	}

	// load point data
	in >> i;
	if (i == -1) {
		// mocap data
		int ofs;
		int step, maxFrame;
		in >> numComp >> numPDComp >> ofs >> step >> maxFrame;
		in >> fname;
		MarkerSet mrefs;
		mrefs.load(fname);
		
		in >> fname;
		if (mocap.loadC3D(fname)) {
			cout << "file has " << mocap.numFrames << " frames " << 
				" and " << mocap.numMarkers << " markers" << endl;
			cout << "using " << (maxFrame-ofs)/step << " frames" << endl;
			dataSet.init((maxFrame-ofs)/step, 1);
			dataSet.charIndex[0] = 0;
			for (i=0;  i < mocap.numFrames * mocap.numMarkers; i++) {
				Vec3d v = mocap.data[i];
				mocap.data[i] = Vec3d(v[2], v[0], v[1] - 1);
			}
			for (i=0; i < (maxFrame-ofs)/step; i++) {
				dataSet.examples[i].character = 0;
				dataSet.examples[i].numPts = mocap.numMarkers;
				dataSet.examples[i].ptsConf = new double[mocap.numMarkers];
				dataSet.examples[i].pts = &mocap.data[(ofs+i*step)*mocap.numMarkers];
				for (j=0; j < mocap.numMarkers; j++) {
					if (dataSet.examples[i].pts[j] == Vec3d(0,0,-1))
						dataSet.examples[i].ptsConf[j] = 0;
					else
						dataSet.examples[i].ptsConf[j] = 1;
				}
			}
			uSolver.markerAssignments = new double[mocap.numMarkers * uMesh->numPts()];
			uSolver.numMarkers = mocap.numMarkers;
			for (i=0; i < mocap.numMarkers * uMesh->numPts(); i++)
				uSolver.markerAssignments[i] = 0; //1.0 / uMesh->numPts();
			for (i=0; i < mocap.numMarkers; i++) {
				uSolver.markerAssignments[i * uMesh->numPts() + mrefs.markers[i].baryVerts[0]] = 1;
			}
//			dataSet.initLookup(mrefs);
		}
	}
	else {
		in >> j >> numComp >> numPDComp;
		dataSet.init(i, j);
		cout << dataSet.numExamples << " examples; " << dataSet.numCharacters << " characters" << endl;

		MarkerSet m;
		lastID = -1;
		for (ex=0; ex < dataSet.numExamples; ex++) {
			in >> charID;
			in >> fname;

			int strInd = (int)strlen(fname) - 1;
			while (strInd > 1 && fname[strInd] != '.') {
				strInd--;
			}
			dataSet.examples[ex].charName[0] = fname[strInd - 5];
			dataSet.examples[ex].charName[1] = fname[strInd - 4];
			dataSet.examples[ex].charName[2] = fname[strInd - 3];
			dataSet.examples[ex].charName[3] = fname[strInd - 2];
			dataSet.examples[ex].charName[4] = 0;
			dataSet.examples[ex].poseCh = fname[strInd - 1];

			if (dataSet.examples[ex].poseCh == 'd') {
				dataSet.examples[ex].minPose = true;
			}

			strcpy(dataSet.examples[ex].fname, fname);
			dataSet.examples[ex].character = charID;
			if (lastID != charID) {
				dataSet.charIndex[charID] = ex;
				lastID = charID;
			}

			if (fname[0] == '-') {
				// build example from initial mesh
				dataSet.examples[ex].init(uMesh->numPts());
				for (j=0; j < uMesh->numPts(); j++) {
					dataSet.examples[ex].ptsConf[j] = 1; // 10
					dataSet.examples[ex].pts[j] = uMesh->getPt(j);
				}
			}
			else if (fname[0] == '*') {
				// null example
				dataSet.examples[ex].init(uMesh->numPts());
				for (j=0; j < uMesh->numPts(); j++) {
					dataSet.examples[ex].ptsConf[j] = 0;
					dataSet.examples[ex].pts[j] = Vec3d();
				}
			}
			else if (fname[strlen(fname)-1] == 'x') {
				// load example file (matched mesh)
				dataSet.examples[ex].load(fname);

				// set to full confidence on 'a' poses
				if (dataSet.examples[ex].poseCh == 'a')
					for (j=0; j < skin.numPts; j++)
						dataSet.examples[ex].ptsConf[j] = 1;
			}
			else if (fname[strlen(fname)-1] == 'y') {
				// load ply file
				dataSet.examples[ex].init(fname);
			}
			else if (fname[strlen(fname)-1] == 't') {
				// load dat (points) file
				dataSet.examples[ex].init(skin.numPts);

				FILE *f;
				if (openFile(&f, fname, "rb", "displacements")) {
					char version;
					fread(&version, sizeof(char), 1, f);

					if (version != '0') {
						cout << "unsupported version" << endl;
						return;
					}

					int fileSize;
					fread(&fileSize, sizeof(int), 1, f);
					if (fileSize != skin.numPts * 3) {
						cout << "wrong filesize: " << fileSize << endl;
						return;
					}

					fread(dataSet.examples[ex].pts, sizeof(double), fileSize, f);
					fclose(f);
					for (j=0; j < skin.numPts; j++)
						dataSet.examples[ex].ptsConf[j] = 10;
				}
			}
			else {
				// load marker file
				if (fname[strlen(fname)-1] == 'd')
					m.loadFromLandmarks(fname);
				else
					m.loadFromMkr(fname);

				// save marker names for later
				if (!markerNames) {
					markerNames = new char80[m.numMarkers];
					for (j=0; j < m.numMarkers; j++)
						strcpy(markerNames[j], m.markers[j].name);
				}

				// hack: eliminate right dactylion on 'C' scans (hand is closed)
				if (dataSet.examples[ex].poseCh == 'c')
					m.markers[38].pos = Vec3d();

				dataSet.examples[ex].init(m.numMarkers);
				for (j=0; j < m.numMarkers; j++) {
					if (!m.markers[j].pos.iszero()) {
						dataSet.examples[ex].ptsConf[j] = 1;
						dataSet.examples[ex].pts[j] = m.markers[j].pos;
					}
					else {
						dataSet.examples[ex].ptsConf[j] = 0;
						dataSet.examples[ex].pts[j] = Vec3d();
					}
				}
			}
		}
	}

	in.close();

	// initialize stuff
	uSolver.init(&dataSet, &skin, uMesh, mirrorFN, numComp, numPDComp);
	uniUI->initSkel(skin.skel);
}

void uShow(int toShow) {
	int i;

	// update shape from the solver
	if (toShow != -1000000)
		visSkel = toShow;
	showMesh = -1;
	showData = true;

	if (visSkel >= 0) {
		int curChar = dataSet.examples[visSkel].character;

		for (i = 0; i < min(uSolver.numComponents, NUM_CUR_COMPS); i++)
			curComps[i] = dataSet.charMu[curChar][i];

		uSolver.updateWeights();
		uSolver.updateSkel(dataSet.charIndex[curChar]);
		skin.updateMats();
		skin.updateJoints();
		uSolver.updateSkel(visSkel);
		uSolver.updatePoints(curChar);
//		updateSkel(NULL);
		skin.updateMats();
		skin.updatePts();

		updateMesh();

/*		if (uSolver.mSpringNumEx == dataSet.numCharacters && dataSet.charIndex[dataSet.examples[visSkel].character] == visSkel) {
			extraMarkers.clear();
			int ind;
			for (ind = 0; ind < uSolver.mSpringWeights.size(); ind++) {
				int joint = uSolver.mSpringTrans[ind];
				SkelTransform *curTrans = skin.skel->transforms.getT(joint);
				double w = uSolver.mSpringWeights[ind];
				Vec3d mPt;
				mPt = uSolver.mSpringVerts2[dataSet.examples[visSkel].character * uSolver.mSpringWeights.size() + ind];
				if (!mPt.iszero()) {
					extraMarkers.push_back(mPt);
				}
			}
		}*/
	}
	else if (uMesh && (-visSkel < dataSet.numExamples)) {
		if (dataSet.examples[-visSkel].mesh) {
			showMesh = -visSkel;
		}
		else {
			for (i=0; i < dataSet.examples[-visSkel].numPts; i++) {
				Vec3d exPt;
				double exConf;
				dataSet.examples[-visSkel].getPt(i, &exPt, &exConf);
				exConf = min(1, exConf);
				uMesh->getPt(i) = exPt;
				uMesh->getPtColor(i) = Vec3d(1, exConf, exConf);
			}
			uMesh->calcNormals();
			if (normalMap && tangentSpace) {
				for (i=0; i < NORMAL_MAP_SIZE; i++) {
					normalMap[i] = dataSet.examples[-visSkel].normals[i];
				}
				calcTangentSpace(uMesh, &tangentSpace);
			}
		}
	}
	redrawV();
}

void uSetPose(const char *params) {
	int ex = 0;
	int i;
	params = extractInt(params, &ex);

	int pp = skin.skel->numPoseDofs * ex;
	for (i=0; i < skin.skel->transforms.size(); i++) {
		SkelTransform *curTrans = skin.skel->transforms.getT(i);
		if (!curTrans->isIntrinsic) {
			curTrans->unloadDofs(uSolver.vPose.n + pp);
			pp += curTrans->numDofs();
		}
	}
}

void uUpdateFromSkel() {
	int i;

//	uSolver.updatePoints(curComps, NUM_CUR_COMPS);
//	updateSkel(NULL);
//	skin.updateMats();
	skin.updatePts();

	updateMesh();
	redrawV();
}

void uShowComps(bool updateSkin) {
	static Skeleton *curSkel = NULL;
	if (!curSkel) {
		curSkel = new Skeleton();
		curSkel->copyFrom(skin.skel);
	}
	else
		curSkel->copyVals(skin.skel);

	showData = false;

	if (updateSkin) {
		uSolver.updateWeights();
		uSolver.updateSkel(curComps, NUM_CUR_COMPS, -1);
		uSolver.updatePoints(curComps, NUM_CUR_COMPS);
		skin.updateMats();
		skin.updateJoints();
	}

	skin.skel->copyVals(curSkel, COPY_POSE);
	skin.skel->updateCoords();
	skin.updateMats();
//	uSolver.updateSkel(curComps, NUM_CUR_COMPS);

	skin.updatePts();

	updateMesh();
	redrawV();
}

void uShow(const char *params) {
	int i = -1;
	params = extractInt(params, &i);
	uShow(i);
}

void uColor(const char *params) {
	// change the colors: -1 = light gray; 0 = all skinning; >0 = single-joint skinning
	int i, j;
	int cMode = -1;
	params = extractInt(params, &cMode);
	
	if (cMode < 0) {
		for (i=0; i < uMesh->numPts(); i++)
			uMesh->getPtColor(i) = Vec3d(1, 1, 1);
	}
	else if (cMode == 0) {
		for (i=0; i < uMesh->numPts(); i++) {
			Vec3d color;
			for (j=0; j < skin.maxInf; j++) 
				if (skin.infJoints[i*skin.maxInf + j] >= 0)
					color += skin.infWeights[i*skin.maxInf + j] *
						skin.skel->transforms.getT(skin.infJoints[i*skin.maxInf + j])->color;
			uMesh->getPtColor(i) = color;
		}
	}
	else {
		for (i=0; i < uMesh->numPts(); i++) {
			Vec3d color = Vec3d(0.5, 0.5, 0.5);
			for (j=0; j < skin.maxInf; j++) 
				if (skin.infJoints[i*skin.maxInf + j] == cMode) {
					color = Vec3d(1, 1.0 - skin.infWeights[i*skin.maxInf + j], 1.0 - skin.infWeights[i*skin.maxInf + j]);
				}
			uMesh->getPtColor(i) = color;
		}
	}
	redrawV();
}

void uColorPdd(const char *params) {
	int i, j;
	int cMode = -1, ind = 0;
	params = extractInt(params, &cMode);
	params = extractInt(params, &ind);
	
	if (cMode >= 0) {
		for (i=0; i < uMesh->numPts(); i++) {
			Vec3d color = Vec3d(0.3, 0.3, 0.3);
			for (j=0; j < skin.pddPtIndex[i].size(); j += 2) {
				int tr = skin.pddPtIndex[i][j];
				if (tr == cMode) {
					//RBF *rbf = skin.pdds[skin.pddPtIndex[i][j]]->rbf;
					int ofs = skin.pddPtIndex[i][j+1];
					color = skin.pddPtKeys[ofs + ind];

					int comp;
					for (comp=0; comp < 3; comp++) {
						color[comp] = min(1.0, max(0.0, (color[comp]/0.1) + 0.5));
					}
				}
			}
			uMesh->getPtColor(i) = color;
		}
	}
	else {
		// test neighbor table
		for (i=0; i < uMesh->numPts(); i++) {
			uMesh->getPtColor(i) = Vec3d();
		}
		for (i=0; i < uSolver.pddNeighTable.size(); i++) {
			NeighborRelation &neigh = uSolver.pddNeighTable[i];
			
			Vec3d color0 = uMesh->getPtColor(neigh.v0);
			if (neigh.ind1 < 0) {
				color0[1] = 1;
			}
			else {
				Vec3d color1 = uMesh->getPtColor(neigh.v1);
				color0[0] += 0.1;
				color1[0] += 0.1;
				uMesh->getPtColor(neigh.v1) = color1;
			}
			uMesh->getPtColor(neigh.v0) = color0;
		}
	}
	redrawV();
}

void uColorPdn(const char *params) {
	int i, j;
	int cMode = -1, ind = 0;
	params = extractInt(params, &cMode);

	for (i=0; i < NORMAL_MAP_SIZE; i++) {
		vector<int> &vi = skin.pddNMIndex[i];

		for (j = 0; j < vi.size(); j+= 2) {
			if (skin.pddNMIndex[i][j] == cMode)
				normalMap[i] = Vec3d();
		}
	}
}

void uColorMarker(const char *params) {
	int i, j;
	int mkr = 0;
	params = extractInt(params, &mkr);
	
	for (i=0; i < uMesh->numPts(); i++) {
		double d = uSolver.markerAssignments[mkr * uMesh->numPts() + i];
		uMesh->getPtColor(i) = hotCold(d, 0, 1);
	}
	redrawV();
}

void uStartSolver(const char *params) {
	int maxIter = 10000;
	char mode[80];
	mode[0] = 0;
	params = extractString(params, mode, 80);
	params = extractInt(params, &maxIter);

	if (strlen(mode) < 1 || mode[0] == '?') {
		cout << "uStartSolver <mode> <max iterations>" << endl;
		cout << "mode options:" << endl;
		cout << "  d  dress" << endl;
		cout << "  b  bones" << endl;
		cout << "  q  poses" << endl;
		cout << "  w  skinning weights" << endl;
		cout << "  p  pose dependent deformations" << endl;
		cout << "  x  PCA weights" << endl;
		cout << "  n  normal maps" << endl;
		cout << "  m  pose dependent normal maps" << endl;
		cout << "  g  (modify global pose only)" << endl;
		cout << "  z  (lock pose zero)" << endl;
		cout << "  i  (ignore point error)" << endl;
		return;
	}

	if (lbfgs) {
		cout << "solver is already running!" << endl;
		return;
	}

	// make sure neighbors are initialized
	if (uSolver.neighTable.size() == 0 && uMesh) {
		uSolver.buildNeighborTable(uMesh);
		uSolver.initWeightStencil();
	}

	uSolver.optDress = (strchr(mode, 'd') != NULL);
	uSolver.optInt = (strchr(mode, 'b') != NULL);
	uSolver.optPose = (strchr(mode, 'q') != NULL);
	uSolver.optWeight = (strchr(mode, 'w') != NULL);
	uSolver.optPDD = (strchr(mode, 'p') != NULL);
	uSolver.optX = (strchr(mode, 'x') != NULL);
	uSolver.globalPoseOnly = (strchr(mode, 'g') != NULL);
	uSolver.lockPoseZero = (strchr(mode, 'z') != NULL);
	uSolver.optNM = (strchr(mode, 'n') != NULL);
	uSolver.optNMP = (strchr(mode, 'm') != NULL);
	uSolver.ignorePoints = (strchr(mode, 'i') != NULL);
	uSolver.initVars(true);
	if (maxIter == -1) {
		int minV = 0;
		int maxV = 10;
		params = extractInt(params, &minV);
		params = extractInt(params, &maxV);
		runTest(&uSolver, uSolver.curVars, minV, maxV);
		uSolver.fromVars(uSolver.curVars);
		return;
	}

	cout << "initializing solver" << endl;
	lbfgs = new LBFGSSolver(&uSolver);
	lbfgs->solve(1e+3, 1e-5, uSolver.curVars, maxIter);
	cout << "final error: " << uSolver.lastErr << endl;
	uSolver.fromVars(uSolver.curVars);

	delete lbfgs;
	lbfgs = NULL;
}

void uStopSolver(const char *params) {
	if (lbfgs)
		lbfgs->stopNow = true;
	else
		cout << "solver is already stopped." << endl;
}

void uRenormalizeX(const char *params) {
	uSolver.renormalizeX();
}

void uAutoWeight(const char *params) {
	uSolver.autoWeight();
}

void uSetMaxSolveEx(const char *params) {
	params = extractInt(params, &uSolver.maxSolveEx);
}

void uSetSingleComponent(const char *params) {
	params = extractInt(params, &uSolver.singleComp);
}

void uNormalizeDofs(const char *params) {
	int i, j;
	for (i=0; i < dataSet.numExamples; i++) {
		uSolver.updateSkel(i);
		for (j=0; j < skin.skel->transforms.size(); j++)
			skin.skel->transforms.getT(j)->normalize();

		// copy back to vPose
		int tr, pp = skin.skel->numPoseDofs * i;
		for (tr=0; tr < skin.skel->transforms.size(); tr++) {
			SkelTransform *curTrans = skin.skel->transforms.getT(tr);
			if (!curTrans->isIntrinsic) {
				curTrans->unloadDofs(uSolver.vPose.n + pp);
				pp += curTrans->numDofs();
			}
		}
	}
}

void uNormalizeNM(const char *params) {
	int i;
	for (i=0; i < NORMAL_MAP_SIZE / 2; i++)
		uSolver.vNM[i].normalize();
}

void uWeightsFromTex(const char *params) {
	char fname[80], outFName[80], s[80];
	ifstream in;
	params = extractString(params, fname, 80);
	params = extractString(params, outFName, 80);
	if (!openIFStream(&in, fname, "texture index"))
		return;
	
	int maxWeights = 0;
	int i, j, ind;
	int numTex;
	int numT = skin.skel->transforms.size();
	int numP = skin.numPts;
	double *weights = new double[numT * numP];
	memset(weights, 0, sizeof(double)*numT * numP);

	in >> numTex;
	for (i=0; i < numTex; i++) {
		in >> s;
		ind = skin.skel->transforms.lookupName(s);
		if (ind < 0) {
			cout << "unknown transform: " << s << endl;
			return;
		}
		in >> fname;
		cout << "loading texture for index " << ind << endl;

		FILE *f;
		if (!openFile(&f, fname, "rb", "texture file"))
			return;

		char version;
		fread(&version, sizeof(char), 1, f);
		if (version != '0') {
			cout << "unknown version; treating as raw data" << endl;
			fseek(f, 0, SEEK_SET);
		}

		int j;
		fread(&j, sizeof(int), 1, f);

		if (j != numP * 3) {
			cout << "point size mismatch; expected " << numP*3 << "; read " << j << endl;
		}
		else {
			for (j=0; j < numP; j++) {
				Vec3d v;
				fread(&v, sizeof(Vec3d), 1, f);
				if (v[2] == 0.8)
					weights[j*numT + ind] = 0.01;
				else if (v[2] > 0.99)
					weights[j*numT + ind] = 0;
				else
					weights[j*numT + ind] = 1.0 - v[2];
			}
		}
		fclose(f);
	}
	in.close();

	FILE *f;
	if (!openFile(&f, outFName, "wb", "influence file"))
		return;

	for (i=0; i < numP; i++) {
		int numWeights = 0;
		for (j=0; j < numT; j++) {
			if (weights[i*numT + j] == 0)
				weights[i*numT + j] = -1;
			else
				numWeights++;
		}
		if (numWeights > maxWeights) maxWeights = numWeights;
		fwrite(weights + i*numT, sizeof(double), numT, f);
	}
	fclose(f);
	cout << "maximum # of weights: " << maxWeights << endl;

	delete []weights;
}

void uColorFromPly(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);

	TriMesh mesh;
	if (!mesh.loadFile(fname))
		return;
	mesh.calcHBB(16);

	int i;
	for (i=0; i < uMesh->numPts(); i++) {
		mesh.calcClosestPoint(uMesh->getPt(i), 1.0);

		Vec3d c;

	}
}

void uSaveExampleMat(const char *params) {
	char fname[80];
	ofstream out;
	params = extractString(params, fname, 80);
	if (!openOFStream(&out, fname, "example matrix"))
		return;

	int ch, pt, comp;
	for (pt=0; pt < dataSet.examples[0].numPts; pt++) {
		for (comp = 0; comp < 3; comp++) {
			for (ch=0; ch < dataSet.numCharacters; ch++) {
				out << dataSet.examples[dataSet.charIndex[ch]].pts[pt][comp] << " ";
			}
			out << endl;
		}
	}

	out.close();
}

void uSaveMarkerSprings(const char *params) {
	int i, ex;
	char fname[80];
	ofstream out;
	params = extractString(params, fname, 80);
	if (!openOFStream(&out, fname, "marker springs"))
		return;
	
	out << (int)uSolver.mSpringTrans.size() << endl;
	for (i=0; i < (int)uSolver.mSpringTrans.size(); i++)
		out << uSolver.mSpringTrans[i] << " " << uSolver.mSpringWeights[i] << endl;
	for (ex = 0; ex < dataSet.numExamples; ex++) {
		for (i=0; i < uSolver.mSpringVerts.size(); i++) {
			Vec3d v;
			dataSet.examples[ex].getPt(uSolver.mSpringVerts[i], &v);
			out << v << endl;
		}
	}

	out.close();
}

void uMakeMarkerSprings(const char *params) {
	int i, ch;
	char fname[80], fname2[80];
	char cmd[256];
	bool saveIntrins = false;
	ofstream out, out2;
	fname2[0] = 0;
	params = extractString(params, fname, 80);
	params = extractString(params, fname2, 80);
	if (!openOFStream(&out, fname, "marker springs"))
		return;
	if (strlen(fname2) > 1) {
		if (!openOFStream(&out2, fname2, "intrinsics matrix"))
			return;
		saveIntrins = true;
	}
	
	out << skin.skel->transforms.size() << endl;
	for (i=0; i < skin.skel->transforms.size(); i++)
		out << skin.skel->transforms.getT(i)->name << " 0" << endl;
	out << dataSet.numCharacters << endl;
	for (ch = 0; ch < dataSet.numCharacters; ch++) {
		int ex = dataSet.charIndex[ch];
		sprintf(cmd, "set ID %s", dataSet.examples[ex].charName);
		processCommand(cmd);
		processCommand("!makeSkel");
		sprintf(cmd, "uPoseToEx %d", ex);
		processCommand(cmd);

		for (i=0; i < skin.skel->transforms.size(); i++) {
			SkelTransform *trans = skin.skel->transforms.getT(i);
			out << trans->globalCoord.v << endl;
			if (trans->isIntrinsic && saveIntrins) {
				int v;
				for (v=0; v < trans->numDofs(); v++) {
					out2 << trans->getDofAddr(v) << " ";
				}
			}
		}
		if (saveIntrins)
			out2 << endl;
	}

	out.close();
	if (saveIntrins)
		out2.close();
}

void uTestRBFs(const char *params) {
/*	char fname[80];
	params = extractString(params, fname, 80);
	ofstream out;

	if (!openOFStream(&out, fname, "rbf test file"))
		return;

	uShow(0);
	int i, j;
	int ind = skin.skel->transforms.lookupName("lElbowA");
	if (ind < 0) {
		cout << "unknown joint, aborting" << endl;
		return;
	}
	SkelEulerRotation *rot = (SkelEulerRotation*)skin.skel->transforms.getT(ind);
	RBF *rbf = skin.transRBFs[ind];
	if (!rbf) {
		cout << "no RBF, aborting" << endl;
		return;
	}
	for (i=0; i < 200; i++) {
		rot->curAngle = -(i / 199.0) * 3.0 + 0.5;
		updateSkel(NULL);
		for (j=0; j < rbf->mN; j++) {
			out << rbf->curWeights[j] << " ";
		}
		out << endl;
	}

	out.close();*/
}

void uLoadTex(const char *params) {
	char fname[80], s[256];
	double x, y, z;
	int ind = 0;
	bool initNM = false;
	params = extractString(params, fname, 80);
	params = extractBool(params, &initNM);
	
	ifstream in;
	if (!openIFStream(&in, fname, "texture OBJ")) {
		return;
	}

	uMesh->gsTriTexCoords(1);

	while (in.good()) {
		in.getline(s, 255);
		if (s[0] == 'v' && s[1] == 't') {
			sscanf(s+3, "%lf %lf %lf", &x, &y, &z);
			if (x < 0) x = 0; if (x >= 1) x = 0.9999;
			if (y < 0) y = 0; if (y >= 1) y = 0.9999;
			if (z < 0) z = 0; if (z >= 1) z = 0.9999;
			x = (x / 256.0) * 252 + (2.0/256);
			y = (y / 256.0) * 252 + (2.0/256);
			z = (z / 256.0) * 252 + (2.0/256);

			uMesh->getTriTexCoords(ind / 3, ind % 3) = Vec3d(x, y, z);
			ind++;
		}
	}
	int numOrig = ind / 3;
	const int map[3] = {0, 2, 1};
	for (; ind < uMesh->numTris() * 3; ind++) {
		uMesh->getTriTexCoords(ind / 3, ind % 3) = uMesh->getTriTexCoords(ind / 3 - numOrig, map[ind % 3]);
		uMesh->getTriTexCoords(ind / 3, ind % 3)[1] = 1.0 - uMesh->getTriTexCoords(ind / 3, ind % 3)[1];
	}

	in.close();

	// make the texture
	glGenTextures(1, &uMesh->textureID);
	glBindTexture(GL_TEXTURE_2D, uMesh->textureID);
	// select modulate to mix texture with color for shading
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	// when texture area is small, bilinear filter the closest mipmap
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST );
	// when texture area is large, bilinear filter the original
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	// the texture wraps over at the edges (repeat)
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	BYTE *data = new BYTE[256*256*3];
	int i, j;
	for (i=0; i < 256; i++) {
		for (j=0; j < 256; j++) {
			if (((i / 8) + (j / 8)) % 2 == 0) {
				data[(i*256+j)*3 + 0] = 255;
			}
			else {
				data[(i*256+j)*3 + 0] = 0;
			}
			data[(i*256+j)*3 + 1] = i;
			data[(i*256+j)*3 + 2] = j;
		}
	}
	gluBuild2DMipmaps(GL_TEXTURE_2D, 3, 256, 256, GL_RGB, GL_UNSIGNED_BYTE, data);
	delete []data;

//	showTex = true;

	cout << "texture loaded" << endl;

	if (initNM) {
		if (!normalMap)
			normalMap = new Vec3d[NORMAL_MAP_SIZE];
		if (!curTex)
			curTex = new unsigned char[NORMAL_MAP_SIZE * 3];
		if (!texBary) {
			texBary = new TexBary[NORMAL_MAP_SIZE];
			calcTexBary(uMesh, texBary);
		}
		for (i=0; i < NORMAL_MAP_SIZE; i++)
			normalMap[i] = Vec3d(0, 0, -1);
		calcTangentSpace(uMesh, &tangentSpace);
	}
}

void uTangentSpaceVis(const char *params) {
	bool enable = true;
	params = extractBool(params, &enable);
	showTangentSpace = enable;

	if (showTangentSpace) {
		calcTangentSpace(uMesh, &tangentSpace);
	}
	redrawV();
}

void uBuildNormalMap(const char *params) {
	int i;
	char fname[80], datFName[80];
	params = extractString(params, datFName, 80);
	params = extractString(params, fname, 80);

	TriMesh *target = new TriMesh();
	if (!target->loadFile(fname))
		return;

	loadPoints(datFName);

	if (!normalMap)
		normalMap = new Vec3d[NORMAL_MAP_SIZE];
	if (!curTex)
		curTex = new unsigned char[NORMAL_MAP_SIZE * 3];
	if (!texBary) {
		texBary = new TexBary[NORMAL_MAP_SIZE];
		calcTexBary(uMesh, texBary);
	}
	calcTangentSpace(uMesh, &tangentSpace);
	buildNormalMap(uMesh, texBary, target, normalMap, tangentSpace);
	delete target;
/*	symmetrifyNormalMap(normalMap);
	diffuseNormalMap(normalMap);

	renderNormalMap(normalMap, curTex, texBary, tangentSpace,
		mainWin->viewer->camera.lightRot, mainWin->viewer->curProjectionMat);
	gluBuild2DMipmaps(GL_TEXTURE_2D, 3, 256, 256, GL_RGB, GL_UNSIGNED_BYTE, curTex);

	// save TGA
	FILE *fptr = fopen("data/lightmap.tga", "wb");
	if (fptr == NULL) {
		cout << "can't open " << fname << endl;
		return;
	}
	// swap R and B
	for (i=0; i < 256*256; i++) {
		curTex[i * 3 + 0] ^= curTex[i * 3 + 2];
		curTex[i * 3 + 2] ^= curTex[i * 3 + 0];
		curTex[i * 3 + 0] ^= curTex[i * 3 + 2];
	}
	// create tga header
	putc(0,fptr);
	putc(0,fptr);
	putc(2,fptr);                         // uncompressed RGB
	putc(0,fptr); putc(0,fptr);
	putc(0,fptr); putc(0,fptr);
	putc(0,fptr);
	putc(0,fptr); putc(0,fptr);           // X origin
	putc(0,fptr); putc(0,fptr);           // y origin
	putc((256 & 0x00FF),fptr);
	putc((256 & 0xFF00) / 256,fptr);
	putc((256 & 0x00FF),fptr);
	putc((256 & 0xFF00) / 256,fptr);
	putc(24,fptr);                        // 24 bit bitmap
	putc(0,fptr);
	// write the data
	fwrite(curTex, 256*256*3*sizeof(char), 1, fptr);
	fclose(fptr);*/
}

void uBuildNormalMapEx(const char *params) {
	int ex = 0;
	char str[1024];
	params = extractInt(params, &ex);

/*	sprintf(str, "%s data/csr/csr%s%c.ply",
		dataSet.examples[ex].fname, // data/csr-4k/csr%s%c.ex dataSet.examples[ex].charName, dataSet.examples[ex].poseCh,
		dataSet.examples[ex].charName, dataSet.examples[ex].poseCh);
*/
	char s[3];
	s[0] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname) - 5];
	if (s[0] == 'o') s[0] = '0';
	s[1] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname) - 4];
	s[2] = 0;
	sprintf(str, "%s data/drago/scans/mesh%d.ply",
		dataSet.examples[ex].fname, atoi(s));
	cout << "building normal map: " << s << ", " << str << endl;

	uBuildNormalMap(str);

	memcpy(dataSet.examples[ex].normals, normalMap, sizeof(Vec3d) * NORMAL_MAP_SIZE);
//	sprintf(str, "data/csr-7k/csr%s%c.ex",
//		dataSet.examples[ex].charName, dataSet.examples[ex].poseCh);
//	cout << "saved " << str << endl;
	dataSet.examples[ex].save(dataSet.examples[ex].fname);
	cout << "saved " << dataSet.examples[ex].fname << endl;

	redrawV();
	uiWait();
}

void uSaveSolve(const char *params) {
	char fname[80];
	FILE *f;
	params = extractString(params, fname, 80);
	if (!openFile(&f, fname, "wb", "solve file"))
		return;
	uSolver.save(f);
	fclose(f);
}

void uLoadSolve(const char *params) {
	char fname[80];
	FILE *f;
	params = extractString(params, fname, 80);
	if (!openFile(&f, fname, "rb", "solve file"))
		return;
	uSolver.load(f);
	fclose(f);
}

void uSavePoses(const char *params) {
	// allowable modes:
	//  0 default
	//  1 intrinsics only
	//  2 single-character
	int mode = 0;
	char fname[80];
	params = extractString(params, fname, 80);
	params = extractInt(params, &mode);
	uSolver.savePoses(fname, mode);
}

void uLoadPoses(const char *params) {
	bool intOnly = false;
	char fname[80];
	params = extractString(params, fname, 80);
	params = extractBool(params, &intOnly);
	uSolver.loadPoses(fname, intOnly);
}

void uSaveWeights(const char *params) {
	int i, j;
	char fname[80];
	FILE *f;
	params = extractString(params, fname, 80);
	if (!openFile(&f, fname, "wb", "weight file"))
		return;
	int numT = skin.skel->transforms.size();
	double *weights = new double[numT];
	for (i=0; i < skin.numPts; i++) {
		for (j=0; j < numT; j++)
			weights[j] = -1;
		for (j=0; j < skin.maxInf; j++) {
			if (skin.infJoints[i*skin.maxInf + j] < 0)
				continue;
			weights[skin.infJoints[i*skin.maxInf + j]] = skin.infWeights[i*skin.maxInf + j];
		}
		fwrite(weights, sizeof(double), numT, f);
	}
	delete []weights;
	fclose(f);
}

void uLoadWeights(const char *params) {
	char fname[80], fname2[80];
	FILE *f, *f2 = NULL;
	int i, j, k, m;

	fname2[0] = 0;
	params = extractString(params, fname, 80);
	params = extractString(params, fname2, 80);
	if (!openFile(&f, fname, "rb", "weight file"))
		return;
	if (fname2[0] != 0 && !openFile(&f2, fname2, "rb", "map file"))
		return;

	int numSmallPts, numBigPts;
	int *mapPts;
	float *mapWeights;
	int numT = skin.skel->transforms.size();
	double *weights = new double[numT], *weights2 = new double[numT];
	int *topInf = new int[skin.maxInf * 2];
	numSmallPts = skin.numPts;

	if (fname2[0] != 0 && !loadMapping(fname2, numSmallPts, numBigPts, skin.numPts, mapPts, mapWeights))
		return;

	memset(skin.infWeights, 0, sizeof(double)* skin.numPts * skin.maxInf);
	const double NULL_WEIGHT = -1000000;

	for (i = 0; i < uSolver.numOrigPts; i++) {
		if (f2) {
			for (j=0; j < numT; j++)
				weights[j] = NULL_WEIGHT;
			for (j=0; j < MAX_MAPPING_SIZE; j++) {
				if (mapPts[i*MAX_MAPPING_SIZE + j] >= 0 && mapWeights[i*MAX_MAPPING_SIZE + j] != 0) {
					fseek(f, sizeof(double)*numT*mapPts[i*MAX_MAPPING_SIZE + j], SEEK_SET);
					fread(weights2, sizeof(double), numT, f);
					for (k=0; k < numT; k++) {
						if (weights2[k] != -1 && weights2[k] != NULL_WEIGHT) {
							if (weights[k] == NULL_WEIGHT)
								weights[k] = 0;
//							if (k == 1)
//								cout << "base" << endl;
							weights[k] += mapWeights[i*MAX_MAPPING_SIZE + j] * weights2[k];
						}
					}
				}
			}
		}
		else {
			fread(weights, sizeof(double), numT, f);
			for (j=0; j < numT; j++)
				// for backwards compatibility
				if (weights[j] == -1)
					weights[j] = NULL_WEIGHT;
		}

		for (j=0; j < skin.maxInf*2; j++)
			topInf[j] = -1;

		for (j=0; j < numT; j++) {
			if (weights[j] != NULL_WEIGHT) {
//				if (j == 1)
//					cout << "baseQ" << endl;
				for (k=0; k < skin.maxInf*2; k++) {
					if (topInf[k] < 0 || weights[topInf[k]] < weights[j]) {
						for (m=skin.maxInf*2-1; m > k; m--)
							topInf[m] = topInf[m-1];
						topInf[k] = j;
						break;
					}
				}
			}
		}

		double sum = 0;
		for (j=0; j < skin.maxInf; j++) {
			if (topInf[j] != -1) {
				SkelTransform *curT = skin.skel->transforms.getT(topInf[j]);
				if (strcmp(curT->className, "SkelTranslation") == 0 && curT->parentPtr) {
					topInf[j] = curT->parentPtr->index;
				}
			}
			skin.infJoints[i*skin.maxInf + j] = topInf[j];
			if (topInf[j] == -1)
				skin.infWeights[i*skin.maxInf + j] = 0;
			else {
				skin.infWeights[i*skin.maxInf + j] = max(weights[topInf[j]], 0.001);
				sum += skin.infWeights[i*skin.maxInf + j];
			}
		}
//		if (topInf[j] >= 0) {
//			cout << "warning: not enough influences" << endl;
//		}
//		if (sum <= 0 || !_finite(sum))
//			cout << "bad times!" << endl;
		for (j=0; j < skin.maxInf; j++)
			skin.infWeights[i*skin.maxInf + j] /= sum;

		if (uSolver.mirrorMap[i] != i) {
			for (j=0; j < skin.maxInf; j++) {
				int trans = skin.infJoints[i*skin.maxInf + j];
				if (trans >= 0)
					trans = uSolver.mirrorTrans[trans];
				skin.infJoints[uSolver.mirrorMap[i]*skin.maxInf + j] = trans;
				skin.infWeights[uSolver.mirrorMap[i]*skin.maxInf + j] =
					skin.infWeights[i*skin.maxInf + j];
			}
		}
	}

	if (f2) {
		delete []mapPts;
		delete []mapWeights;
	}

	delete []topInf;
	delete []weights;
	delete []weights2;
	fclose(f);

	uSolver.weightsToVars();
}

void uSavePdds(const char *params) {
	int i, j;
	char fname[80];
	FILE *f;
	params = extractString(params, fname, 80);
	if (!openFile(&f, fname, "wb", "PDD file"))
		return;

	fwrite(uSolver.vPDD.n, sizeof(double), uSolver.vPDD.size(), f);

	fclose(f);
}

void uLoadPdds(const char *params) {
	int i, j;
	char fname[80];
	fname[0] = 0;
	FILE *f;
	params = extractString(params, fname, 80);

	if (fname[0] == 0) {
		uSolver.vPDD.zeroElements();
		return;
	}

	if (!openFile(&f, fname, "rb", "PDD file"))
		return;

	fread(uSolver.vPDD.n, sizeof(double), uSolver.vPDD.size(), f);

	fclose(f);
}

void uSaveDress(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);
	uSolver.saveDress(fname);
}

void uLoadDress(const char *params) {
	char fname[80], fname2[80];
	params = extractString(params, fname, 80);
	fname2[0] = 0;
	params = extractString(params, fname2, 80);

	uSolver.loadDress(fname, fname2);
}

void uBuildMapping(const char *params) {
	char fname1[80], fname2[80], fnameO[80];
	TriMesh newMesh;
	int i, j, k, numOrigPts;
	TriMesh *origMesh = NULL;
	MarkerSet *origMS = NULL;

	params = extractString(params, fname1, 80);
	params = extractString(params, fname2, 80);
	params = extractString(params, fnameO, 80);

	// load original data
	if (fname1[strlen(fname1)-1] == 'r') {
		origMS = new MarkerSet();
		origMS->loadFromMkr(fname1);
		numOrigPts = origMS->numMarkers;
	}
	else {
		origMesh = new TriMesh();
		if (!origMesh->loadFile(fname1))
			return;
		numOrigPts = origMesh->numPts();
	}
	if (!newMesh.loadFile(fname2))
		return;
	newMesh.calcHBB(16);
	newMesh.closestRestrictNormal = false;

	KNNStruct knn(newMesh.numPts(), 4);

	cout << "initializing knn" << endl;
	for (i=0; i < numOrigPts; i++) {
		Vec3d v;
		if (origMS)
			v = origMS->markers[i].pos;
		else
			v = newMesh.getPt(i);
		if (!newMesh.calcClosestPoint(v, 1.10))
			cout << "WARNING: no closest point for vert " << i << endl;
		else {
			for (j=0; j < 3; j++) {
				int pt = newMesh.closestTri[j];
				if (pt > -1) {
					knn.insertDist(pt, (v-newMesh.getPt(pt)).length2(), i);
				}
			}
		}
	}

	cout << "propagating knn" << endl;
	bool isNew = false;
	for (i=0; i < 100; i++) {
		for (j=0; j < newMesh.numTris(); j++) {
			for (k=0; k < 3; k++) {
				int v0 = newMesh.getTri(j, k);
				int v1 = newMesh.getTri(j, (k+1)%3);
				double edgeLen = (newMesh.getPt(v1) - newMesh.getPt(v0)).length2();
				int index;
				for (index=0; index < knn.k+1; index++) {
					if (knn.dists[v0*(knn.k+1) + index] < DBL_MAX)
						isNew |= knn.insertDist(v1, knn.dists[v0*(knn.k+1) + index] + edgeLen, knn.ids[v0*(knn.k+1) + index]);
				}
			}
		}
		if (!isNew)
			break;
	}
	
	cout << "writing map..." << endl;
	FILE *f;
	if (openFile(&f, fnameO, "wb", "map file")) {
		fwrite(&numOrigPts, sizeof(int), 1, f);
		i = newMesh.numPts();
		fwrite(&i, sizeof(int), 1, f);
		for (i=0; i < newMesh.numPts(); i++) {
			int ids[4];
			double weights[4];
			knn.getWeights(i, ids, weights);
			if (!_finite(weights[0]))
				cout << "problem!" << endl;
			fwrite(ids, sizeof(int), 4, f);
			fwrite(weights, sizeof(double), 4, f);
		}
	}
	fclose(f);

	if (origMS) delete origMS;
	if (origMesh) delete origMesh;
}

void uDressFromEx(const char *params) {
	int ex = 0;
	int pt, inf;

	extractInt(params, &ex);

	uShow(ex);
	for (pt = 0; pt < skin.numPts; pt++) {
		Mat3d m;
		Vec3d base;
		skin.getSkinMat(pt, m, base);

		Vec3d v;
		dataSet.examples[ex].getPt(pt, &v);
		if (v.iszero())
			skin.dressPts[pt] = Vec3d();
		else
			skin.dressPts[pt] = m.inverse() * (v - base);
	}

	uSolver.dressToVars(dataSet.examples[ex].character);
//	uSolver.dressToVars();

	uShow(ex);
}

void saveDisplacements(const char *fname, Vecd &v) {
	FILE *f;
	if (!openFile(&f, fname, "wb", "displacements")) {
		return;
	}

	char version = '0';
	fwrite(&version, sizeof(char), 1, f);
	int size = v.size();
	fwrite(&size, sizeof(int), 1, f);
	fwrite(v.n, sizeof(double), size, f);
	fclose(f);
}

void uSaveSkinDisp(const char *params) {
	int ex, firstEx = 0, lastEx = 1000;
	char fname[80], fmask[80];
	int pt, i, inf;
	Vecd vars(skin.numPts * 12);
	Mat4d *origMats = new Mat4d[skin.numPts];

	params = extractInt(params, &firstEx);
	params = extractInt(params, &lastEx);
	params = extractString(params, fmask);

	for (ex=firstEx; ex < min(lastEx, dataSet.numExamples); ex++) {
		if (dataSet.examples[ex].poseCh == 'a')
			continue;

		uShow(dataSet.charIndex[dataSet.examples[ex].character]);
		for (pt = 0; pt < skin.numPts; pt++) {
			Vec3d ofs;
			origMats[pt] = Mat4d(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
			for (inf=0; inf < skin.maxInf; inf++) {
				int joint = skin.infJoints[pt*skin.maxInf + inf];
				if (joint >= 0) {
					double w = skin.infWeights[pt*skin.maxInf + inf];

					origMats[pt] += w * skin.curMats[joint];
					ofs += w * (Vec3d(skin.curMats[joint][0][3], skin.curMats[joint][1][3],
						skin.curMats[joint][2][3]) - vec4to3(skin.curMats[joint] * 
						vec3to4z(skin.dressJoints[joint])));
				}
			}
			origMats[pt][0][3] = ofs[0];
			origMats[pt][1][3] = ofs[1];
			origMats[pt][2][3] = ofs[2];
		}


		uShow(ex);
		for (pt = 0; pt < skin.numPts; pt++) {
			Mat4d m(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
			Vec3d ofs;

			for (inf=0; inf < skin.maxInf; inf++) {
				int joint = skin.infJoints[pt*skin.maxInf + inf];
				if (joint >= 0) {
					double w = skin.infWeights[pt*skin.maxInf + inf];

					m += w * skin.curMats[joint];
					ofs += w * (Vec3d(skin.curMats[joint][0][3], skin.curMats[joint][1][3],
						skin.curMats[joint][2][3]) - vec4to3(skin.curMats[joint] * 
						vec3to4z(skin.dressJoints[joint])));
				}
			}
			m[0][3] = ofs[0];
			m[1][3] = ofs[1];
			m[2][3] = ofs[2];
			m = m * origMats[pt].inverse();

			for (i=0; i < 12; i++)
				vars[pt*12 + i] = m.n[i];
		}

		sprintf(fname, fmask, dataSet.examples[ex].charName, dataSet.examples[ex].poseCh);
		cout << "saving " << fname << endl;
		saveDisplacements(fname, vars);	
	}

	delete []origMats;
}

void uZeroVars(const char *params) {
	char toZero[80];
	params = extractString(params, toZero, 80);

	if (strchr(toZero, 'q')) {
		cout << "zeroing pose variables" << endl;

		skin.skel->zero();
		int i, pp = 0;
		for (i=0; i < skin.skel->transforms.size(); i++) {
			SkelTransform *curTrans = skin.skel->transforms.getT(i);
			if (!curTrans->isIntrinsic) {
				int ex;
				for (ex=0; ex < dataSet.numExamples; ex++)
					curTrans->unloadDofs(uSolver.vPose.n + pp + skin.skel->numPoseDofs * ex);
				pp += curTrans->numDofs();
			}
		}
	}

	if (strchr(toZero, 'd')) {
		cout << "zeroing dress variables" << endl;
		uSolver.vDress.zeroElements();
	}

	if (strchr(toZero, 'p')) {
		cout << "zeroing PDD variables" << endl;
		uSolver.vPDD.zeroElements();
	}

	if (strchr(toZero, 'x')) {
		cout << "zeroing recon variables" << endl;
		int ch, comp;
		for (ch=0; ch < dataSet.numCharacters; ch++) {
			dataSet.charMu[ch][0] = 1;
			for (comp = 1; comp < uSolver.numComponents; comp++)
				dataSet.charMu[ch][comp] = 0;
		}
	}

	if (strchr(toZero, 'n')) {
		cout << "zeroing normal map variables" << endl;
		int i;
		for (i=0; i < NORMAL_MAP_SIZE / 2; i++) {
			uSolver.vNM[i] = Vec3d(0, 0, -1);
		}
		for (; i < NORMAL_MAP_SIZE / 2 * uSolver.numComponents; i++) {
			uSolver.vNM[i] = Vec3d(0, 0, 0);
		}
	}

	if (strchr(toZero, 'm')) {
		cout << "zeroing PD normal map variables" << endl;
		uSolver.vNMP.zeroElements();
	}

	if (strchr(toZero, '!')) {
		cout << "zeroing components above zero" << endl;
		int i;
		for (i = uSolver.numOrigPts * 3; i < uSolver.vDress.size(); i++)
			uSolver.vDress[i] = 0;
		for (i = uSolver.skin->skel->numIntrinsicDofs; i < uSolver.vInt.size(); i++)
			uSolver.vInt[i] = 0;
		for (i = uSolver.skin->numPddPtKeys * 3; i < uSolver.vPDD.size(); i++)
			uSolver.vPDD[i] = 0;
		for (i = NORMAL_MAP_SIZE / 2; i < NORMAL_MAP_SIZE / 2 * uSolver.numComponents; i++)
			uSolver.vNM[i] = Vec3d();
		for (i = uSolver.skin->numPddNMKeys * 3; i < uSolver.vNMP.size(); i++)
			uSolver.vNMP[i] = 0;
	}
}

QuatNorm determineRotation(Vec3d q0, Vec3d dq0, Vec3d q1, Vec3d dq1, Vec3d q2 = Vec3d(), Vec3d dq2 = Vec3d()) {
	Vec3d v0[4];
	Vec3d v1[4];
	v0[0] = Vec3d();
	v1[0] = Vec3d();
	v0[1] = q0; //v0[1].normalize();
	v1[1] = dq0; //v1[1].normalize();
	v0[2] = q1; //v0[2].normalize();
	v1[2] = dq1; //v1[2].normalize();
	v0[3] = q2;
	v1[3] = dq2;
	Mat3d m = ptsRot(v0, v1, (q2.iszero())? 3 : 4);
//	cout << v0[1] << " -> " << m * v0[1] << " = " << v1[1] << endl;
//	cout << " -> " << (matToQuat(m)).toMatrixD() * q0 << endl;
//	cout << v0[2] << " -> " << m * v0[2] << " = " << v1[2] << endl;

	return matToQuat(m);
}

void uSkelFromSurf(const char *params) {
	char fname[80], datFName[80];
	int i;
	SkelTranslation *tTrans;
	SkelQuatRotation *tQuat;
	SkelEulerRotation *tEuler;
	double len, angle;
	Vec3d right, left, axis;
	QuatNorm r;
	Mat4d parent;
	MarkerSet ms;

	params = extractString(params, fname, 80);
	params = extractString(params, datFName, 80);
	ms.load(fname);

	if (!loadPoints(datFName))
		return;

/*	// set mesh
	ex = dataSet.charIndex[ch];
	for (i=0; i < dataSet.examples[ex].numPts; i++) {
		Vec3d exPt;
		double exConf;
		dataSet.examples[ex].getPt(i, &exPt, &exConf);
		uMesh->getPt(i) = exPt;
	}
	uMesh->calcNormals();*/

	// connect markers to mesh
	for (i=0; i < ms.numMarkers; i++) {
		ms.markers[i].mesh = uMesh;
	}

	// set root
	tTrans = (SkelTranslation*)skin.skel->transforms.getT("baseT");
	Vec3d rootPos = 0.5 * (ms.markers[0].curPos() + ms.markers[1].curPos());
	tTrans->curVal = rootPos;

	// calculate hip offsets and global rotation
	Vec3d rHip = (0.86 * ms.markers[22].curPos() + 0.14 * ms.markers[23].curPos()) - rootPos;
	Vec3d lHip = (0.14 * ms.markers[22].curPos() + 0.86 * ms.markers[23].curPos()) - rootPos;
	len = 0.5 * (rHip.length() + lHip.length());
	rHip.normalize(); rHip *= len;
	lHip.normalize(); lHip *= len;
	Vec3d hipDiff = lHip - rHip;
	Vec3d upDir = 0.5 * (ms.markers[2].curPos() + ms.markers[3].curPos()) - rootPos;
	r = determineRotation(upDir, Vec3d(0, 0, upDir.length()), hipDiff, Vec3d(0, -hipDiff.length(), 0));
/*	angle = atan2(hipDiff[0], -hipDiff[1]);
	r = QuatNorm(0, 0, angle);
	hipDiff = r.toMatrixD() * hipDiff;
	angle = atan2(-hipDiff[2], -hipDiff[1]);
	r = QuatNorm(angle, 0, 0) * r;
	lHip = r.toMatrixD() * lHip;
	cout << "baseq: " << r << endl;
*/	tQuat = (SkelQuatRotation*)skin.skel->transforms.getT("baseQ");
	tTrans = (SkelTranslation*)skin.skel->transforms.getT("lPelvisT");
	tQuat->curQuat = r;
	tQuat->curQuat.w *= -1;
	tTrans->curVal = r.toMatrixD() * lHip;
	skin.skel->updateCoords();

	// calculate knee offsets and hip rotation
	Vec3d rKnee = skin.skel->transforms.getT("rPelvisT")->globalCoord.mat.inverse() * 
		(0.5 * (ms.markers[24].curPos() + ms.markers[25].curPos()));
	Vec3d lKnee = skin.skel->transforms.getT("lPelvisT")->globalCoord.mat.inverse() * 
		(0.5 * (ms.markers[29].curPos() + ms.markers[30].curPos()));
	len = 0.5 * (rKnee.length() + lKnee.length());
	rKnee.normalize(); rKnee *= len;
	lKnee.normalize(); lKnee *= len;
	angle = atan2(rKnee[1], -rKnee[2]);
	r = QuatNorm(angle, 0, 0);
	rKnee = r.toMatrixD() * rKnee;
	angle = atan2(-rKnee[0], -rKnee[2]);
	r = QuatNorm(0, angle, 0) * r;
	tQuat = (SkelQuatRotation*)skin.skel->transforms.getT("rHipQ");
	tQuat->curQuat = r;
	tQuat->curQuat.w *= -1;
	angle = atan2(lKnee[1], -lKnee[2]);
	r = QuatNorm(angle, 0, 0);
	lKnee = r.toMatrixD() * lKnee;
	angle = atan2(-lKnee[0], -lKnee[2]);
	r = QuatNorm(0, angle, 0) * r;
	tQuat = (SkelQuatRotation*)skin.skel->transforms.getT("lHipQ");
	tQuat->curQuat = r;
	tQuat->curQuat.w *= -1;
	tTrans = (SkelTranslation*)skin.skel->transforms.getT("lThighT");
	tTrans->curVal[2] = -len;
	skin.skel->updateCoords();

	// calculate shin offsets and knee rotation
	Vec3d rAnkle = skin.skel->transforms.getT("rThighT")->globalCoord.mat.inverse() * 
		(0.5 * (ms.markers[26].curPos() + ms.markers[27].curPos()));
	Vec3d lAnkle = skin.skel->transforms.getT("lThighT")->globalCoord.mat.inverse() * 
		(0.5 * (ms.markers[31].curPos() + ms.markers[32].curPos()));
	len = 0.5 * (rAnkle.length() + lAnkle.length());
	rAnkle.normalize(); rAnkle *= len;
	lAnkle.normalize(); lAnkle *= len;
	// right
	// knee carrying angle
	angle = atan2(rAnkle[1], -rAnkle[2]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("rKneeCarryingA");
	tEuler->curAngle = angle;
	r = QuatNorm(angle, 0, 0);
	rAnkle = r.toMatrixD() * rAnkle;
	// knee angle
	angle = atan2(-rAnkle[0], -rAnkle[2]);
	rAnkle = QuatNorm(0, angle, 0).toMatrixD() * rAnkle;
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("rKneeA");
	tEuler->curAngle = angle;
	// left
	// knee carrying angle
	angle = atan2(lAnkle[1], -lAnkle[2]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("lKneeCarryingA");
	tEuler->curAngle = angle;
	r = QuatNorm(angle, 0, 0);
	lAnkle = r.toMatrixD() * lAnkle;
	// knee angle
	angle = atan2(-lAnkle[0], -lAnkle[2]);
	lAnkle = QuatNorm(0, angle, 0).toMatrixD() * lAnkle;
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("lKneeA");
	tEuler->curAngle = angle;
	// shin length
	tTrans = (SkelTranslation*)skin.skel->transforms.getT("lShinT");
	rAnkle[1] *= -1;
	tTrans->curVal = Vec3d(0, 0, 0.5 * (lAnkle[2] + rAnkle[2]));
	skin.skel->updateCoords();

	// calculate foot offsets and ankle rotation
	Vec3d rFoot = skin.skel->transforms.getT("rShinT")->globalCoord.mat.inverse() *
		ms.markers[28].curPos();
	Vec3d lFoot = skin.skel->transforms.getT("lShinT")->globalCoord.mat.inverse() * 
		ms.markers[33].curPos();
	rFoot[1] *= -1;
	tTrans = (SkelTranslation*)skin.skel->transforms.getT("lFootT");
	tTrans->curVal = 0.5 * (lFoot + rFoot);
	skin.skel->updateCoords();

	// calculate abdomen offset and waist rotation
	parent = skin.skel->transforms.getT("baseQ")->globalCoord.mat.inverse();
	Vec3d abdomen = parent * (0.5 * (ms.markers[2].curPos() + ms.markers[3].curPos()));
//	cout << "ab: " << abdomen << endl
//		<< "across: " << (parent * (ms.markers[3].curPos() - ms.markers[2].curPos())) << endl;
	parent[0][3] = 0; parent[1][3] = 0; parent[2][3] = 0;
	axis = parent * (ms.markers[3].curPos() - ms.markers[2].curPos());
	r = determineRotation(abdomen, Vec3d(0, 0, abdomen.length()), 
		axis, Vec3d(0, -axis.length(), 0),
		-axis, Vec3d(0, axis.length(), 0));
	tQuat = (SkelQuatRotation*)skin.skel->transforms.getT("waistQ");
	r.w *= -1;
//	cout << "ab: " << r.toMatrixD() * abdomen << endl
//		<< "across: " << r.toMatrixD() * (parent * (ms.markers[3].curPos() - ms.markers[2].curPos())) << endl;
	tQuat->curQuat = r;
	tTrans = (SkelTranslation*)skin.skel->transforms.getT("abdomenT");
	tTrans->curVal[0] = 0;
	tTrans->curVal[1] = 0;
	tTrans->curVal[2] = abdomen.length();
	skin.skel->updateCoords();

	// calculate chest/clavicle offsets and abdomen rotation
	parent = skin.skel->transforms.getT("abdomenT")->globalCoord.mat.inverse();
	Vec3d neck = parent * (0.5 * (ms.markers[34].curPos() + ms.markers[35].curPos()));
	right = parent * (0.5 * (ms.markers[4].curPos() + ms.markers[5].curPos()));
	left = parent * (0.5 * (ms.markers[13].curPos() + ms.markers[14].curPos()));
	axis = left - right;
	r = determineRotation(neck, Vec3d(0, 0, neck.length()), 
		axis, Vec3d(0, -axis.length(), 0));
	tQuat = (SkelQuatRotation*)skin.skel->transforms.getT("abdomenQ");
	r.w *= -1;
	tQuat->curQuat = r;
	r.w *= -1;
	tTrans = (SkelTranslation*)skin.skel->transforms.getT("chestT");
	neck = r.toMatrixD() * neck;
	neck[1] = 0;
	tTrans->curVal = neck;
	tTrans = (SkelTranslation*)skin.skel->transforms.getT("leftT");
	left = r.toMatrixD() * left;
	right = r.toMatrixD() * right;
	right[1] *= -1;
	tTrans->curVal = 0.5 * (left + right);
	skin.skel->updateCoords();

	// calculate clavicle offsets and rotation
	parent = skin.skel->transforms.getT("rightT")->globalCoord.mat.inverse();
	right = parent * (0.5 * (ms.markers[6].curPos() + ms.markers[7].curPos()));
	parent = skin.skel->transforms.getT("leftT")->globalCoord.mat.inverse();
	left = parent * (0.5 * (ms.markers[15].curPos() + ms.markers[16].curPos()));
	len = 0.5 * (right.length() + left.length());
	right.normalize(); right *= len;
	left.normalize(); left *= len;
	angle = atan2(right[0], right[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("rClavicleZR");
	tEuler->curAngle = -angle;
	r = QuatNorm(0, 0, -angle);
	right = r.toMatrixD() * right;
	angle = atan2(right[2], right[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("rClavicleXR");
	tEuler->curAngle = angle;
	angle = atan2(-left[0], -left[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("lClavicleZR");
	tEuler->curAngle = -angle;
	r = QuatNorm(0, 0, -angle);
	left = r.toMatrixD() * left;
	angle = atan2(-left[2], -left[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("lClavicleXR");
	tEuler->curAngle = angle;
	tTrans = (SkelTranslation*)skin.skel->transforms.getT("lClavicleT");
	tTrans->curVal[0] = 0;
	tTrans->curVal[1] = -len;
	tTrans->curVal[2] = 0;
	skin.skel->updateCoords();

	// calculate upper arm offsets and shoulder rotation
	parent = skin.skel->transforms.getT("rClavicleT")->globalCoord.mat.inverse();
	right = parent * (0.5 * (ms.markers[8].curPos() + ms.markers[9].curPos()));
	parent[0][3] = 0; parent[1][3] = 0; parent[2][3] = 0;
	axis = parent * (ms.markers[8].curPos() - ms.markers[9].curPos());
	parent = skin.skel->transforms.getT("lClavicleT")->globalCoord.mat.inverse();
	left = parent * (0.5 * (ms.markers[17].curPos() + ms.markers[18].curPos()));
	len = 0.5 * (right.length() + left.length());
	right.normalize(); right *= len;
	left.normalize(); left *= len;
	angle = atan2(right[2], right[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("rShoulderXR");
	tEuler->curAngle = angle;
	r = QuatNorm(angle, 0, 0);
	right = r.toMatrixD() * right;
	axis = r.toMatrixD() * axis;
	angle = atan2(-right[0], right[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("rShoulderZR");
	tEuler->curAngle = angle;
	r = QuatNorm(0, 0, angle);
	right = r.toMatrixD() * right;
	axis = r.toMatrixD() * axis;
	angle = atan2(-axis[2], axis[0]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("rShoulderYR");
	tEuler->curAngle = angle;
	r = QuatNorm(0, angle, 0);
	//axis = r.toMatrixD() * axis;
	parent[0][3] = 0; parent[1][3] = 0; parent[2][3] = 0;
	axis = parent * (ms.markers[17].curPos() - ms.markers[18].curPos());
	angle = atan2(-left[2], -left[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("lShoulderXR");
	tEuler->curAngle = angle;
	r = QuatNorm(angle, 0, 0);
	left = r.toMatrixD() * left;
	axis = r.toMatrixD() * axis;
	angle = atan2(left[0], -left[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("lShoulderZR");
	tEuler->curAngle = angle;
	r = QuatNorm(0, 0, angle);
	left = r.toMatrixD() * left;
	axis = r.toMatrixD() * axis;
	angle = atan2(-axis[2], axis[0]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("lShoulderYR");
	tEuler->curAngle = angle;
	r = QuatNorm(0, angle, 0);
	axis = r.toMatrixD() * axis;
	tTrans = (SkelTranslation*)skin.skel->transforms.getT("lUpperArmT");
	tTrans->curVal[1] = -len;
	skin.skel->updateCoords();

	// calculate forearm offsets and elbow rotation
	parent = skin.skel->transforms.getT("rUpperArmT")->globalCoord.mat.inverse();
	right = parent * (0.5 * (ms.markers[10].curPos() + ms.markers[11].curPos()));
	parent[0][3] = 0; parent[1][3] = 0; parent[2][3] = 0;
	axis = parent * (ms.markers[10].curPos() - ms.markers[11].curPos());
	parent = skin.skel->transforms.getT("lUpperArmT")->globalCoord.mat.inverse();
	left = parent * (0.5 * (ms.markers[19].curPos() + ms.markers[20].curPos()));
	len = 0.5 * (right.length() + left.length());
	right.normalize(); right *= len;
	left.normalize(); left *= len;
	// right
	// carrying angle
	angle = atan2(-right[0], right[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("rCarryingA");
	tEuler->curAngle = angle;
	r = QuatNorm(0, 0, angle);
	right = r.toMatrixD() * right;
	axis = r.toMatrixD() * axis;
	// elbow angle
	angle = atan2(right[2], right[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("rElbowA");
	tEuler->curAngle = angle;
	r = QuatNorm(angle, 0, 0);
	right = r.toMatrixD() * right;
	axis = r.toMatrixD() * axis;
	// forearm angle
	angle = atan2(axis[0], axis[2]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("rForearmA");
	tEuler->curAngle = angle;
	r = QuatNorm(0, angle, 0);
	//axis = r.toMatrixD() * axis;
	// left
	parent[0][3] = 0; parent[1][3] = 0; parent[2][3] = 0;
	axis = parent * (ms.markers[19].curPos() - ms.markers[20].curPos());
	// carrying angle
	angle = atan2(left[0], -left[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("lCarryingA");
	tEuler->curAngle = angle;
	r = QuatNorm(0, 0, angle);
	right = r.toMatrixD() * right;
	axis = r.toMatrixD() * axis;
	// elbow angle
	angle = atan2(-left[2], -left[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("lElbowA");
	tEuler->curAngle = angle;
	r = QuatNorm(angle, 0, 0);
	left = r.toMatrixD() * left;
	axis = r.toMatrixD() * axis;
	// forearm angle
	angle = atan2(axis[0], axis[2]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("lForearmA");
	tEuler->curAngle = angle;
	r = QuatNorm(0, angle, 0);
	axis = r.toMatrixD() * axis;
	// trans
	tTrans = (SkelTranslation*)skin.skel->transforms.getT("lForearmT");
	tTrans->curVal[1] = -len;
	skin.skel->updateCoords();

	// calculate hand offsets and wrist rotation
	parent = skin.skel->transforms.getT("rForearm2T")->globalCoord.mat.inverse();
	right = parent * ms.markers[12].curPos();
	parent = skin.skel->transforms.getT("lForearm2T")->globalCoord.mat.inverse();
	left = parent * ms.markers[21].curPos();
	len = 0.5 * (right.length() + left.length());
	right.normalize(); right *= len;
	left.normalize(); left *= len;
	// right
	angle = atan2(-right[0], right[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("rWristA");
	tEuler->curAngle = angle;
	r = QuatNorm(0, 0, angle);
	right = r.toMatrixD() * right;
	// left
	angle = atan2(left[0], -left[1]);
	tEuler = (SkelEulerRotation*)skin.skel->transforms.getT("lWristA");
	tEuler->curAngle = angle;
	r = QuatNorm(0, 0, angle);
	left = r.toMatrixD() * left;
	skin.skel->updateCoords();

	// calculate head offset and neck rotation
	parent = skin.skel->transforms.getT("chestT")->globalCoord.mat.inverse();
	right = parent * ms.markers[33].curPos();
	angle = atan2(right[1], -right[2]);
	r = QuatNorm(angle, 0, 0);
	tQuat = (SkelQuatRotation*)skin.skel->transforms.getT("neckR");
	r.w *= -1;
	tQuat->curQuat = r;
	skin.skel->updateCoords();
/*
	extraMarkers.clear();
	extraMarkers.push_back(0.86 * ms.markers[22].curPos() + 0.14 * ms.markers[23].curPos());
	extraMarkers.push_back(0.14 * ms.markers[22].curPos() + 0.86 * ms.markers[23].curPos());
	extraMarkers.push_back(0.5 * (ms.markers[24].curPos() + ms.markers[25].curPos()));
	extraMarkers.push_back(0.5 * (ms.markers[29].curPos() + ms.markers[30].curPos()));
	extraMarkers.push_back(0.5 * (ms.markers[26].curPos() + ms.markers[27].curPos()));
	extraMarkers.push_back(0.5 * (ms.markers[31].curPos() + ms.markers[32].curPos()));
	extraMarkers.push_back(ms.markers[28].curPos());
	extraMarkers.push_back(ms.markers[33].curPos());
	extraMarkers.push_back(0.5 * (ms.markers[2].curPos() + ms.markers[3].curPos()));
	extraMarkers.push_back(0.5 * (ms.markers[34].curPos() + ms.markers[35].curPos()));
	extraMarkers.push_back(0.5 * (ms.markers[4].curPos() + ms.markers[5].curPos()));
	extraMarkers.push_back(0.5 * (ms.markers[13].curPos() + ms.markers[14].curPos()));
	extraMarkers.push_back(0.5 * (ms.markers[6].curPos() + ms.markers[7].curPos()));
	extraMarkers.push_back(0.5 * (ms.markers[15].curPos() + ms.markers[16].curPos()));
	extraMarkers.push_back(0.5 * (ms.markers[8].curPos() + ms.markers[9].curPos()));
	extraMarkers.push_back(ms.markers[8].curPos());
	extraMarkers.push_back(ms.markers[9].curPos());
	extraMarkers.push_back(0.5 * (ms.markers[17].curPos() + ms.markers[18].curPos()));
	extraMarkers.push_back(ms.markers[17].curPos());
	extraMarkers.push_back(ms.markers[18].curPos());
	extraMarkers.push_back(0.5 * (ms.markers[10].curPos() + ms.markers[11].curPos()));
	extraMarkers.push_back(0.5 * (ms.markers[19].curPos() + ms.markers[20].curPos()));
	extraMarkers.push_back(ms.markers[12].curPos());
	extraMarkers.push_back(ms.markers[21].curPos());
	extraMarkers.push_back(ms.markers[33].curPos());
*/
	redrawVNow();
}

void uSaveSkelMat(const char *params) {
	char msName[80], matName[80], jointName[80], cmd[256];
	double dofs[10];
	params = extractString(params, msName, 80);
	params = extractString(params, matName, 80);
	jointName[0] = 0;
	params = extractString(params, jointName, 80);

	ofstream out, jOut;
	if (!openOFStream(&out, matName, "skeleton matrix")) {
		return;
	}
	if (!(jointName[0] != 0 && openOFStream(&jOut, jointName, "joint matrix"))) {
		return;
	}

	int ch;
	for (ch = 0; ch < dataSet.numCharacters; ch++) {
		int ex = dataSet.charIndex[ch];
		sprintf(cmd, "%s data/csr-4k/%sa.dat", msName, dataSet.examples[ex].charName);
		uSkelFromSurf(cmd);

		if (jointName[0] != 0) {
			int tr;
			for (tr=0; tr < skin.skel->transforms.size(); tr++) {
				SkelTransform *trans = skin.skel->transforms.getT(tr);
				jOut << trans->globalCoord.v << endl;
			}
		}

		sprintf(cmd, "data/skels/%sa.po.txt", dataSet.examples[ex].charName);
		ofstream poseOut;
		if (!openOFStream(&poseOut, cmd, "pose file"))
			return;
		skin.skel->savePose(poseOut);
		poseOut.close();

		int tr, dof;
		for (tr=0; tr < skin.skel->transforms.size(); tr++) {
			SkelTransform *trans = skin.skel->transforms.getT(tr);
			if (trans->isIntrinsic) {
				trans->unloadDofs(dofs);
				for (dof = 0; dof < trans->numDofs(); dof++) {
					out << dofs[dof] << " ";
				}
			}
		}
		out << endl;
	}
	out.close();
}

void uLoadPCASkels(const char *params) {
	ifstream in;
	int i, j;
	double d;
	char s[256];

	if (!openIFStream(&in, "data/skels/mu.txt", "mean matrix"))
		return;

	for (i=0; i < skin.skel->numIntrinsicDofs; i++) {
		in >> uSolver.vInt[i];
	}
	in.close();

	if (!openIFStream(&in, "data/skels/u.txt", "u matrix"))
		return;

	for (i=0; i < skin.skel->numIntrinsicDofs; i++) {
		for (j=1; j < uSolver.numComponents; j++) {
			in >> uSolver.vInt[j * skin.skel->numIntrinsicDofs + i];
	
		}
		// skip over extra components:
		for (; j < dataSet.numCharacters; j++) {
			in >> d;
		}
	}
	in.close();

	if (!openIFStream(&in, "data/skels/x.txt", "x matrix"))
		return;

	for (j=1; j < uSolver.numComponents; j++) {
		for (i=0; i < dataSet.numCharacters; i++) {
			in >> dataSet.charMu[i][j];
		}
	}
	for (i=0; i < dataSet.numCharacters; i++) {
		cout << dataSet.charMu[i] << endl;
	}
	in.close();

	for (i=0; i < dataSet.numCharacters; i++) {
		sprintf(s, "data/skels/%sa.po.txt", dataSet.examples[dataSet.charIndex[i]].charName);
		if (!openIFStream(&in, s, "pose file"))
			return;
		skin.skel->loadPose(in);
		in.close();

		int index = 0;
		for (j=0; j < skin.skel->transforms.size(); j++) {
			SkelTransform *trans = skin.skel->transforms.getT(j);
			if (!trans->isIntrinsic) {
				trans->unloadDofs(uSolver.vPose.n + dataSet.charIndex[i] * skin.skel->numPoseDofs + index);
				index += trans->numDofs();
			}
		}
	}
}

void uLoadPoseFiles(const char *params) {
	int ch, j;
	char s[1024];
	ifstream in;
	int mode = 0;
	params = extractInt(params, &mode);

	if (mode == 0) {
		// per-character mu mode (store intrinsics)
		for (ch = 0; ch < dataSet.numCharacters; ch++) {
			int ex = dataSet.charIndex[ch];

			sprintf(s, "data/skels/%sa.po.txt", dataSet.examples[ex].charName);
			if (!openIFStream(&in, s, "pose file"))
				return;
			skin.skel->loadPose(in);
			in.close();

			int index = 0;
			int intIndex = 0;
			for (j=0; j < skin.skel->transforms.size(); j++) {
				SkelTransform *trans = skin.skel->transforms.getT(j);
				if (trans->numDofs() == 0)
					continue;

				if (!trans->isIntrinsic) {
					trans->unloadDofs(uSolver.vPose.n + ex * skin.skel->numPoseDofs + index);
					index += trans->numDofs();
				}
				else {
					trans->unloadDofs(uSolver.vInt.n + ch * skin.skel->numIntrinsicDofs + intIndex);
					intIndex += trans->numDofs();
				}
			}

			// fix dress pose
			uShow(ex);
			int pt;
			for (pt = 0; pt < skin.numPts; pt++) {
				Mat3d m;
				Vec3d base;
				skin.getSkinMat(pt, m, base);

				Vec3d v;
				dataSet.examples[ex].getPt(pt, &v);
				if (v.iszero())
					skin.dressPts[pt] = Vec3d();
				else
					skin.dressPts[pt] = m.inverse() * (v - base);
			}
			uSolver.dressToVars(ch);
		}
	}
}

void uLoadPosesEx(const char *params) {
	int first = 0;
	int last = dataSet.numExamples - 1;
	int i, j;
	ifstream in;
	char s[256], mask[256];
	sprintf(mask, "data/skels/%s%c.po.txt");
	if (params) {
		params = extractInt(params, &first);
		params = extractInt(params, &last);
		params = extractString(params, mask, 256);
	}

	for (i=first; i <= last; i++) {
		sprintf(s, mask, 
			dataSet.examples[i].charName, dataSet.examples[i].poseCh);
		if (!openIFStream(&in, s, "pose file"))
			continue;
		skin.skel->loadPose(in);
		in.close();

		int index = 0;
		for (j=0; j < skin.skel->transforms.size(); j++) {
			SkelTransform *trans = skin.skel->transforms.getT(j);
			if (!trans->isIntrinsic) {
				trans->unloadDofs(uSolver.vPose.n + i * skin.skel->numPoseDofs + index);
				index += trans->numDofs();
			}
		}
	}
}

void uSavePosesEx(const char *params) {
	int first = 0;
	int last = dataSet.numExamples - 1;
	int i, j;
	ofstream out;
	char s[256], mask[256];
	sprintf(mask, "data/skels/%s%c.po.txt");
	if (params) {
		params = extractInt(params, &first);
		params = extractInt(params, &last);
		params = extractString(params, mask, 256);
	}

	for (i=first; i <= last; i++) {
		sprintf(s, mask, 
			dataSet.examples[i].charName, dataSet.examples[i].poseCh);
		if (!openOFStream(&out, s, "pose file"))
			continue;
		cout << "saving " << s << endl;
		uShow(i);
		skin.skel->savePose(out);
		out.close();
	}
}

void uPoseToEx(const char *params) {
	int ex = 0;
	params = extractInt(params, &ex);

	int i;
	int ip = 0;
	int pp = 0;
	uSolver.vInt.zeroElements();
	for (i=0; i < skin.skel->transforms.size(); i++) {
		SkelTransform *curTrans = skin.skel->transforms.getT(i);
		if (curTrans->isIntrinsic) {
			curTrans->unloadDofs(uSolver.vInt.n + ip);
			ip += curTrans->numDofs();
		}
		else {
			curTrans->unloadDofs(uSolver.vPose.n + pp + skin.skel->numPoseDofs * ex);
			pp += curTrans->numDofs();
		}
	}
}

void uFixMirrorLine(const char *params) {
	int pt, comp;
	for (pt = 0; pt < uSolver.numOrigPts; pt++) {
		if (uSolver.mirrorMap[pt] == pt) {
			for (comp = 0; comp < uSolver.numComponents; comp++) {
				cout << uSolver.vDress[(comp * uSolver.numOrigPts + pt) * 3 + 1] << endl;
				uSolver.vDress[(comp * uSolver.numOrigPts + pt) * 3 + 1] = 0;
			}
		}
	}
}

void uAdjustWeights(const char *params) {
	uSolver.adjustWeights();
}

const double MAX_NO_INF = 0.95;

void uLoadPddDefs(const char *params) {
	int i;
	char fname[80];
	UPoseDepDef *curPDD, *mirrorPDD;
	params = extractString(params, fname, 80);
	ifstream in;
	if (!openIFStream(&in, fname, "PDD defs"))
		return;

	while (1) {
		// create a new PDD
		curPDD = new UPoseDepDef();
		int curIndex = (int)skin.pdds.size();

		// load the color data
		Vec3d *colorData = curPDD->load(&skin, uMesh, in);
		if (!colorData)
			break;

		// make a mirror PDD if needed
		mirrorPDD = NULL;
		int mirrorIndex = -1;
		if (uSolver.mirrorTrans[curPDD->transInd] > -1 && 
			 uSolver.mirrorTrans[curPDD->transInd] != curPDD->transInd) {
			mirrorIndex = curIndex + 1;
			mirrorPDD = new UPoseDepDef();
			mirrorPDD->mirrorFrom(curPDD, uSolver.mirrorTrans[curPDD->transInd]);
		}

		int oldNumPddKeys = skin.numPddPtKeys;
		for (i=0; i < skin.numPts; i++) {
			if (colorData[i][1] < MAX_NO_INF) {
				skin.pddPtIndex[i].push_back(curIndex);
				skin.pddPtIndex[i].push_back(skin.numPddPtKeys);

				if (mirrorPDD) {
					int i2 = uSolver.mirrorMap[i];
					skin.pddPtIndex[i2].push_back(mirrorIndex);
					skin.pddPtIndex[i2].push_back(skin.numPddPtKeys);
				}

				skin.numPddPtKeys += curPDD->rbf->mN - 1;
			}
		}
		if (skin.numPddPtKeys != oldNumPddKeys) {
			Vec3d *oldKeys = skin.pddPtKeys;
			skin.pddPtKeys = new Vec3d[skin.numPddPtKeys];
			memcpy(skin.pddPtKeys, oldKeys, sizeof(Vec3d) * oldNumPddKeys);
			delete []oldKeys;
		}
		cout << "now have " << skin.numPddPtKeys << " point keys" << endl;

		oldNumPddKeys = skin.numPddNMKeys;
		for (i=0; i < NORMAL_MAP_SIZE; i++) {
			int tri = texBary[i].tri;
			if (tri >= 0 && (colorData[uMesh->getTri(tri, 0)][1] < MAX_NO_INF ||
				colorData[uMesh->getTri(tri, 1)][1] < MAX_NO_INF ||
				colorData[uMesh->getTri(tri, 2)][1] < MAX_NO_INF)) {
				skin.pddNMIndex[i].push_back(curIndex);
				skin.pddNMIndex[i].push_back(skin.numPddNMKeys);

				if (mirrorPDD) {
					int i2 = (NORMAL_MAP_H - 1 - (i/NORMAL_MAP_W)) * NORMAL_MAP_W + (i % NORMAL_MAP_W);
					skin.pddNMIndex[i2].push_back(mirrorIndex);
					skin.pddNMIndex[i2].push_back(skin.numPddNMKeys);
				}

				skin.numPddNMKeys += curPDD->rbf->mN - 1;
			}
		}
		if (skin.numPddNMKeys != oldNumPddKeys) {
			Vec3d *oldKeys = skin.pddNMKeys;
			skin.pddNMKeys = new Vec3d[skin.numPddNMKeys];
			memcpy(skin.pddNMKeys, oldKeys, sizeof(Vec3d) * oldNumPddKeys);
			delete []oldKeys;
		}
		cout << "now have " << skin.numPddNMKeys << " normal map keys" << endl;

		delete []colorData;

		skin.pdds.push_back(curPDD);
		if (mirrorPDD)
			skin.pdds.push_back(mirrorPDD);
	}
	in.close();

	uSolver.vPDD.resize(skin.numPddPtKeys * 3 * uSolver.numPDComponents, true);
	uSolver.vNMP.resize(skin.numPddNMKeys * 3 * uSolver.numPDComponents, true);

	uSolver.buildPDNeighTable();
}

void uSaveCurPose(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);
	ofstream out;
	if (!openOFStream(&out, fname, "pose file"))
		return;
	skin.skel->savePose(out);
	out.close();
}

void uLoadCurPose(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);
	ifstream in;
	if (!openIFStream(&in, fname, "pose file"))
		return;
	skin.skel->loadPose(in);
	in.close();

	updateSkel(NULL);
}

void uSaveAllPoses(const char *params) {
	char fname[80];
	int i;

	for (i=0; i < dataSet.numExamples; i++) {
		strcpy(fname, dataSet.examples[i].fname);
		if (fname[0] == '-')
			strcpy(fname, "base.mkr");
		int ext = (int)strlen(fname) - 3;
		fname[ext+0] = 'p';
		fname[ext+1] = 'o';
		fname[ext+2] = '.';
		fname[ext+3] = 't';
		fname[ext+4] = 'x';
		fname[ext+5] = 't';
		fname[ext+6] = 0;

		ofstream out;
		if (!openOFStream(&out, fname, "pose file"))
			return;

		uSolver.updateSkel(i);
		skin.skel->savePose(out);
		out.close();
	}
}

void uShowDress(const char *params) {
	showDress = false;
	params = extractBool(params, &showDress);
	updateSkel(NULL);
}

void uDumpWeights(const char *params) {
/*	int i, j, k;
	double weights[20];
//	Vec3d v[KNN_MAX_SAMPLES];
	int vert = 0;
	params = extractInt(params, &vert);

	for (j=0; j < skin.maxInf; j++) {
		int transInd = skin.infJoints[vert*skin.maxInf + j];
		if (transInd < 0)
			continue;

		SkelTransform *curTrans = skin.skel->transforms.getT(transInd);
		cout << "transform: " << curTrans->name << " (" << transInd << ")" << endl;

		if (skin.transRBFs && skin.transRBFs[transInd]) {
			for (i=0; i < dataSet.numExamples; i++) {
				uSolver.updateSkel(i);
				QuatNorm q = curTrans->curCoord.q;
				Vec4d v(q.x, q.y, q.z, q.w);
				skin.transRBFs[transInd]->eval(v.n, weights);
				int numSamples = skin.transPddOfs[transInd+1]-skin.transPddOfs[transInd];

				cout << "(" << skin.transRBFs[transInd]->mN << "): ";
				for (k=0; k < skin.transRBFs[transInd]->mN; k++) {
					cout << weights[k] << " ";
				}
				cout << endl;
			}
		}
	}*/
}

void uDumpStats(const char *params) {
	int pt, inf;
	int numWeights = 0;

	cout << "# of points: " << uMesh->numPts() << endl;

	for (pt=0; pt < uMesh->numPts(); pt++) {
		for (inf=0; inf < skin.maxInf; inf++) {
			int joint = skin.infJoints[pt*skin.maxInf + inf];
			if (joint >= 0) {
				if (fabs(skin.infWeights[pt*skin.maxInf + inf]) > 0.001)
					numWeights++;
			}
		}
	}

	cout << "# of weights: " << numWeights << endl;
	cout << "# of point PDDs: " << skin.numPddPtKeys << endl;
	cout << "# of normal map PDDs: " << skin.numPddNMKeys << endl;

	/*
	int ex;
	int numEx = 0;
	for (ex=0; ex < dataSet.numExamples; ex++) {
		for (pt=0; pt < dataSet.examples[ex].numPts; pt++) {
			double conf;
			Vec3d v;
			dataSet.examples[ex].getPt(pt, &v, &conf);
			if (conf > 0.01)
				numEx++;
		}
	}*/
}

void uDumpPdds(const char *params) {
	int i, j;
	int vert = 0;
	params = extractInt(params, &vert);

	cout << "PDDs for point " << vert << ":" << endl;

	for (j=0; j < skin.pddPtIndex[vert].size(); j += 2) {
		UPoseDepDef *pdd = skin.pdds[skin.pddPtIndex[vert][j]];
		int ofs = skin.pddPtIndex[vert][j+1];
		cout << pdd->transform->name << " (" << ofs << ")" << endl;

		for (i=0; i < pdd->rbf->mN - 1; i++) {
			if (pdd->isMirror)
				cout << i << ": " << skin.pddPtKeys[ofs+i][0] << " " << -skin.pddPtKeys[ofs+i][1] << " " << 
					skin.pddPtKeys[ofs+i][0] << endl;
			else
				cout << i << ": " << skin.pddPtKeys[ofs+i] << endl;
		}
	}
/*
	cout << "PDNs for point " << vert << ":" << endl;

	for (j=0; j < skin.pddNMIndex[vert].size(); j += 2) {
		UPoseDepDef *pdd = skin.pdds[skin.pddNMIndex[vert][j]];
		int ofs = skin.pddNMIndex[vert][j+1];
		cout << pdd->transform->name << " (" << ofs << ")" << endl;

		for (i=0; i < pdd->rbf->mN - 1; i++) {
			if (pdd->isMirror)
				cout << i << ": " << skin.pddNMKeys[ofs+i][0] << " " << -skin.pddNMKeys[ofs+i][1] << " " << 
					skin.pddNMKeys[ofs+i][0] << endl;
			else
				cout << i << ": " << skin.pddNMKeys[ofs+i] << endl;
		}
	}*/
}

void uDumpPddWeights(const char *params) {
	int tr, i;
	for (tr=0; tr < skin.pdds.size(); tr++) {
		cout << skin.pdds[tr]->transform->name << ": " << endl;
		for (i=0; i < skin.pdds[tr]->rbf->mN; i++) {
			cout << skin.pdds[tr]->rbf->curWeights[i] << " ";
		}
		cout << endl;
	}
}

void uDumpPtInfo(int pt) {
	int i;
	cout << "vertex " << pt;
	if (uSolver.mirrorMap)
		cout << "(mirror = " << uSolver.mirrorMap[pt] << ")";
	cout << endl;
	for (i=0; i < skin.maxInf; i++) {
		if (skin.infJoints[pt*skin.maxInf + i] >= 0) {
			SkelTransform *trans = skin.skel->transforms.getT(skin.infJoints[pt*skin.maxInf+i]);
			cout << trans->name << " " << skin.infWeights[pt*skin.maxInf + i] << endl;
		}
	}
	cout << "dress position: " << skin.dressPts[pt] << endl;
}

void uDumpPtInfo(const char *params) {
	int vert = 0;
	params = extractInt(params, &vert);
	uDumpPtInfo(vert);
}

void uDumpJoints(const char *params) {
	int tr;
	for (tr = 0; tr < skin.skel->transforms.size(); tr++) {
		SkelTransform *curTrans = skin.skel->transforms.getT(tr);
		cout << tr << ": " << curTrans->name << ", " <<
			curTrans->curCoord.v << ", " << curTrans->curCoord.q << endl;
	}
}

void uDumpExJoints(const char *params) {
	char jointName[256];
	params = extractString(params, jointName, 256);
	int ex, tr, dof;

	for (ex = 0; ex < dataSet.numExamples; ex++) {
		uShow(ex);
		for (tr = 0; tr < skin.skel->transforms.size(); tr++) {
			SkelTransform *curTrans = skin.skel->transforms.getT(tr);
			if (strcmp(curTrans->name, jointName) == 0) {
				cout << ex << ": ";
				for (dof = 0; dof < curTrans->numDofs(); dof++) {
					cout << curTrans->getDofAddr(dof) << " ";
				}
				cout << endl;
			}
		}
	}
}

void uDumpJointDofs(const char *params) {
	int tr, dof = 0, iDof = 0;
	for (tr = 0; tr < skin.skel->transforms.size(); tr++) {
		SkelTransform *curTrans = skin.skel->transforms.getT(tr);
		if (curTrans->numDofs() > 0) {
			if (!curTrans->isIntrinsic) {
				cout << tr << ": " << curTrans->name << " pose " <<
					dof << "-" << (dof+curTrans->numDofs()-1) << endl;
				dof += curTrans->numDofs();
			}
			else {
				cout << tr << ": " << curTrans->name << " bones " <<
					iDof << "-" << (iDof+curTrans->numDofs()-1) << endl;
				iDof += curTrans->numDofs();
			}
		}
		else
			cout << tr << ": " << curTrans->name << endl;
	}
}

void uDumpJointEx(const char *params) {
	char jointName[256];
	int vert, ex, mode = 0;
	params = extractString(params, jointName, 256);
	params = extractInt(params, &vert);
	params = extractInt(params, &mode);

	SkelTransform *trans = skin.skel->transforms.getT(jointName);
	if (!trans) {
		cout << "joint " << jointName << " not found" << endl;
		return;
	}

	ofstream out;
	if (!openOFStream(&out, "jointData.txt", "joint dump"))
		return;

	if (mode == 0) {
		for (ex = 0; ex < dataSet.numExamples; ex++) {
			uShow(ex);
			out << dataSet.examples[ex].character << " ";
			out << ((SkelEulerRotation*)trans)->curAngle << " ";
			Vec3d v;
			dataSet.examples[ex].getPt(vert, &v);

			Mat3d m3;
			Vec3d v2;
/*			skin.getSkinMat(vert, m3, v2);
			Mat4d m;
			m[0][0] = m3[0][0];
			m[0][1] = m3[0][1];
			m[0][2] = m3[0][2];
			m[0][3] = v2[0];
			m[1][0] = m3[1][0];
			m[1][1] = m3[1][1];
			m[1][2] = m3[1][2];
			m[1][3] = v2[1];
			m[2][0] = m3[2][0];
			m[2][1] = m3[2][1];
			m[2][2] = m3[2][2];
			m[2][3] = v2[2];
			m[3][3] = 1;
			m = m.inverse();*/
			Mat4d m = trans->parentPtr->globalCoord.mat;
			m = m.inverse();
			out << m * v << " " << m * uMesh->getPt(vert) << endl;
		}
	}
	else {
		double angle;
		for (angle = -2.5; angle < 0.05; angle += 0.05) {
			((SkelEulerRotation*)trans)->curAngle = angle;
			updateSkel(NULL);

			out << ((SkelEulerRotation*)trans)->curAngle << " ";
			Vec3d v = uMesh->getPt(vert);

						Mat3d m3;
			Vec3d v2;
/*			skin.getSkinMat(vert, m3, v2);
			Mat4d m;
			m[0][0] = m3[0][0];
			m[0][1] = m3[0][1];
			m[0][2] = m3[0][2];
			m[0][3] = v2[0];
			m[1][0] = m3[1][0];
			m[1][1] = m3[1][1];
			m[1][2] = m3[1][2];
			m[1][3] = v2[1];
			m[2][0] = m3[2][0];
			m[2][1] = m3[2][1];
			m[2][2] = m3[2][2];
			m[2][3] = v2[2];
			m = m.inverse();*/
			Mat4d m = trans->parentPtr->globalCoord.mat;
			m = m.inverse();
			out << m * v << endl;
		}
	}

	out.close();
}

void uDumpMaps(const char *params) {
	unsigned char *data = new unsigned char[256*256*3];
	FILE *fptr;
	int i;

	// save TGA
	fptr = fopen("normalmap.tga", "wb");
	if (fptr == NULL) {
		cout << "can't open " << "normalmap.tga" << endl;
		return;
	}
	for (i=0; i < 256*256; i++) {
		Vec3d v = normalMap[i];
		v.normalize();
		data[i*3+0] = min(max((v[2] + 1.0) * 128, 0), 255);
		data[i*3+1] = min(max((v[1] + 1.0) * 128, 0), 255);
		data[i*3+2] = min(max((v[0] + 1.0) * 128, 0), 255);
	}
	// create tga header
	putc(0,fptr);
	putc(0,fptr);
	putc(2,fptr);                         // uncompressed RGB
	putc(0,fptr); putc(0,fptr);
	putc(0,fptr); putc(0,fptr);
	putc(0,fptr);
	putc(0,fptr); putc(0,fptr);           // X origin
	putc(0,fptr); putc(0,fptr);           // y origin
	putc((256 & 0x00FF),fptr);
	putc((256 & 0xFF00) / 256,fptr);
	putc((256 & 0x00FF),fptr);
	putc((256 & 0xFF00) / 256,fptr);
	putc(24,fptr);                        // 24 bit bitmap
	putc(0,fptr);
	// write the data
	fwrite(data, 256*256*3*sizeof(char), 1, fptr);
	fclose(fptr);

	fptr = fopen("lightmap.tga", "wb");
	if (fptr == NULL) {
		cout << "can't open " << "lightmap.tga" << endl;
		return;
	}
	memcpy(data, curTex, 256*256*3);
	// swap R and B
	for (i=0; i < 256*256; i++) {
		data[i * 3 + 0] ^= data[i * 3 + 2];
		data[i * 3 + 2] ^= data[i * 3 + 0];
		data[i * 3 + 0] ^= data[i * 3 + 2];
	}
	// create tga header
	putc(0,fptr);
	putc(0,fptr);
	putc(2,fptr);                         // uncompressed RGB
	putc(0,fptr); putc(0,fptr);
	putc(0,fptr); putc(0,fptr);
	putc(0,fptr);
	putc(0,fptr); putc(0,fptr);           // X origin
	putc(0,fptr); putc(0,fptr);           // y origin
	putc((256 & 0x00FF),fptr);
	putc((256 & 0xFF00) / 256,fptr);
	putc((256 & 0x00FF),fptr);
	putc((256 & 0xFF00) / 256,fptr);
	putc(24,fptr);                        // 24 bit bitmap
	putc(0,fptr);
	// write the data
	fwrite(data, 256*256*3*sizeof(char), 1, fptr);
	fclose(fptr);

	delete []data;
}

void uZapPdds(const char *params) {
	int i;
	for (i=0; i < skin.numPddPtKeys; i++)
		skin.pddPtKeys[i].zeroElements();
	for (i=0; i < skin.numPddNMKeys; i++)
		skin.pddNMKeys[i].zeroElements();
//	memcpy(skin.dressPts, skin.pddDressPts, sizeof(Vec3d)*skin.numPts);
//	memcpy(skin.baseNormalMap, normalMap, sizeof(Vec3d)*NORMAL_MAP_SIZE);
	skin.updatePts();

	updateMesh();
	redrawV();
}

void uCopyPdds(const char *params) {
	int ch, comp, var;
	params = extractInt(params, &ch);

	uSolver.updatePoints(dataSet.charMu[ch].Ref(), uSolver.numComponents, true);
	updateSkel(NULL);
}

void uDumpVars(const char *params) {
	int start = 0, end = 1;
	params = extractInt(params, &start);
	params = extractInt(params, &end);

	int i;
	for (i=start; i < end; i++) {
		cout << i << ": " << uSolver.curVars[i] << " (" << uSolver.grad[i] << ")" << endl;
	}
}

void uLoadPddMask(const char *params) {
	Vec3d v;
	char fname[80];
	params = extractString(params, fname, 80);

	FILE *f;
	if (!openFile(&f, fname, "rb", "points"))
		return;

	char version;
	fread(&version, sizeof(char), 1, f);

	if (!uSolver.pddMask)
		uSolver.pddMask = new double[uMesh->numPts()];

	if (version == '0') {
		int i;
		fread(&i, sizeof(int), 1, f);

		if (i != uMesh->numPts() * 3) {
			cout << "point size mismatch; expected " << uMesh->numPts()*3 << "; read " << i << endl;
		}
		else {
			for (i=0; i < uMesh->numPts(); i++) {
				fread(v.n, sizeof(double), 3, f);
				uSolver.pddMask[i] = 1.0 - v[1];
			}
		}
		fclose(f);
	}
	else
		cout << "format not supported" << endl;
}

void uGreedyPdds(const char *params) {
	int i, j, k;
	double maxDist, minDist;
	QuatNorm maxQuat, minQuat;
	SkelTransform *joint, *mirrorJoint;
	char jointName[80], mirrorJointName[80];
	mirrorJointName[0] = 0;
	params = extractString(params, jointName, 80);
	params = extractString(params, mirrorJointName, 80);
	
	joint = skin.skel->transforms.getT(jointName);
	mirrorJoint = skin.skel->transforms.getT(mirrorJointName);

	if (!joint) {
		cout << "unknown joint: " << jointName << endl;
		return;
	}

	vector<QuatNorm> quats;
	quats.push_back(QuatNorm(0, 0, 0, 1));

	for (i=0; i < 15; i++) {
		maxDist = 0;

		for (j=0; j < dataSet.numExamples; j++) {
			uSolver.updateSkel(j);
			double dist;
			minDist = 10;

			for (k=0; k < quats.size(); k++) {
				dist = quatDist(quats[k], joint->curCoord.q);
				if (dist < minDist) {
					minDist = dist;
					minQuat = joint->curCoord.q;
				}

				QuatNorm q;
				if (mirrorJoint) {
					q = mirrorJoint->curCoord.q;

					q.y *= -1;
					q.w *= -1;
					dist = quatDist(quats[k], q);
					if (dist < minDist) {
						minDist = dist;
						minQuat = q;
					}
				}
			}

			if (minDist > maxDist) {
				maxDist = minDist;
				maxQuat = minQuat;
			}
		}

		quats.push_back(maxQuat);
		cout << maxQuat << " - " << maxDist << endl;
/*		if (!mirrorJoint) {
			maxQuat.y *= -1;
			maxQuat.w *= -1;
			quats.push_back(maxQuat);
			cout << maxQuat << " (mirror)" << endl;
		}*/
	}
}

/*
void uFixShoulder(const char *params) {
	int ex, ofs;
	QuatNorm q;

	for (ex=0; ex < dataSet.numExamples; ex++) {
		ofs = uSolver.poseDofPos + ex*skin.skel->numPoseDofs + 
			skin.dofIndex[skin.skel->transforms.lookupName("lShoulderQ")];
		q = QuatNorm(uSolver.vars[ofs + 0], uSolver.vars[ofs + 1], uSolver.vars[ofs + 2],uSolver.vars[ofs + 3]);
		q.normalize();
		q = slerp(q, QuatNorm(0,0,0,1), 0.5);
        uSolver.vars[ofs + 0] = q.x;
		uSolver.vars[ofs + 1] = q.y;
		uSolver.vars[ofs + 2] = q.z;
		uSolver.vars[ofs + 3] = q.w;

		ofs = uSolver.poseDofPos + ex*skin.skel->numPoseDofs + 
		skin.dofIndex[skin.skel->transforms.lookupName("lClavicleQ")];
		q = QuatNorm(uSolver.vars[ofs + 0], uSolver.vars[ofs + 1], uSolver.vars[ofs + 2],uSolver.vars[ofs + 3]);
		q.normalize();
		q = slerp(q, QuatNorm(0,0,0,1), 0.5);
        uSolver.vars[ofs + 0] = q.x;
		uSolver.vars[ofs + 1] = q.y;
		uSolver.vars[ofs + 2] = q.z;
		uSolver.vars[ofs + 3] = q.w;

		ofs = uSolver.poseDofPos + ex*skin.skel->numPoseDofs + 
			skin.dofIndex[skin.skel->transforms.lookupName("rShoulderQ")];
		q = QuatNorm(uSolver.vars[ofs + 0], uSolver.vars[ofs + 1], uSolver.vars[ofs + 2],uSolver.vars[ofs + 3]);
		q.normalize();
		q = slerp(q, QuatNorm(0,0,0,1), 0.5);
        uSolver.vars[ofs + 0] = q.x;
		uSolver.vars[ofs + 1] = q.y;
		uSolver.vars[ofs + 2] = q.z;
		uSolver.vars[ofs + 3] = q.w;

		ofs = uSolver.poseDofPos + ex*skin.skel->numPoseDofs + 
			skin.dofIndex[skin.skel->transforms.lookupName("rClavicleQ")];
		q = QuatNorm(uSolver.vars[ofs + 0], uSolver.vars[ofs + 1], uSolver.vars[ofs + 2],uSolver.vars[ofs + 3]);
		q.normalize();
		q = slerp(q, QuatNorm(0,0,0,1), 0.5);
        uSolver.vars[ofs + 0] = q.x;
		uSolver.vars[ofs + 1] = q.y;
		uSolver.vars[ofs + 2] = q.z;
		uSolver.vars[ofs + 3] = q.w;
	}
}*/

void uLoadMocap(const char *params) {
	char fname[256];
	params = extractString(params, fname, 256);
	mocapPoses.load(fname);
	uniUI->frameVS->maximum(mocapPoses.numFrames-1);

	mocapPoses.frameToSkel(uniUI->frameVS->value(), skin.skel);
	uShowComps();
}

void uFixMocap(const char *params) {
	int ex = 0;
	params = extractInt(params, &ex);
	
	Skeleton *tempSkel = new Skeleton();
	tempSkel->copyFrom(skin.skel);

	mocapPoses.zeroFix();
	if (ex < 0) {
		uShow(-ex);
		mocapPoses.calcFix(skin.skel, tempSkel);
	}
	else {
		mocapPoses.frameToSkel(ex, skin.skel);
		uShowComps();
		mocapPoses.calcFix(tempSkel, skin.skel);
	}
}

void uLoadFix(const char *params) {
	char fname[256];
	params = extractString(params, fname, 255);
	mocapPoses.loadFix(fname);
}

void uSaveFix(const char *params) {
	char fname[256];
	params = extractString(params, fname, 255);
	mocapPoses.saveFix(fname);
}

void uSaveMocap(const char *params) {
	char fName[256], mask[256];
	char format;
	mask[0] =0;
	int start = 0;
	int end = mocapPoses.numFrames;
	int step = 4;
	int special = 0;
	params = extractString(params, mask, 255);
	params = extractInt(params, &start);
	params = extractInt(params, &end);
	params = extractInt(params, &step);
	params = extractInt(params, &special);

//	if (mask[0] == 0)
//		strcpy(mask, "anim/frame%04d.tga");
	format = mask[strlen(mask)-3];

	double *sPCA, *tPCA;
	int comp;
	if (special == 1) {
		sPCA = new double[uSolver.numComponents];
		tPCA = new double[uSolver.numComponents];
		for (comp=1; comp < uSolver.numComponents; comp++) {
			tPCA[comp] = curComps[comp];
		}
	}

	int i, frame = 0;
	for (i=start; i < end; i += step) {
		if (special == 1) {
			int pos = (i - start) % 30;
			if (pos == 0) {
				int newGuy = boundedRand(0, dataSet.numCharacters-0.0001);
				for (comp=1; comp <uSolver.numComponents; comp++) {
					sPCA[comp] = tPCA[comp];
					tPCA[comp] = dataSet.charMu[newGuy][comp];
//					tPCA[comp] = sqrt(-2 * log(boundedRand(0, 1) + 0.0000001)) * cos(boundedRand(0, PI*2)) * 0.6;
				}
			}
			for (comp=1; comp <uSolver.numComponents; comp++) {
				curComps[comp] = sPCA[comp] * (30 - pos) / 30.0 + tPCA[comp] * pos / 30.0;
			}
			uShowComps(true);
		}
		else if (special == 2) {
			float fVals[2];
			fVals[0] = 66;
			fVals[1] = 100 + 250.0 * ((1.0*end - i) / (1.0*end - start));
			setFeatures(fVals, curComps);
			uShowComps(true);
		}

		mocapPoses.frameToSkel(i, skin.skel);
		uShowComps(false);
		redrawVNow();
		uiWait();

		if (strlen(mask) > 1) {
			sprintf(fName, mask, frame);
			cout << "saving " << fName << endl;
			if (format == 't')
				uiScreenshot(fName);
			else
				uMesh->saveFile(fName);
		}
		frame++;
	}
}

void uExToMocap(const char *params) {
	int frame;
	char fname[256];
	fname[0] = 0;
	params = extractString(params, fname, 256);

	mocapPoses.init(dataSet.numExamples, mocapPoses.numVars);
	for (frame = 0; frame < dataSet.numExamples; frame++) {
		uSolver.updateSkel(frame);
		mocapPoses.skelToFrame(frame, skin.skel);
	}

	if (strlen(fname) > 0) {
		mocapPoses.savePoses(fname);
	}
}

void uSaveCurPDs(const char *params) {
	FILE *f;
	char fname[256];
	params = extractString(params, fname, 256);

	if (!openFile(&f, fname, "wb", "pd file"))
		return;

	fwrite(skin.pddPtKeys, sizeof(double), skin.numPddPtKeys * 3, f);
	fwrite(skin.pddNMKeys, sizeof(double), skin.numPddNMKeys * 3, f);

	fclose(f);
}

void uLoadCurPDs(const char *params) {
	FILE *f;
	char fname[256];
	params = extractString(params, fname, 256);

	if (!openFile(&f, fname, "rb", "pd file"))
		return;

	fread(skin.pddPtKeys, sizeof(double), skin.numPddPtKeys * 3, f);
	fread(skin.pddNMKeys, sizeof(double), skin.numPddNMKeys * 3, f);

	fclose(f);
}

void uMirrorEx(const char *params) {
	char inMesh[256], inMkr[256], outMesh[256], outMkr[256];
	params = extractString(params, inMesh, 256);
	params = extractString(params, inMkr, 256);
	params = extractString(params, outMesh, 256);
	params = extractString(params, outMkr, 256);

	int i;
	TriMesh tm;
	MarkerSet mkr;

	if (!tm.loadFile(inMesh))
		return;
	mkr.loadFromMkr(inMkr);

	// mirror the mesh
	for (i=0; i < tm.numPts(); i++) {
		tm.getPt(i)[1] *= -1;
	}
	for (i=0; i < tm.numTris(); i++) {
		int v0 = tm.getTri(i, 0);
		int v1 = tm.getTri(i, 1);
		tm.getTri(i, 1) = v0;
		tm.getTri(i, 0) = v1;
	}

	// mirror the vertices
	for (i=0; i < mkr.numMarkers; i++) {
		mkr.markers[i].pos[1] *= -1;
	}
	// swapping!
	swapVec(mkr.markers[2].pos, mkr.markers[3].pos);
	swapVec(mkr.markers[5].pos, mkr.markers[7].pos);
	swapVec(mkr.markers[6].pos, mkr.markers[8].pos);
	swapVec(mkr.markers[10].pos, mkr.markers[12].pos);
	swapVec(mkr.markers[13].pos, mkr.markers[14].pos);
	swapVec(mkr.markers[16].pos, mkr.markers[18].pos);
	swapVec(mkr.markers[17].pos, mkr.markers[19].pos);
	swapVec(mkr.markers[20].pos, mkr.markers[22].pos);
	swapVec(mkr.markers[21].pos, mkr.markers[23].pos);
	swapVec(mkr.markers[26].pos, mkr.markers[27].pos);
	swapVec(mkr.markers[29].pos, mkr.markers[41].pos);
	swapVec(mkr.markers[30].pos, mkr.markers[42].pos);
	swapVec(mkr.markers[31].pos, mkr.markers[43].pos);
	swapVec(mkr.markers[32].pos, mkr.markers[44].pos);
	swapVec(mkr.markers[33].pos, mkr.markers[45].pos);
	swapVec(mkr.markers[34].pos, mkr.markers[46].pos);
	swapVec(mkr.markers[35].pos, mkr.markers[47].pos);
	swapVec(mkr.markers[36].pos, mkr.markers[48].pos);
	swapVec(mkr.markers[37].pos, mkr.markers[49].pos);
	swapVec(mkr.markers[38].pos, mkr.markers[50].pos);
	swapVec(mkr.markers[39].pos, mkr.markers[51].pos);
	swapVec(mkr.markers[40].pos, mkr.markers[52].pos);
	swapVec(mkr.markers[53].pos, mkr.markers[63].pos);
	swapVec(mkr.markers[54].pos, mkr.markers[64].pos);
	swapVec(mkr.markers[55].pos, mkr.markers[65].pos);
	swapVec(mkr.markers[56].pos, mkr.markers[66].pos);
	swapVec(mkr.markers[57].pos, mkr.markers[67].pos);
	swapVec(mkr.markers[58].pos, mkr.markers[68].pos);
	swapVec(mkr.markers[59].pos, mkr.markers[69].pos);
	swapVec(mkr.markers[60].pos, mkr.markers[70].pos);
	swapVec(mkr.markers[61].pos, mkr.markers[71].pos);
	swapVec(mkr.markers[62].pos, mkr.markers[72].pos);
	if (mkr.numMarkers > 77)
		swapVec(mkr.markers[76].pos, mkr.markers[77].pos);
	if (mkr.numMarkers > 79)
		swapVec(mkr.markers[78].pos, mkr.markers[79].pos);


	tm.savePly(outMesh);
	mkr.saveToMkr(outMkr);
}

void uMatchAll(const char *params) {
	char res[80], pose[80], csrID[80], cmd[256];
	int ex, i, ch = 0;
	int iter = 25;
	int startAt = 0, stopAt = dataSet.numExamples;
	params = extractString(params, res, 80);
	params = extractInt(params, &ch);
	params = extractInt(params, &startAt);
	params = extractInt(params, &stopAt);

	startAt += dataSet.charIndex[ch];
	stopAt += dataSet.charIndex[ch];

	for (ex=startAt; ex < stopAt; ex++) {
		if (dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname)-5] == 'z') {
			pose[0] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname)-5];
			pose[1] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname)-4];
			pose[2] = 0;
			csrID[0] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname)-9];
			csrID[1] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname)-8];
			csrID[2] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname)-7];
			csrID[3] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname)-6];
			csrID[4] = 0;
		}
		else {
			pose[0] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname)-4];
			pose[1] = 0;
			csrID[0] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname)-8];
			csrID[1] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname)-7];
			csrID[2] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname)-6];
			csrID[3] = dataSet.examples[ex].fname[strlen(dataSet.examples[ex].fname)-5];
			csrID[4] = 0;
		}

		if (csrID[0] == '9') {
			sprintf(cmd, "uEdgeMatchInit %d data/csr/csr%s%s.ply data/csr/csr%s%s.mkr data/j%s.mkr",
				ex, csrID, pose, csrID, pose, res);
		}
		else {
			sprintf(cmd, "uEdgeMatchInit %d data/csr/csr%s%s.ply data/csr/csr%s%s.lnd data/j%s.mkr",
				ex, csrID, pose, csrID, pose, res);
		}
		processCommand(cmd);

		sprintf(cmd, "uEdgeMatchLoadNW data/bend-%s.tex", res);
		processCommand(cmd);
		sprintf(cmd, "uEdgeMatchLoadSW data/lock-%s.dat", res);
		processCommand(cmd);
//		sprintf(cmd, "uEdgeMatchLoadRestrict data/arm.tex", res);
//		sprintf(cmd, "uEdgeMatchSolve 0 1 0 10 %d", iter/5);
//		processCommand(cmd);
//		processCommand("uEdgeMatchNorm");
//		sprintf(cmd, "uEdgeMatchSolve 0 1 0 10 %d", iter);
//		processCommand(cmd); 
		//uEdgeMatchSolve 1 0.1 0 1000 %d
		sprintf(cmd, "uEdgeMatchSolve 1 1 1 100 %d", iter); // 0.02
		for (i=0; i < 6; i++) {
			processCommand(cmd);
			processCommand("uEdgeMatchNorm");
		}
//		sprintf(cmd, "uEdgeMatchSolve 10 0.1 0 10 %d", iter);
		for (i=0; i < 16; i++) {
			sprintf(cmd, "uEdgeMatchSolve %f 0.1 %f 0.1 %d", 
				15.0 * (i+1), 1.0 / (i+1), iter);

			processCommand(cmd);
			processCommand("uEdgeMatchNorm");
		}

		sprintf(cmd, "uEdgeMatchSaveEx data/csr-%s/csr%s%s.ex", res, csrID, pose);
		processCommand(cmd);
	}
}

void uAutoNW(const char *params) {
/*	int i, j;
	int ch = 0;
	char fname[80];
	fname[0] = 0;
	params = extractString(params, fname, 80);

	double *weights = new double[uMesh->numPts()];
	memset(weights, 0, sizeof(double)*uMesh->numPts());
	QuatNorm *quats = new QuatNorm[uMesh->numPts()];

	uShow(dataSet.charIndex[ch]);

	for (i=0; i < uMesh->numPts(); i++) {
		Vec3d exPt;
		dataSet.examples[dataSet.charIndex[ch]].getPt(i, &exPt);
		uMesh->getPt(i) = exPt;
		skin.dressPts[i] = exPt;
	}
	uMesh->calcNormals();
	TriMesh *origMesh = new TriMesh(uMesh);
	origMesh->calcNormals();

	vector<int> *tmNeigh = findTMNeighbors(uMesh);

	int ex = 0;
	while (dataSet.examples[ex].character == 0) {
		// skin to new pose
		uSolver.updateSkel(ex);
	//	don't uSolver.updatePoints(ch) because we're using the dress points, not the reconstructed points
		skin.skel->updateCoords();
	//	updateSkel(NULL);
		skin.updateMats();
		skin.updatePts();
		// copy points into the TriMesh
		if (uMesh) {
			for (i=0; i < uMesh->numPts(); i++) {
				uMesh->getPt(i) = skin.curPts[i];
				uMesh->getPtColor(i) = Vec3d(0.8, 0.8, 0.8);
			}
			uMesh->calcNormals();
		}

		// calculate rotations
		const int MAX_NEIGH = 10;
		for (i=0; i < uMesh->numPts(); i++) {
			int numInvolved = 0;
			Vec3d v0[MAX_NEIGH], v1[MAX_NEIGH];
			int vInd[MAX_NEIGH];
			vInd[numInvolved++] = i;
			for (j=0; j < (int)tmNeigh[i].size(); j++) {
				vInd[numInvolved++] = tmNeigh[i][j];
				if (numInvolved >= MAX_NEIGH)
					break;
			}
			for (j=1; j < numInvolved; j++) {
				int k;
				for (k=0; k < (int)tmNeigh[vInd[j]].size(); k++) {
					int neigh = tmNeigh[vInd[j]][k];
					int x;
					for (x=0; x < numInvolved; x++)
						if (vInd[x] == neigh) break;
					if (x == numInvolved)
						vInd[numInvolved++] = neigh;
					if (numInvolved >= MAX_NEIGH)
						break;
				}
				if (numInvolved >= MAX_NEIGH)
					break;
			}


			for (j=0; j < numInvolved; j++) {
				v0[j] = origMesh->getPt(vInd[j]);
				v1[j] = uMesh->getPt(vInd[j]);
			}

			Mat3d m = ptsRot(v0, v1, numInvolved);
			quats[i] = matToQuat(m);
		}

		ex++;
	}
*/	
}

void uEdgeMatchInit(const char *params) {
	int i;
	int example = 0;
	char plyName[255], markerName[255], mRefsName[255], s[80];

	params = extractInt(params, &example);
	params = extractString(params, plyName, 255);
	params = extractString(params, markerName, 255);
	params = extractString(params, mRefsName, 255);

	static TriMesh *targetMesh = NULL;
	
	if (targetMesh)
		delete targetMesh;
	targetMesh = new TriMesh();
	if (!targetMesh->loadFile(plyName))
		return;
	targetMesh->calcNormals();
	targetMesh->calcHBB(16);

	// set up dress pose
/*	int ch = dataSet.examples[example].character;
	uShow(dataSet.charIndex[ch]);
	uMesh->calcNormals();
	TriMesh *origMesh = new TriMesh(uMesh);
	origMesh->calcNormals();
*/
	// skin to new pose
	uSolver.updateSkel(example);
//	don't uSolver.updatePoints(ch) because we're using the dress points, not the reconstructed points
	skin.skel->updateCoords();
//	updateSkel(NULL);
	skin.updateMats();
	skin.updatePts();
	// copy points into the TriMesh
	updateMesh();

	vector<int> *neighbors = findTMNeighbors(uMesh);

	// initialize goal function
	if (!edgeMatchGF) {
//		edgeMatchGF = new SkinMatchGF(uMesh, neighbors, origMesh);
		edgeMatchGF = new SkinMatchGF(uMesh, neighbors, &skin);

		// load marker refs
		edgeMatchGF->srcMarkers = new MarkerSet();
		edgeMatchGF->srcMarkers->load(mRefsName);
		int i;
		for (i=0; i < edgeMatchGF->srcMarkers->numMarkers; i++)
			edgeMatchGF->srcMarkers->markers[i].mesh = uMesh;

//		for (i=0; i < cvGF->cv->tm->numPts(); i++) {
//			edgeMatchGF->origVerts[i] = matchTM->getPt(i);
//		}
	}
	else {
//		edgeMatchGF->newMatch(uMesh, origMesh);
		edgeMatchGF->newMatch();
	}
	
//	cout << "calculating geodesics..." << endl;
/*	int numGeo = 0;
	for (i=0; i < origMesh->numPts(); i++) {
		int j;
		for (j=i+1; j < origMesh->numPts(); j++) {
			double dist = (origMesh->getPt(i) - origMesh->getPt(j)).length();
			if (dist > 0.2)
				continue;

			double compat = 0;
			int inf1, inf2;
			for (inf1 = 0; inf1 < skin.maxInf; inf1++) {
				for (inf2 = 0; inf2 < skin.maxInf; inf2++) {
					if ((skin.infJoints[i*skin.maxInf+inf1] == skin.infJoints[j*skin.maxInf+inf2]) && 
						skin.infJoints[i*skin.maxInf+inf1] > 0) {
						compat += skin.infWeights[i*skin.maxInf+inf1] * skin.infWeights[j*skin.maxInf+inf2];
					}
				}
			}
			if (compat > 0.2) {
				edgeMatchGF->geodesics[i].push_back(TMNeigh(i, j, dist, compat));
				numGeo++;
			}
		}
	}
	cout << "added " << numGeo << " geodesics" << endl;*/
/*	edgeMatchGF->geodesics = findTMGeodesics(origMesh, 0.2);
	for (i=0; i < origMesh->numPts(); i++) {
		int j;
		for (j=0; j < edgeMatchGF->geodesics[i].size(); j++) {
			baAssert(edgeMatchGF->geodesics[i][j].vert >= 0, "bad vertex", true);
			baAssert(edgeMatchGF->geodesics[i][j].vert2 >= 0, "bad vertex", true);
			baAssert(edgeMatchGF->geodesics[i][j].vert < origMesh->numPts(), "bad vertex", true);
			baAssert(edgeMatchGF->geodesics[i][j].vert2 < origMesh->numPts(), "bad vertex", true);
			baAssert(edgeMatchGF->geodesics[i][j].vert2 != edgeMatchGF->geodesics[i][j].vert, "bad vertex", true);

			edgeMatchGF->geodesics[i][j].weight = 0.5*(0.2 - edgeMatchGF->geodesics[i][j].dist);
			edgeMatchGF->geodesics[i][j].dist = 
				(origMesh->getPt(edgeMatchGF->geodesics[i][j].vert) - 
				origMesh->getPt(edgeMatchGF->geodesics[i][j].vert2)).length();
		}
	}
	calcInnerSprings(origMesh);*/

	// load markers
	edgeMatchGF->markers = new MarkerSet();
	if (markerName[strlen(markerName)-3] == 'm')
		edgeMatchGF->markers->loadFromMkr(markerName);
	else {
		cout << "loading from landmarks..." << endl;
		edgeMatchGF->markers->loadFromLandmarks(markerName);
		// wipe out extraneous markers
		//int i;
		//for (i=73; i < 85; i++) {
			//edgeMatchGF->markers->markers[i].pos = Vec3d();
		//}
	}

	// ignore right dactylion in c pose
	if (markerName[strlen(markerName)-5] == 'c')
		cout << "ignoring right dactylion" << endl;
		edgeMatchGF->markers->markers[38].pos = Vec3d();
	
	// ignore patella
//	matchGF->markers->markers[73].pos = Vec3d();
//	matchGF->markers->markers[74].pos = Vec3d();

//	matchMeshes[example].calcNormals();
//	matchMeshes[example].calcHBB(16);
	edgeMatchGF->prepareTriMesh(targetMesh);
	edgeMatchGF->zeroDeformation();

	/*
	// color based on # of edges
	int *counts = new int[edgeMatchGF->curMesh->numPts()];
	memset(counts, 0, sizeof(int)*edgeMatchGF->curMesh->numPts());
	for (i=0; i < edgeMatchGF->numEdges; i++) {
		counts[edgeMatchGF->templateEdges[i].v0]++;
		counts[edgeMatchGF->templateEdges[i].v1]++;
	}
	for (i=0; i < edgeMatchGF->curMesh->numPts(); i++) {
		if (counts[i] < 2)
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d();
		else if (counts[i] == 3)
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d(1,0,0);
		else if (counts[i] == 4)
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d(0,1,0);
		else if (counts[i] == 5)
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d(0,0,1);
		else if (counts[i] == 6)
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d(1,1,0);
		else if (counts[i] == 7)
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d(0,1,1);
		else 
			edgeMatchGF->curMesh->getPtColor(i) = Vec3d(1,1,1);
	}*/

	// initialization
/*	cvGF->cv->updateTM(cvGF->examples[example].trans);
	for (i=0; i < edgeMatchGF->curMesh->numPts(); i++) {
		edgeMatchGF->vars[i*3 + 0] = cvGF->cv->tm->getPt(i)[0];
		edgeMatchGF->vars[i*3 + 1] = cvGF->cv->tm->getPt(i)[1];
		edgeMatchGF->vars[i*3 + 2] = cvGF->cv->tm->getPt(i)[2];
	}*/
}

void uEdgeMatchLoadNW(const char *params) {
	Vec3d v;
	char fname[80];
	params = extractString(params, fname, 80);

	FILE *f;
	if (!openFile(&f, fname, "rb", "points"))
		return;

	char version;
	fread(&version, sizeof(char), 1, f);

	if (version == '0') {
		int i;
		fread(&i, sizeof(int), 1, f);

		if (i != uMesh->numPts() * 3) {
			cout << "point size mismatch; expected " << uMesh->numPts()*3 << "; read " << i << endl;
		}
		else {
			for (i=0; i < uMesh->numPts(); i++) {
				fread(v.n, sizeof(double), 3, f);
				if (edgeMatchGF) {
					if (v[1] < 0)
						v[1] = 0;
					if (v[1] > 1)
						v[1] = 1;
					edgeMatchGF->neighWeights[i] = v[1];
				}
			}
		}
		fclose(f);
	}
	else
		cout << "format not supported" << endl;
}

void uEdgeMatchLoadSW(const char *params) {
	Vec3d v;
	char fname[80];
	params = extractString(params, fname, 80);

	FILE *f;
	if (!openFile(&f, fname, "rb", "points"))
		return;

	char version;
	fread(&version, sizeof(char), 1, f);

	if (version == '0') {
		int i;
		fread(&i, sizeof(int), 1, f);

		if (i != uMesh->numPts() * 3) {
			cout << "point size mismatch; expected " << uMesh->numPts()*3 << "; read " << i << endl;
		}
		else {
			for (i=0; i < uMesh->numPts(); i++) {
				fread(v.n, sizeof(double), 3, f);
				if (edgeMatchGF)
					edgeMatchGF->surfWeights[i] = v[1];
			}
		}
		fclose(f);
	}
	else
		cout << "format not supported" << endl;
}

void uEdgeMatchLoadRestrict(const char *params) {
	Vec3d v;
	char fname[80];
	params = extractString(params, fname, 80);

	if (!edgeMatchGF) {
		cout << "edgeMatchGF not initialized; aborting uEdgeMatchLoadRestrict" << endl;
	}
	edgeMatchGF->restrict = new bool[skin.numPts];

	FILE *f;
	if (!openFile(&f, fname, "rb", "points"))
		return;

	char version;
	fread(&version, sizeof(char), 1, f);

	if (version == '0') {
		int i;
		fread(&i, sizeof(int), 1, f);

		if (i != uMesh->numPts() * 3) {
			cout << "point size mismatch; expected " << uMesh->numPts()*3 << "; read " << i << endl;
		}
		else {
			for (i=0; i < uMesh->numPts(); i++) {
				fread(v.n, sizeof(double), 3, f);
				edgeMatchGF->restrict[i] = (v[2] > 0.98);
			}
		}
		fclose(f);
	}
	else
		cout << "format not supported" << endl;
}

void uEdgeMatchMap(const char *params) {
	char fname[256];
	FILE *f;
	int i, j;

	if (!edgeMatchGF) {
		cout << "edgeMatchGF not initialized; aborting uEdgeMatchMap" << endl;
		return;
	}

	params = extractString(params, fname, 256);
	if (!openFile(&f, fname, "rb", "mapping file"))
		return;

	int numSmallPts, numBigPts;
	int *mapPts;
	float *mapWeights;

	// load mapping header
	fread(&numSmallPts, sizeof(int), 1, f);
	fread(&numBigPts, sizeof(int), 1, f);
	if (numBigPts != skin.numPts) {
		cout << "point size mismatch: mapping is for " << numBigPts << "; skin has " << skin.numPts << endl;
		return;
	}

	// load mapping data
	mapPts = new int[numBigPts * 4];
	mapWeights = new float[numBigPts * 4];
	for (i=0; i < numBigPts; i++) {
		int numSources;
		fread(&numSources, sizeof(int), 1, f);

		for (j=0; j < numSources; j++) {
			fread(mapPts + 4*i + j, sizeof(int), 1, f);
			fread(mapWeights + 4*i + j, sizeof(float), 1, f);
		}
		for (; j < 4; j++) {
			mapPts[4*i+j] = -1;
			mapWeights[4*i+j] = 0;
		}
	}

	int *lo2hi = new int[numSmallPts];
	for (j=0; j < numSmallPts; j++)
		lo2hi[j] = -1;
	for (i=0; i < numBigPts; i++) {
		for (j=0; j < 4; j++) {
			if (fabs(mapWeights[4*i+j]-1.0) < 0.01) {
				lo2hi[mapPts[4*i+j]] = i;
			}
		}
	}
	for (j=0; j < numSmallPts; j++) {
		if (lo2hi[j] == -1)
			cout << "warning: no mapping found for " << j << endl;
	}

	edgeMatchGF->baryRecon = new double[numBigPts * 3];
	edgeMatchGF->baryPts = new int[numBigPts * 3];
	for (i=0; i < numBigPts; i++) {
		for (j=0; j < 3; j++) {
			edgeMatchGF->baryPts[i*3+j] = lo2hi[mapPts[4*i+j]];
			if (edgeMatchGF->baryPts[i*3+j] < 0)
				cout << "warning: not enough points for " << i << endl;
			edgeMatchGF->baryRecon[i*3+j] = mapWeights[4*i+j];
		}
	}
	for (i=0; i < numSmallPts; i++) {
		edgeMatchGF->baryPts[lo2hi[i]*3+0] = -1;
	}

	for (i=0; i < numBigPts; i++) {
		if (edgeMatchGF->baryPts[i*3+0] < 0)
			continue;

		Vec3d ofs;
		for (j=0; j < 3; j++) {
			ofs += edgeMatchGF->baryRecon[j] * Vec3d(
				edgeMatchGF->vars[edgeMatchGF->baryPts[j*3+0]*3+0],
				edgeMatchGF->vars[edgeMatchGF->baryPts[j*3+0]*3+1],
				edgeMatchGF->vars[edgeMatchGF->baryPts[j*3+0]*3+2]);
		}
		edgeMatchGF->vars[i*3+0] -= ofs[0];
		edgeMatchGF->vars[i*3+1] -= ofs[1];
		edgeMatchGF->vars[i*3+2] -= ofs[2];
	}
}

void uEdgeMatchNorm(const char *params) {
	int i, j;
	if (!edgeMatchGF)
		return;
/*	for (i=0; i < uMesh->numPts(); i++) {
		double sum = 0;
		for (j=3; j < 7; j++) {
			sum += sqr(edgeMatchGF->vars[i*7+j]);
		}
		sum = sqrt(sum);
		for (j=3; j < 7; j++) {
			edgeMatchGF->vars[i*7+j] /= sum;
		}
	}*/
}

static LBFGSSolver *solver = NULL;

void uEdgeMatchSolve(const char *params) {

	if (!edgeMatchGF) {
		cout << "warning: goal function not initialized in uEdgeMatchSolve" << endl;
		return;
	}
	if (solver) {
		cout << "warning: solver already running in uEdgeMatchSolve" << endl;
		return;
	}

	edgeMatchGF->findTargets = true;

	params = extractDouble(params, &(edgeMatchGF->surfaceMatchWeight));
	params = extractDouble(params, &(edgeMatchGF->smoothnessWeight));
	params = extractDouble(params, &(edgeMatchGF->bendWeight));
	params = extractDouble(params, &(edgeMatchGF->markerMatchWeight));
	int maxIter = 1000;
	params = extractInt(params, &maxIter);
	int minV = 0, maxV = 10;
	params = extractInt(params, &minV);
	params = extractInt(params, &maxV);

	if (maxIter == -1) {
		runTest(edgeMatchGF, edgeMatchGF->vars, minV, maxV);
		return;
	}

	solver = new LBFGSSolver(edgeMatchGF);
	solver->solve(1e+3, 1e-7, edgeMatchGF->vars, maxIter);
	edgeMatchGF->applyDef(edgeMatchGF->vars);
	redrawV();
//	solver->solve(1e+3, 1e-5, edgeMatchGF->vars, maxIter);

	delete solver;
	solver = NULL;
}

void uEdgeMatchStopSolve(const char *params) {
	if (solver)
		solver->stopNow = true;
}

void uEdgeMatchSaveEx(const char *params) {
/*	char fname[80];
	params = extractString(params, fname);

	if (!edgeMatchGF) {
		cout << "WARNING: edgeMatchGF not initialized; aborting" << endl;
		return;
	}

	CVExample newEx;
	newEx.init(uMesh->numPts());
	int i;
	edgeMatchGF->applyDef(edgeMatchGF->vars);
	for (i=0; i < uMesh->numPts(); i++) {
		newEx.points[i] = edgeMatchGF->curMesh->getPt(i);
		newEx.conf[i] = edgeMatchGF->conf[i];
	}
	newEx.save(fname);*/
}

void uEdgeMatchColor(const char *params) {
	int i;
	if (!edgeMatchGF)
		return;
	for (i=0; i < uMesh->numPts(); i++) {
		double v = 0.5+0.5*(edgeMatchGF->vars[i*edgeMatchGF->varsPerVert+6]);
		uMesh->getPtColor(i) = Vec3d(1,v,v);
	}
	redrawV();
}

void uSpandex(const char *params) {
	params = extractDouble(params, &spandex);
}

void uZeroPose(const char *params) {
	int tr;
	for (tr = 0; tr < skin.skel->transforms.size(); tr++) {
		SkelTransform *st = skin.skel->transforms.getT(tr);
		if (!st->isIntrinsic) {
			st->zero();
		}
	}
	updateSkel(NULL);
}

void uAnimPoses(const char *params) {
	vector<int> poses;
	vector<int> frames;
	int i, j, mode;
	char fname[256];
	bool saveImages = false;

//	bool updateLocal = true;
//	params = extractBool(params, &updateLocal);

/*	uSolver.updateWeights();
	uSolver.updateSkel(0);
	uSolver.updatePoints(0);
	skin.updateMats();
	skin.updateJoints();
*/

	while (1) {
		i = -1;
		params = extractInt(params, &i);
		if (i < 0)
			break;
		poses.push_back(i);
		i = -1;
		params = extractInt(params, &i);
		if (i < 0)
			break;
		frames.push_back(i);
	}
	mode = i;
	if (mode < -1)
		saveImages = true;

	Skeleton *firstSkel = new Skeleton();
	Skeleton *lastSkel = new Skeleton();
	firstSkel->copyFrom(skin.skel);
	lastSkel->copyFrom(skin.skel);

	int frame = 0;
	uSolver.updateSkel(poses[0]);
	lastSkel->copyVals(skin.skel, COPY_POSE);
	for (i=0; i < frames.size(); i++) {
//		cout << "animating " << frames[i] << " frames from pose " << poses[i] << " to " << poses[i+1] << endl;
		firstSkel->copyVals(lastSkel, COPY_POSE);
		uSolver.updateSkel(poses[i+1]);
		lastSkel->copyVals(skin.skel, COPY_POSE);

		for (j=0; j < frames[i]; j++) {
			double interp = 1.0 * j / (frames[i]);

			if (mode == -3) {
				int firstCh = frame / 70;
				double w = (frame - firstCh * 70) / 70.0;
				int k;
				for (k = 0; k < uSolver.numComponents; k++)
					curComps[k] = dataSet.charMu[firstCh][k] * (1.0-w) + 
						dataSet.charMu[firstCh+1][k] * w;
				uShowComps();
			}

			skin.skel->interpVals(firstSkel, lastSkel, interp, true);
			updateSkel(NULL);

/*			double weights[20];
			int joint = skin.skel->transforms.lookupName("lShoulderQ");
			QuatNorm q = skin.skel->transforms.getT(joint)->curCoord.q;
			int numSamp = skin.transPddOfs[joint+1]-skin.transPddOfs[joint];
			Vec4d q2 =  *(skin.pddQuats + skin.transPddOfs[joint]);
			cout << q2 << endl;
			cout << q << endl;
			cout << quatDist(q, QuatNorm(q2[0], q2[1], q2[3], q2[3])) <<
				" " << quatDist(q, QuatNorm(0,0,0,1)) << endl;
			knnQuatInterp(q, skin.pddQuats + skin.transPddOfs[joint], 
				numSamp, skin.pddKeys, weights);
			int k;
			for (k=0; k < numSamp; k++) cout << weights[k] << " ";
			cout << endl;
*/
			redrawVNow();

			if (saveImages) {
				sprintf(fname, "anim/frame%04d.tga", frame);
				cout << "saving " << fname << endl;
				uiScreenshot(fname);
			}
			uiWait();
			frame++;
		}
	}
}

void uRepose(const char *params) {
	int newPose = 0;
	params = extractInt(params, &newPose);

	Skeleton *tempSkel = new Skeleton();
	tempSkel->copyFrom(skin.skel);
	uSolver.updateSkel(newPose);
	tempSkel->copyVals(skin.skel, COPY_POSE);
	skin.skel->copyVals(tempSkel);

	updateSkel(NULL);
	redrawV();

	delete tempSkel;
}

void uAnimSeq(const char *params) {
	bool save = false;
	int startF = 0, endF = dataSet.numExamples;
	params = extractBool(params, &save);
	params = extractInt(params, &startF);
	params = extractInt(params, &endF);
	int f;

	Skeleton *origSkel = new Skeleton();
	origSkel->copyFrom(skin.skel);

	for (f=startF; f < endF; f++) {
		uShow(f);

/*		uSolver.updateSkel(f);
		origSkel->copyVals(skin.skel, Skeleton::COPY_INT);
		skin.skel->copyVals(origSkel);
		((SkelTranslation*)skin.skel->transforms.getT("baseT"))->curVal = Vec3d();
		updateSkel(NULL);
*/
		redrawVNow();
		if (save) {
			char fname[80];
			sprintf(fname, "anim/frame%04d.tga", f);
			uiScreenshot(fname);
			uiWait();
		}
	}
}

void uMixChar(const char *params) {
	int i, ch;
	double weight;

	memset(curComps, 0, sizeof(double)*NUM_CUR_COMPS);
	curComps[0] = 1;

	while (1) {
		ch = -1;
		params = extractInt(params, &ch);
		if (ch < 0)
			break;
		params = extractDouble(params, &weight);
		
		for (i=1; i < uSolver.numComponents; i++)
			curComps[i] += weight * dataSet.charMu[ch][i];
	}

	uShowComps();
}

void uSaveSkinPts(const char *params) {
	int ch = -1;
	char fname[80], fmask[80];
	int ex, pt;
	FILE *f;
	char version = '0';
	int size = uMesh->numPts() * 3;

	params = extractInt(params, &ch);
	params = extractString(params, fmask);

	int minCh = 0, maxCh = dataSet.numCharacters - 1;
	if (ch >= 0) {
		minCh = ch; 
		maxCh = ch;
	}

	for (ch = minCh; ch <= maxCh; ch++) {
		ex = dataSet.charIndex[ch];
		uShow(ex);
		sprintf(fname, "%d", ex);
		uDressFromEx(fname);
		ex++;

		for (; (ex < dataSet.numExamples) && (dataSet.examples[ex].character == ch); ex++) {
			uShow(ex);
			redrawVNow();
			uiWait();

			sprintf(fname, fmask, dataSet.examples[ex].charName, dataSet.examples[ex].poseCh);
			cout << "saving " << fname << endl;
			if (!openFile(&f, fname, "wb", "points file"))
				continue;

			fwrite(&version, sizeof(char), 1, f);
			fwrite(&size, sizeof(int), 1, f);
			for (pt=0; pt < uMesh->numPts(); pt++)
				fwrite(uMesh->getPt(pt).n, sizeof(double), 3, f);

			fclose(f);
		}
	}

	/* old uSaveDressMat
	int ex, var;

	ofstream out;
	char fname[256];
	strcpy(fname, "dressmat.txt");
	params = extractString(params, fname, 256);

	if (!openOFStream(&out, fname, "dress matrix"))
		return;

	for (ex=0; ex < dataSet.numExamples; ex++) {
		if (!dataSet.examples[ex].isDress)
			continue;

		if (dataSet.examples[ex].numPts < 100) {
			// marker set: only save certain points
			for (var = 1; var < 73; var++) {
				if (var == 16 || var == 19)
					continue;	// skip 10th rib markers
				Vec3d exPt;
				dataSet.examples[ex].getPt(var, &exPt);
				out << exPt << " ";
			}
		}
		else {
			for (var = 0; var < dataSet.examples[ex].numPts; var++) {
				Vec3d exPt;
				dataSet.examples[ex].getPt(var, &exPt);
				out << exPt << " ";
			}
		}
		out << endl;
	}

	out.close();*/
}

void uLoadSkinPts(const char *params) {
	char fname[80];
	FILE *f;
	char version;
	int size = uMesh->numPts() * 3;
	int i, pt;

	params = extractString(params, fname, 80);

	if (!openFile(&f, fname, "rb", "points file"))
		return;

	fread(&version, sizeof(char), 1, f);
	fread(&i, sizeof(int), 1, f);
	if (i != size) {
		cout << "point size mismatch" << endl;
		return;
	}
	for (pt=0; pt < uMesh->numPts(); pt++)
		fread(uMesh->getPt(pt).n, sizeof(double), 3, f);

	uMesh->calcNormals();
	redrawV();

	fclose(f);
}

/*
void uSaveTemplates(const char *params) {
	char mask[256], fname[256], s[256];
	strcpy(mask, "../ganger/tia/disp/%d%c.disp");
	params = extractString(params, mask, 256);

	int ch, ex;
	for (ch = 0; ch < dataSet.numCharacters; ch++) {
		ex = dataSet.charIndex[ch];
		uShow(ex);
		sprintf(s, "%d", ex);
		uDressFromEx(s);
		
		
	}
}*/

void uSaveMesh(const char *params) {
	char fname[256];
	params = extractString(params, fname, 256);
	uMesh->saveFile(fname);
}


void uLoadMu(const char *params) {
	char fname[256], fname2[256];
	ifstream in;
	int ch, comp, comp2;

	fname2[0] = 0;
	params = extractString(params, fname, 256);
	params = extractString(params, fname2, 256);
	if (!openIFStream(&in, fname, "mu"))
		return;

	for (comp=0; comp < uSolver.numComponents; comp++) {
		for (ch=0; ch < dataSet.numCharacters; ch++) {
			if (comp == 0)
				dataSet.charMu[ch][comp] = 1;
			else
				in >> dataSet.charMu[ch][comp];
		}
	}

	in.close();

	if (strlen(fname2) > 0 && openIFStream(&in, fname2, "mu variance")) {
		for (comp = 1; comp < uSolver.numComponents; comp++) {
			double scale;
			in >> scale;
			scale = sqrt(scale);

			for (ch=0; ch < dataSet.numCharacters; ch++) {
				dataSet.charMu[ch][comp] /= scale;
			}
		}
	}

/*	// fill in phi
	for (ch = 0; ch < dataSet.numCharacters; ch++) {
		dataSet.charPhi[ch].MakeZero();
		
		for (comp = 0; comp < uSolver.numComponents; comp++) {
			for (comp2 = 0; comp2 < uSolver.numComponents; comp2++) {
				dataSet.charPhi[ch][comp][comp2] = dataSet.charMu[ch][comp] * dataSet.charMu[ch][comp2];
			}
		}
	}*/
	dataSet.phiLogEntropy = 0;		// actually, it's -infinity
}

void uSaveMu(const char *params) {
	char fname[256];
	ofstream out;
	int ch, comp, comp2;

	params = extractString(params, fname, 256);
	if (!openOFStream(&out, fname, "mu"))
		return;

	for (comp=1; comp < uSolver.numComponents; comp++) {
		for (ch=0; ch < dataSet.numCharacters; ch++) {
			out << dataSet.charMu[ch][comp] << " ";
		}
		out << endl;
	}

	out.close();
}

void uDumpMu(const char *params) {
	int first = 0;
	int last = 5;
	params = extractInt(params, &first);
	params = extractInt(params, &last);

	int i, j;
	for (i=first; i < min(last, dataSet.numCharacters); i++) {
		cout << "character " << i << ": ";
		for (j=0; j < uSolver.numComponents; j++) {
			cout << dataSet.charMu[i][j] << " ";
		}
		cout << endl;
	}
}

void uPerCharMu(const char *params) {
	// set reconstruction weights so that each character is independent from the rest
	int ch;

	for (ch = 0; ch < min(uSolver.numComponents, dataSet.numCharacters); ch++) {
		dataSet.charMu[ch] = vl_zero;
		dataSet.charMu[ch][ch] = 1;
	}
}

void uFreeEnergy(const char *params) {
/*	double a = uSolver.evaluateFunctionExp() / sqr(uSolver.vertSigma);
	double b = 0;
	double c = 
		uSolver.numVert * log(uSolver.vertSigma) + 
		uSolver.numWeightNeigh * log(uSolver.weightNeighSigma) + 
		uSolver.numPddNeigh * log(uSolver.pddNeighSigma) +
		uSolver.numDressMuNeigh * log(uSolver.dressMuNeighSigma) + 
		uSolver.numDressWNeigh * log(uSolver.dressWNeighSigma);
	double d = dataSet.phiLogEntropy;
	int ch, comp;
	for (ch=0; ch < dataSet.numCharacters; ch++) {
		for (comp=1; comp < uSolver.numComponents; comp++) {
			b += sqr(dataSet.charMu[ch][comp]);
		}
	}
	cout << "free energy: " << (a+b+c-d) << " (" << a << " + " << b << " + " << c << " - " << d << ")" << endl;*/
}

void uEStep(const char *params) {
	uSolver.eStep();
	uShow();
}

void uTest(const char *params) {
	// check derivatives
	int i;
/*	SkelCombinedTransform *shoulder = (SkelCombinedTransform*)skin.skel->transforms.getT("lShoulderCR");
	SkelEulerRotation *euler[3];
	euler[0] = (SkelEulerRotation*)shoulder->orig[0];
	euler[1] = (SkelEulerRotation*)shoulder->orig[1];
	euler[2] = (SkelEulerRotation*)shoulder->orig[2];

	euler[0]->curAngle = boundedRand(-3.14, 3.14);
	euler[1]->curAngle = boundedRand(-3.14, 3.14);
	euler[2]->curAngle = boundedRand(-3.14, 3.14);
	skin.skel->updateCoords();
	QuatNorm orig = shoulder->curCoord.q;

	double cx = cos(euler[0]->curAngle/2);
	double sx = sin(euler[0]->curAngle/2);
	double cy = cos(euler[2]->curAngle/2);
	double sy = sin(euler[2]->curAngle/2);
	double cz = cos(euler[1]->curAngle/2);
	double sz = sin(euler[1]->curAngle/2);
	double delta = 0.00001;

	cout << 
		( 0.5*cy*cz*cx + 0.5*sy*sz*sx) << " " <<
		(-0.5*sy*cz*sx - 0.5*cy*sz*cx) << " " <<
		(-0.5*cy*sz*sx + 0.5*sy*cz*cx) << " " <<
		( 0.5*cy*cz*sx - 0.5*sy*sz*cx) << " " << endl;
	euler[0]->curAngle += delta;
	skin.skel->updateCoords();
	for (i=0; i < 4; i++) {
		cout << (shoulder->curCoord.q[i] - orig[i])/delta << " ";
	}
	cout << endl;
	euler[0]->curAngle -= delta;

	cout << 
		(-0.5*cy*sz*sx - 0.5*sy*cz*cx) << " " <<
		(-0.5*sy*sz*cx - 0.5*cy*cz*sx) << " " <<
		( 0.5*cy*cz*cx - 0.5*sy*sz*sx) << " " <<
		( 0.5*cy*sz*cx - 0.5*sy*cz*sx) << " " << endl;
	euler[1]->curAngle += delta;
	skin.skel->updateCoords();
	for (i=0; i < 4; i++) {
		cout << (shoulder->curCoord.q[i] - orig[i])/delta << " ";
	}
	cout << endl;
	euler[1]->curAngle -= delta;

	cout << 
		(-0.5*sy*cz*sx - 0.5*cy*sz*cx) << " " <<
		( 0.5*cy*cz*cx + 0.5*sy*sz*sx) << " " <<
		(-0.5*sy*sz*cx + 0.5*cy*cz*sx) << " " <<
		(+0.5*sy*cz*cx - 0.5*cy*sz*sx) << " " << endl;
	euler[2]->curAngle += delta;
	skin.skel->updateCoords();
	for (i=0; i < 4; i++) {
		cout << (shoulder->curCoord.q[i] - orig[i])/delta << " ";
	}
	cout << endl;
	euler[2]->curAngle -= delta;
*/

	SkelCombinedTransform *clav = (SkelCombinedTransform*)skin.skel->transforms.getT("lClavicleCQ");
	SkelEulerRotation *euler[2];
	euler[0] = (SkelEulerRotation*)clav->orig[0];
	euler[1] = (SkelEulerRotation*)clav->orig[1];

	euler[0]->curAngle = boundedRand(-3.14, 3.14);
	euler[1]->curAngle = boundedRand(-3.14, 3.14);
	skin.skel->updateCoords();
	QuatNorm orig = clav->curCoord.q;

	double cx = cos(euler[1]->curAngle/2);
	double sx = sin(euler[1]->curAngle/2);
	double cz = cos(euler[0]->curAngle/2);
	double sz = sin(euler[0]->curAngle/2);
	double delta = 0.00001;

	cout << 
		( 0.5*sx*sz) << " " <<
		(-0.5*sx*cz) << " " <<
		(-0.5*cx*cz) << " " <<
		(-0.5*cx*sz) << " " << endl;
	euler[0]->curAngle += delta;
	skin.skel->updateCoords();
	for (i=0; i < 4; i++) {
		cout << (clav->curCoord.q[i] - orig[i])/delta << " ";
	}
	cout << endl;
	euler[0]->curAngle -= delta;

	cout << 
		(-0.5*cx*cz) << " " <<
		(-0.5*cx*sz) << " " <<
		( 0.5*sx*sz) << " " <<
		(-0.5*sx*cz) << " " << endl;
	euler[1]->curAngle += delta;
	skin.skel->updateCoords();
	for (i=0; i < 4; i++) {
		cout << (clav->curCoord.q[i] - orig[i])/delta << " ";
	}
	cout << endl;
	euler[1]->curAngle -= delta;


/*	int pt, tri, tInd;
	int newNumOrigPts = uSolver.numOrigPts;
	int preSubPts = uMesh->numPts();
	int *mappingVerts;
	double *mappingWeights;

	// subdivide the mesh
    butterflySubdiv(uMesh, &mappingVerts, &mappingWeights);

	// initialize re-ordering
	int *vertMap = new int[uMesh->numPts()];
	int *invMap = new int[uMesh->numPts()];
	for (pt=0; pt < uMesh->numPts(); pt++) {
		vertMap[pt] = pt;
	}

	// re-order points so that non-mirrored are first
	for (pt=preSubPts; pt < uMesh->numPts(); pt++) {
		if (uMesh->getPt(pt)[0] < 1e-5) {
			swap(vertMap[newNumOrigPts], vertMap[pt]);
			newNumOrigPts++;
		}
	}
	for (pt=0; pt < uMesh->numPts(); pt++) {
		invMap[vertMap[pt]] = pt;
	}
	cout << "before reordering: " << uSolver.numOrigPts << 
		"; after reordering: " << newNumOrigPts << endl;

	// save low-high vertex mapping
	FILE *fMap;
	if (openFile(&fMap, "data/low-hi.map", "wb", "mapping file")) {
		cout << "saving mapping" << endl;
		// save mapping header
//		char version = '0';
//		fwrite(&version, sizeof(char), 1, mapF);
		pt = uSolver.numPts;
		fwrite(&pt, sizeof(int), 1, fMap);
		pt = uMesh->numPts();
		fwrite(&pt, sizeof(int), 1, fMap);

		// save mapping data
		for (pt=0; pt < uMesh->numPts(); pt++) {
			int src;
			if (vertMap[pt] < uSolver.numPts) {
				// original points stay in place
				src = 1;
				fwrite(&src, sizeof(int), 1, fMap);
				fwrite(&vertMap[pt], sizeof(int), 1, fMap);
				float w = 1.0;
				fwrite(&w, sizeof(float), 1, fMap);
			}
			else {
				src = 8;
				fwrite(&src, sizeof(int), 1, fMap);
				for (src=0; src < 8; src++) {
					fwrite(&mappingVerts[8*(vertMap[pt] - uSolver.numPts) + src], sizeof(int), 1, fMap);
					float w = mappingWeights[8*(vertMap[pt] - uSolver.numPts) + src];
					fwrite(&w, sizeof(float), 1, fMap);
				}
			}
		}
		fclose(fMap);
	}

	// apply the reordering to the trimesh
	Vec3d *origVerts = new Vec3d[uMesh->numPts()];
	for (pt=0; pt < uMesh->numPts(); pt++) {
		origVerts[pt] = uMesh->getPt(pt);
	}
	for (pt=0; pt < uMesh->numPts(); pt++) {
		uMesh->getPt(pt) = origVerts[vertMap[pt]];
	}
	delete []origVerts;
	for (tri=0; tri < uMesh->numTris(); tri++) {
		for (tInd = 0; tInd < 3; tInd++) {
			uMesh->getTri(tri, tInd) = invMap[uMesh->getTri(tri, tInd)];
		}
	}
	uMesh->savePly("data/subdiv.ply");

	FILE *f;
	if (openFile(&f, "data/mirror-map.dat", "wb", "mirror map")) {
		cout << "saving mirror map" << endl;
		int i;
		i = 0;
		fwrite(&i, sizeof(int), 1, f);	// version
		i = newNumOrigPts;
		fwrite(&i, sizeof(int), 1, f);
		i = uMesh->numPts();
		fwrite(&i, sizeof(int), 1, f);
		i = 0; //numOrigTris;
		fwrite(&i, sizeof(int), 1, f);
		i = 0; //meshes[ind]->numTris();
		fwrite(&i, sizeof(int), 1, f);

		for (pt=0; pt < uMesh->numPts(); pt++) {
			Vec3d v = uMesh->getPt(pt);
			v[0] = -v[0];

			if (fabs(v[0]) < 1e-5) {
				fwrite(&pt, sizeof(int), 1, f);
			}
			else {
				int pt2;
				for (pt2=0; pt2 < uMesh->numPts(); pt2++) {
					if ((uMesh->getPt(pt2) - v).length2() < 1e-9)
						break;
				}
				if (pt2 >= uMesh->numPts()) {
					cout << "can't find mirror point for vertex " << pt << ": " << v << endl;
				}
				fwrite(&pt2, sizeof(int), 1, f);
			}
		}
		//fwrite(triMirrorMap, sizeof(int), meshes[ind]->numTris(), f);
		fclose(f);
	}

	// everything's messed up now; let's exit
	exit(0);

	redrawV();*/
}

void smInitMesh(const char *params) {
	smInitMesh(uMesh);
}

void smInitTarget(const char *params) {
	smInitTarget(uMesh);
}

void smRun(const char *params) {
	int maxIter = -1;
	params = extractInt(params, &maxIter);
	smRun(maxIter);
}

void updatePDDW() {
	int pdd, inf;
	double w;

	if (pddWJoint < 0)
		return;


	for (pdd=0; pdd < skin.pdds.size(); pdd++) {
		if (skin.pdds[pdd]->transInd == pddWJoint) {
			w = 0;
			if (skin.pdds[pdd]->rbf->mN > 0)
				w = skin.pdds[pdd]->rbf->curWeights[0];
			uniUI->pdd0VS->value(w);

			w = 0;
			if (skin.pdds[pdd]->rbf->mN > 1)
				w = skin.pdds[pdd]->rbf->curWeights[1];
			uniUI->pdd1VS->value(w);

			w = 0;
			if (skin.pdds[pdd]->rbf->mN > 2)
				w = skin.pdds[pdd]->rbf->curWeights[2];
			uniUI->pdd2VS->value(w);

			w = 0;
			if (skin.pdds[pdd]->rbf->mN > 3)
				w = skin.pdds[pdd]->rbf->curWeights[3];
			uniUI->pdd3VS->value(w);
		}
	}
}

void setPDDWJoint(const char *name) {
	pddWJoint = skin.skel->transforms.lookupName(name);
	updatePDDW();
}

void tmLoad(const char* params) {
	char fname[256];
	fname[0] = 0;
	params = extractString(params, fname, 256);

	delete dispMesh;
	dispMesh = NULL;

	if (strlen(fname) < 1)
		return;
	dispMesh = new TriMesh();
	dispMesh->loadFile(fname);

	int i;
	for (i=0; i < dispMesh->numPts(); i++)
		dispMesh->getPtColor(i) = Vec3d(1, 1, 1);
	dispMesh->calcNormals();

	redrawV();
}

void uRandom(const char *params) {
	int comp;
	for (comp=1; comp <uSolver.numComponents; comp++) {
		curComps[comp] = sqrt(-2 * log(boundedRand(0, 1) + 0.0000001)) * cos(boundedRand(0, PI*2)) * 0.6;
	}
	uShowComps(true);
}

void uLoadSliders(const char *params) {
	char fname[256];
	params = extractString(params, fname, 256);

	initFeatures(uSolver.numComponents - 1);
	loadFeatures(fname);
}

void uSetSliders(const char *params) {
	float fVals[10];
	double d;
	int i;

	for (i=0; i < 10; i++) {
		params = extractDouble(params, &d);
		fVals[i] = d;
	}
	setFeatures(fVals, curComps);
	uShowComps(true);
	redrawV();
}

void uSetSliderBase(const char *params) {
	float fVals[10];
	double d;
	int i;

	for (i=0; i < 10; i++) {
		params = extractDouble(params, &d);
		fVals[i] = d;
	}
	setFeatureBase(fVals, curComps);
}

void uAnimSliders(const char *params) {
	float fVals[10], f2Vals[10], iVals[10];
	double d;
	int i, frame, numS, numFrames, special = 0;
	char fName[256], mask[256];

	params = extractString(params, mask, 256);
	params = extractInt(params, &numFrames);
	params = extractInt(params, &numS);
	if (mask[0] = '=') {
		strcpy(mask, "anim/frame%04d.tga");
	}

	for (i=0; i < numS; i++) {
		params = extractDouble(params, &d);
		fVals[i] = d;
	}
	for (i=0; i < numS; i++) {
		params = extractDouble(params, &d);
		f2Vals[i] = d;
	}
	params = extractInt(params, &special);


	for (frame=0; frame < numFrames; frame++) {
		double w = 1.0 * frame / (numFrames-1);
		for (i=0; i < numS; i++) {
			iVals[i] = (1.0 - w) * fVals[i] + w  *f2Vals[i];
		}
		setFeatures(iVals, curComps);
		uShowComps(true);
		if (special == 1)
			uZapPdds("");
		else if (special == 2)
			uCopyPdds("0");

		redrawVNow();
		uiWait();
		if (strlen(mask) > 1) {
			sprintf(fName, mask, frame);
			cout << "saving " << fName << endl;
			uiScreenshot(fName);
		}
	}
}

void uAnimElbow(const char *params) {
	int i, frame, numFrames = 20, special = 0;
	char fName[256], mask[256];

	params = extractString(params, mask, 256);
	params = extractInt(params, &numFrames);
	if (mask[0] = '=') {
		strcpy(mask, "anim/frame%04d.tga");
	}
	params = extractInt(params, &special);

	for (frame=0; frame < numFrames; frame++) {
		double w = 1.0 * frame / (numFrames-1);
		SkelEulerRotation *trans = (SkelEulerRotation*)skin.skel->transforms.getT("lElbowA");
		trans->curAngle = -w * 1.88;

		updateSkel(NULL);
		if (special == 1)
			uZapPdds("");
		else if (special == 2)
			uCopyPdds("0");

		redrawVNow();
		uiWait();
		if (strlen(mask) > 1) {
			sprintf(fName, mask, frame);
			cout << "saving " << fName << endl;
			uiScreenshot(fName);
		}
	}
}

void initUMaster() {
/*	RBF testR(3, 4, 0.4, true);
	QuatNorm q, q1(0.5,0,0),q2(1.0,0,0);
	double d[4],weights[3];
	double quats[8] = {q1.x,q1.y,q1.z,q1.w,q2.x,q2.y,q2.z,q2.w};
	testR.init(quats,true);

	int i, j;
	for (i=0; i < 11; i++) {
		q = QuatNorm(i*1.0/10.0,0,0);
		d[0] = q.x; d[1] = q.y; d[2] = q.z; d[3] = q.w;
		testR.eval(d, weights);
		for (j=0; j < 3; j++)
			cout << weights[j] << " ";
		cout << endl;
	}
*/

	registerFunction(uLoadExampleSet, "uLoadExampleSet");
	registerFunction(uShow, "uShow");
	registerFunction(uSetPose, "uSetPose");
	registerFunction(uColor, "uColor");
	registerFunction(uColorPdd, "uColorPdd");
	registerFunction(uColorPdn, "uColorPdn");
	registerFunction(uColorMarker, "uColorMarker");
	registerFunction(uStartSolver, "uStartSolver");
	registerFunction(uStopSolver, "uStopSolver");
	registerFunction(uRenormalizeX, "uRenormalizeX");
	registerFunction(uAutoWeight, "uAutoWeight");
	registerFunction(uSetMaxSolveEx, "uSetMaxSolveEx");
	registerFunction(uSetSingleComponent, "uSetSingleComponent");
	registerFunction(uNormalizeDofs, "uNormalizeDofs");
	registerFunction(uNormalizeNM, "uNormalizeNM");
	registerFunction(uWeightsFromTex, "uWeightsFromTex");
	registerFunction(uColorFromPly, "uColorFromPly");
	registerFunction(uSaveExampleMat, "uSaveExampleMat");
	registerFunction(uSaveMarkerSprings, "uSaveMarkerSprings");
	registerFunction(uMakeMarkerSprings, "uMakeMarkerSprings");
	registerFunction(uTestRBFs, "uTestRBFs");
	registerFunction(uLoadTex, "uLoadTex");
	registerFunction(uTangentSpaceVis, "uTangentSpaceVis");
	registerFunction(uBuildNormalMap, "uBuildNormalMap");
	registerFunction(uBuildNormalMapEx, "uBuildNormalMapEx");

	registerFunction(uSaveSolve, "uSaveSolve");
	registerFunction(uLoadSolve, "uLoadSolve");
	registerFunction(uSavePoses, "uSavePoses");
	registerFunction(uLoadPoses, "uLoadPoses");
	registerFunction(uSavePosesEx, "uSavePosesEx");
	registerFunction(uLoadPosesEx, "uLoadPosesEx");
	registerFunction(uSaveWeights, "uSaveWeights");
	registerFunction(uLoadWeights, "uLoadWeights");
	registerFunction(uSavePdds, "uSavePdds");
	registerFunction(uLoadPdds, "uLoadPdds");
	registerFunction(uSaveDress, "uSaveDress");
	registerFunction(uLoadDress, "uLoadDress");
	registerFunction(uBuildMapping, "uBuildMapping");
	registerFunction(uDressFromEx, "uDressFromEx");
	registerFunction(uSaveSkinDisp, "uSaveSkinDisp");
	registerFunction(uZeroVars, "uZeroVars");
	registerFunction(uSkelFromSurf, "uSkelFromSurf");
	registerFunction(uSaveSkelMat, "uSaveSkelMat");
	registerFunction(uLoadPCASkels, "uLoadPCASkels");
	registerFunction(uLoadPoseFiles, "uLoadPoseFiles");
	registerFunction(uPoseToEx, "uPoseToEx");
	registerFunction(uFixMirrorLine, "uFixMirrorLine");
	registerFunction(uAdjustWeights, "uAdjustWeights");

	registerFunction(uLoadPddDefs, "uLoadPddDefs");

	registerFunction(uSaveCurPose, "uSaveCurPose");
	registerFunction(uLoadCurPose, "uLoadCurPose");
	registerFunction(uSaveAllPoses, "uSaveAllPoses");
	registerFunction(uShowDress, "uShowDress");
	registerFunction(uDumpWeights, "uDumpWeights");
	registerFunction(uDumpStats, "uDumpStats");
	registerFunction(uDumpVars, "uDumpVars");
	registerFunction(uDumpPdds, "uDumpPdds");
	registerFunction(uDumpPddWeights, "uDumpPddWeights");
	registerFunction(uDumpPtInfo, "uDumpPtInfo");
	registerFunction(uDumpJoints, "uDumpJoints");
	registerFunction(uDumpExJoints, "uDumpExJoints");
	registerFunction(uDumpJointDofs, "uDumpJointDofs");
	registerFunction(uDumpJointEx, "uDumpJointEx");
	registerFunction(uDumpMaps, "uDumpMaps");
	registerFunction(uZapPdds, "uZapPdds");
	registerFunction(uCopyPdds, "uCopyPdds");
	registerFunction(uLoadPddMask, "uLoadPddMask");
	registerFunction(uGreedyPdds, "uGreedyPdds");
//	registerFunction(uFixShoulder, "uFixShoulder");
	registerFunction(uLoadMocap, "uLoadMocap");
	registerFunction(uFixMocap, "uFixMocap");
	registerFunction(uLoadFix, "uLoadFix");
	registerFunction(uSaveFix, "uSaveFix");
	registerFunction(uSaveMocap, "uSaveMocap");
	registerFunction(uExToMocap, "uExToMocap");

	registerFunction(uSaveCurPDs, "uSaveCurPDs");
	registerFunction(uLoadCurPDs, "uLoadCurPDs");

	registerFunction(uMirrorEx, "uMirrorEx");

	registerFunction(uMatchAll, "uMatchAll");
	registerFunction(uAutoNW, "uAutoNW");
	registerFunction(uEdgeMatchInit, "uEdgeMatchInit");
	registerFunction(uEdgeMatchLoadNW, "uEdgeMatchLoadNW");
	registerFunction(uEdgeMatchLoadSW, "uEdgeMatchLoadSW");
	registerFunction(uEdgeMatchLoadRestrict, "uEdgeMatchLoadRestrict");
	registerFunction(uEdgeMatchMap, "uEdgeMatchMap");
	registerFunction(uEdgeMatchNorm, "uEdgeMatchNorm");
	registerFunction(uEdgeMatchSolve, "uEdgeMatchSolve");
	registerFunction(uEdgeMatchStopSolve, "uEdgeMatchStopSolve");
	registerFunction(uEdgeMatchSaveEx, "uEdgeMatchSaveEx");
	registerFunction(uEdgeMatchColor, "uEdgeMatchColor");
	registerFunction(uSpandex, "uSpandex");

	registerFunction(uZeroPose, "uZeroPose");
	registerFunction(uAnimPoses, "uAnimPoses");
	registerFunction(uRepose, "uRepose");
	registerFunction(uAnimSeq, "uAnimSeq");
	registerFunction(uMixChar, "uMixChar");

	registerFunction(uTest, "uTest");
	registerFunction(uRandom, "uRandom");

	registerFunction(uSaveSkinPts, "uSaveSkinPts");
	registerFunction(uLoadSkinPts, "uLoadSkinPts");
//	registerFunction(uSaveTemplates, "uSaveTemplates");
	registerFunction(uSaveMesh, "uSaveMesh");
	registerFunction(uLoadMu, "uLoadMu");
	registerFunction(uSaveMu, "uSaveMu");
	registerFunction(uDumpMu, "uDumpMu");
	registerFunction(uPerCharMu, "uPerCharMu");
	registerFunction(uFreeEnergy, "uFreeEnergy");
	registerFunction(uEStep, "uEStep");

	registerFunction(uLoadSliders, "uLoadSliders");
	registerFunction(uSetSliders, "uSetSliders");
	registerFunction(uSetSliderBase, "uSetSliderBase");
	registerFunction(uAnimSliders, "uAnimSliders");

	registerFunction(uAnimElbow, "uAnimElbow");

	registerFunction(smInitMesh, "smInitMesh");
	registerFunction(smInitTarget, "smInitTarget");
	registerFunction(smRun, "smRun");

	registerFunction(tmLoad, "tmLoad");
}