#include "doppel2.h"
#include "uSolver.h"
#include "uSkin.h"
#include "uMaster.h"
#include "skeleton.h"
#include "knnQuat.h"
#include "trimesh.h"
#include "edgelist.h"
#include "markers.h"
#include "normalMap.h"
#include "trimesh_util.h"
#include "vl/VLd.h"

extern int visSkel;

const double NORMAL_TOL = cos(20 * DEG_TO_RAD);

const double NUM_EX = 85;
//const double NUM_EX = 215.0;
//const double NUM_EX = 18.0;
const double POSE_DOF_SCALE = NUM_EX;
const double WEIGHT_DOF_SCALE = 1; //NUM_EX;

#define ICM

// UExample ===============================================

UExample::UExample() {
	character = 0;
	numPts = 0;
	minPose = false;
	ptsConf = NULL;
	pts = NULL;
	mesh = NULL;
	edgeList = NULL;
	vertList = NULL;
	lookup = NULL;
	normals = new Vec3d[NORMAL_MAP_SIZE];
}

void UExample::init(int iNumPts) {
	numPts = iNumPts;
	ptsConf = new double[numPts];
	memset(ptsConf, 0, sizeof(double)*numPts);
	pts = new Vec3d[numPts];
}

void UExample::init(const char *meshName) {
	mesh = new TriMesh;
	mesh->loadFile(fname);
	mesh->calcHBB(16);
	mesh->calcNormals();
	int i;
	for (i=0; i < mesh->numPts(); i++)
		mesh->getPtColor(i) = Vec3d(0.8, 0.8, 0.8);
	
	edgeList = new EdgeList();
	edgeList->buildFromTriMesh(*mesh);
	vertList = new char[edgeList->numVerts];
	edgeList->markVerts(vertList);
}

void UExample::initLookup(int *lu, int lus) {
	lookup = lu;
	luSize = lus;
	luOffsets = new Vec3d[numPts];
}

bool UExample::save(const char *fname) {
	FILE *f;

	if (openFile(&f, fname, "wb", "example file")) {
		fwrite(&numPts, sizeof(int), 1, f);
		fwrite(pts, sizeof(Vec3d), numPts, f);
		fwrite(ptsConf, sizeof(double), numPts, f);
		fwrite(normals, sizeof(Vec3d), NORMAL_MAP_SIZE, f);
		fclose(f);
		return true;
	}
	else
		return false;
}

bool UExample::load(const char *fname) {
	FILE *f;
	int j;

	if (openFile(&f, fname, "rb", "example file")) {
		int numPoints;
		fread(&numPoints, sizeof(int), 1, f);
		init(numPoints);
		fread(pts, sizeof(Vec3d), numPoints, f);
		fread(ptsConf, sizeof(double), numPoints, f);
		fread(normals, sizeof(Vec3d), NORMAL_MAP_SIZE, f);
		for (j=0; j < numPts; j++) {
			ptsConf[j] = min(1.0, /*0.1+*/ ptsConf[j]);
			if (!_finite(ptsConf[j]))
				ptsConf[j] = 0;
			if (ptsConf[j] > 0)
				ptsConf[j] = 0.5 + 0.5 * ptsConf[j];
		}
		fclose(f);
		return true;
	}
	else {
		init(0);
		return false;
	}
}

bool UExample::getPt(int ind, Vec3d *v, double *conf) {
	if (lookup) {
		if (ind >= luSize || ind < 0 || lookup[ind] < 0) {
			if (conf) *conf = 0;
			return false;
		}
		*v = pts[lookup[ind]] + luOffsets[lookup[ind]];
		if (conf) {
			if (ptsConf)
				*conf = ptsConf[lookup[ind]];
			else
				*conf = 1;
		}
		return true;
	}
	else {
		*v = pts[ind];
		if (conf) {
			if (ptsConf)
				*conf = ptsConf[ind];
			else
				*conf = 1;
		}
		return true;
	}
}


// UDataSet ===============================================

UDataSet::UDataSet() {
	numExamples = 0;
	numCharacters = 0;
	examples = NULL;
	charIndex = NULL;
	charMu = NULL;
//	charPhi = NULL;
}

void UDataSet::init(int iNumExamples, int iNumCharacters) {
	numExamples = iNumExamples;
	examples = new UExample[numExamples];

	numCharacters = iNumCharacters;
	charIndex = new int[numCharacters];
	memset(charIndex, 0, sizeof(int) * numCharacters);
	charMu = new VLVecd[numCharacters];
//	charPhi = new VLMatd[numCharacters];
}

void UDataSet::initLookup(MarkerSet &mrefs) {
	int i;
	int n = 0;
	for (i=0; i < mrefs.numMarkers; i++)
		if (mrefs.markers[i].kind == mkrPT && mrefs.markers[i].baryVerts[0] >= n)
			n = mrefs.markers[i].baryVerts[0]+1;

	cout << n << endl;
	int *lu = new int[n];
	for (i=0; i < n; i++)
		lu[i] = -1;
	for (i=0; i < mrefs.numMarkers; i++) {
		if (mrefs.markers[i].kind == mkrPT)
			lu[mrefs.markers[i].baryVerts[0]] = i + 3;
	}
	int ex;
	for (ex = 0; ex < numExamples; ex++) {
		examples[ex].initLookup(lu, n);
	}
}


// USolver ================================================

USolver::USolver() {
	dataSet = NULL;
	skin = NULL;
	lastErr = 0;
	optDress = optInt = optPose = optWeight = optPDD = optX = optNM = optNMP = false;
	globalPoseOnly = lockPoseZero = ignorePoints = false;
	maxSolveEx = -1;

	vertSigma = 1000000;
	normSigma = 10000;
	pddSigma = 1000;
	minPoseSigma = 10000;
	nmSigma = 10000;
	markerSpringSigma = 100000000;

	weightNeighSigma = 10000;
	nmNeighSigma = 10;
	pddNeighSigma = 10000.0;
	pdnNeighSigma = 10.0;
	dressMuNeighSigma = 10000;
	dressWNeighSigma = 0; //0.01;

	markerAssignments = NULL;
	pddMask = NULL;
	dfdv = NULL;

	singleComp = -1;
}

void USolver::init(UDataSet *iDataSet, USkin *iSkin, TriMesh *iMesh, const char *mirrorFN, int iNumComponents, int iNumPDComponents) {
	int i, j, index, ch, ex;

	dataSet = iDataSet;
	skin = iSkin;
	uMesh = iMesh;
	numComponents = iNumComponents;
	numPDComponents = iNumPDComponents;

	dfdv = new Vec3d[skin->numPts];

	if (iMesh) {
		// load mirror map
		FILE *f;
		if (openFile(&f, mirrorFN, "rb", "mirror map")) {
			fread(&i, sizeof(int), 1, f);	// version
			if (i != 0) {
				cout << "unknown mirror-map version in " << mirrorFN << endl;
				exit(0);
			}
			fread(&numOrigPts, sizeof(int), 1, f);
			fread(&numPts, sizeof(int), 1, f);
			if (numPts != uMesh->numPts()) {
				cout << "mismatched # of points in mirror-map " << mirrorFN << endl;
				exit(0);
			}
			fread(&i, sizeof(int), 1, f);	// # of original triangles (unused for now)
			fread(&i, sizeof(int), 1, f);	// # of triangles (unused for now)
			mirrorMap = new int[numPts];
			fread(mirrorMap, sizeof(int), numPts, f);
			//fread(triMirrorMap, sizeof(int), meshes[ind]->numTris(), f);
			fclose(f);
		}
		else {
			exit(0);
		}
	}
	else {
		numPts = skin->numPts;
		numOrigPts = numPts;

		ifstream in;
		mirrorMap = new int[numOrigPts];
		if (openIFStream(&in, mirrorFN, "mirror map")) {
			mirrorMap[0] = 0;
			for (i=1; i < numOrigPts; i++) {
				char s[80];
				in >> s;
				in >> mirrorMap[i];
			}
			in.close();
		}
		else {
			for (i=0; i < numOrigPts; i++)
				mirrorMap[i] = -1;
		}
	}

	// find the mirror of each transform
	mirrorTrans = new int[skin->skel->transforms.size()];
	for (i=0; i < skin->skel->transforms.size(); i++) {
		char name[80];
		strcpy(name, skin->skel->transforms.getT(i)->name);
		if (name[0] == 'l' || name[0] == 'r') {
			if (name[0] == 'l')
				name[0] = 'r';
			else
				name[0] = 'l';
			int mInd = skin->skel->transforms.lookupName(name);
			if (mInd > -1)
				mirrorTrans[i] = mInd;
			else
				mirrorTrans[i] = i;
		}
		else
			mirrorTrans[i] = i;
	}

	for (ch=0; ch < dataSet->numCharacters; ch++) {
		dataSet->charMu[ch].SetSize(numComponents);
//		dataSet->charPhi[ch].SetSize(numComponents, numComponents);
	}

	vDress.resize(numOrigPts * 3 * numComponents, true);
    vInt.resize(skin->skel->numIntrinsicDofs * numComponents, true);
	vPose.resize(skin->skel->numPoseDofs * dataSet->numExamples, true); 
	vWeight.resize(numOrigPts * (skin->maxInf-1), true);
	vPDD.resize(skin->numPddPtKeys * 3 * numPDComponents, true);
	if (uMesh)
		vNM = new Vec3d[NORMAL_MAP_SIZE / 2 * numComponents];
	vNMP.resize(skin->numPddNMKeys * 3 * numPDComponents, true);

	index = 0;

	// initialize skeleton vars
	if (dataSet->examples[0].lookup != NULL) {
		// mocap
		ifstream in;
		cout << "loading mocap.po.txt" << endl;
		if (openIFStream(&in, "data/mocap.po.txt", "pose file")) {
			skin->skel->loadPose(in);
			in.close();

			int ip = 0;
			int pp = 0;
			for (i=0; i < skin->skel->transforms.size(); i++) {
				SkelTransform *curTrans = skin->skel->transforms.getT(i);
				if (curTrans->isIntrinsic) {
					curTrans->unloadDofs(vInt.n + ip);
					ip += curTrans->numDofs();
				}
				else {
					int ex;
					for (ex=0; ex < dataSet->numExamples; ex++)
						curTrans->unloadDofs(vPose.n + pp + skin->skel->numPoseDofs * ex);
					pp += curTrans->numDofs();
				}
			}
		}
	}
	else {
		// example shapes
		for (ex = 0; ex < dataSet->numExamples; ex++) {
			ifstream in;
			char poseFN[256];
			sprintf(poseFN, "data/mocap/%c.po.txt", dataSet->examples[ex].poseCh);
			if (openIFStream(&in, poseFN, "pose file")) {
				skin->skel->loadPose(in);
				in.close();
			}
			else {
//				skin->skel->zero();
			}

			int ip = 0;
			int pp = 0;
			for (i=0; i < skin->skel->transforms.size(); i++) {
				SkelTransform *curTrans = skin->skel->transforms.getT(i);
				if (curTrans->isIntrinsic) {
					if (ex == 0) {
						curTrans->unloadDofs(vInt.n + ip);
						ip += curTrans->numDofs();
					}
				}
				else {
					curTrans->unloadDofs(vPose.n + pp + skin->skel->numPoseDofs * ex);
					pp += curTrans->numDofs();
				}
			}
		}
	}

/* old way
	for (ch = 0; ch < dataSet->numCharacters; ch++) {
//		skin->transInit(&dataSet->examples[dataSet->charIndex[ch]]);
		skin->skel->updateCoords();

		int ip = intDofPos;
		int pp = poseDofPos;
		for (i=0; i < skin->skel->transforms.size(); i++) {
			SkelTransform *curTrans = skin->skel->transforms.getT(i);
			if (curTrans->isIntrinsic) {
				if (ch == 0) {
					curTrans->unloadDofs(vars.n + ip);
					ip += curTrans->numDofs();
				}
			}
			else {
				j = dataSet->charIndex[ch];
				while (j < dataSet->numExamples && dataSet->examples[j].character == ch) {
					curTrans->unloadDofs(vars.n + pp + skin->skel->numPoseDofs * j);
					j++;
				}
				pp += curTrans->numDofs();
			}
		}
	}*/

	weightsToVars();

/*	if (pddDofPos >= 0) {
		for (i=pddDofPos; i < numVars; i++)
			vars[i] = boundedRand(-0.001, 0.001);
	}*/

	skin->skel->allocDerivs(skin->skel->numDofs);

	cout << "# of pose dofs: " << skin->skel->numPoseDofs << endl;
}

void USolver::buildNeighborTable(TriMesh *tm) {
	int pt, i, j;
	int *infMap = new int[skin->skel->transforms.size()];
	int *infInd0 = new int[skin->skel->transforms.size()];
	int *infInd1 = new int[skin->skel->transforms.size()];
	NeighborRelation nr;
	int numReg = 0, numZero = 0;

	// build a list of neighbor relationship (count each pair only once)
	neighbors = new vector<int>[tm->numPts()];
	for (i=0; i < tm->numTris(); i++) {
		int v0 = tm->getTri(i, 0);
		int v1 = tm->getTri(i, 1);
		int v2 = tm->getTri(i, 2);
		if (v0 < v1)
			neighbors[v0].push_back(v1);
		if (v1 < v2)
			neighbors[v1].push_back(v2);
		if (v2 < v0)
			neighbors[v2].push_back(v0);
	}

	// build a list of neighbor/influence relations
	for (pt=0; pt < tm->numPts(); pt++) {
		for (i=0; i < neighbors[pt].size(); i++) {
			int pt2 = neighbors[pt][i];
			// only add each edge once
			if (pt2 < pt) 
				continue;
			memset(infMap, 0, sizeof(int)*skin->skel->transforms.size());
			for (j=0; j < skin->maxInf; j++) {
				if (skin->infJoints[pt*skin->maxInf + j] > -1) {
					infMap[skin->infJoints[pt*skin->maxInf + j]] += 1;
					infInd0[skin->infJoints[pt*skin->maxInf + j]] = j;
				}
				if (skin->infJoints[pt2*skin->maxInf + j] > -1) {
					infMap[skin->infJoints[pt2*skin->maxInf + j]] += 2;
					infInd1[skin->infJoints[pt2*skin->maxInf + j]] = j;
				}
			}

			double dist = (tm->getPt(pt) - tm->getPt(pt2)).length();

			for (j=0; j < skin->skel->transforms.size(); j++) {
				if (infMap[j] == 1) {
					nr.set(dist, pt, infInd0[j]);
					neighTable.push_back(nr);
					numZero++;
				}
				else if (infMap[j] == 2) {
					nr.set(dist, pt2, infInd1[j]);
					neighTable.push_back(nr);
					numZero++;
				}
				else if (infMap[j] == 3) {
					nr.set(dist, pt, infInd0[j], pt2, infInd1[j]);
					neighTable.push_back(nr);
					numReg++;
				}
			}
		}
	}

	/*
	double avg = 0;
	for (i = 0; i < neighTable.size(); i++)
		avg += neighTable[i].dist;
	avg /= neighTable.size();
	for (i = 0; i < neighTable.size(); i++) {
		if (neighTable[i].dist <= 1e-6)
			neighTable[i].dist = avg / 1e-6;
		else
			neighTable[i].dist = avg / neighTable[i].dist;
	}

	cout << "found " << numReg << " neighbor pairs, and " << numZero << " zero assignments" << endl;*/

	delete []infMap;
	delete []infInd0;
	delete []infInd1;
}

void USolver::buildPDNeighTable() {
	int pt, pdd, pdd2, neigh, key;
	int curPDD, curOfs, curOfs2;
	bool hasNeigh;
	NeighborRelation nr;

	vector<int> *tmNeigh = findTMNeighbors(uMesh);

	// PDDs
	for (pt = 0; pt < skin->numPts; pt++) {
		for (pdd = 0; pdd < skin->pddPtIndex[pt].size(); pdd += 2) {
			curPDD = skin->pddPtIndex[pt][pdd];
			curOfs = skin->pddPtIndex[pt][pdd + 1];

			for (neigh = 0; neigh < tmNeigh[pt].size(); neigh++) {
				int pt2 = tmNeigh[pt][neigh];
				hasNeigh = false;

				double dist = (uMesh->getPt(pt) - uMesh->getPt(pt2)).length();

				for (pdd2 = 0; pdd2 < skin->pddPtIndex[pt2].size(); pdd2 += 2) {
					if (skin->pddPtIndex[pt2][pdd2] == curPDD) {
						hasNeigh = true;
						if (pt < pt2) {		// (only add each pair once)
							curOfs2 = skin->pddPtIndex[pt2][pdd2 + 1];
							for (key = 0; key < skin->pdds[curPDD]->numSamples-1; key++) {
								nr.set(dist, pt, curOfs + key, pt2, curOfs2 + key);
								pddNeighTable.push_back(nr);
							}
						}
					}
				}

				if (!hasNeigh) {
					// force to zero
					for (key = 0; key < skin->pdds[curPDD]->numSamples-1; key++) {
						nr.set(dist, pt, curOfs + key, -1, -1);
						pddNeighTable.push_back(nr);
					}
				}
			}
		}
	}
	// normalize (and invert) distances
	double avg = 0;
	for (pt = 0; pt < pddNeighTable.size(); pt++)
		avg += pddNeighTable[pt].dist;
	avg /= pddNeighTable.size();
	for (pt = 0; pt < pddNeighTable.size(); pt++) {
		if (pddNeighTable[pt].dist <= 1e-6)
			pddNeighTable[pt].dist = avg / 1e-6;
		else
			pddNeighTable[pt].dist = avg / pddNeighTable[pt].dist;
	}
	cout << "pdd neighbor table size: " << (int)pddNeighTable.size() << endl;

	// PDNs
	for (pt = 0; pt < NORMAL_MAP_SIZE; pt++) {
		for (pdd = 0; pdd < skin->pddNMIndex[pt].size(); pdd += 2) {
			curPDD = skin->pddNMIndex[pt][pdd];
			curOfs = skin->pddNMIndex[pt][pdd + 1];

			for (neigh = 0; neigh < 4; neigh++) {
				int pt2;
				if (neigh == 0) {
					// above
					if (pt / NORMAL_MAP_W > 0)
						pt2 = pt - NORMAL_MAP_W;
					else
						continue;
				}
				else if (neigh == 1) {
					// right
					if (pt % NORMAL_MAP_W < NORMAL_MAP_W - 1)
						pt2 = pt + 1;
					else
						continue;
				}
				else if (neigh == 2) {
					// below
					if (pt / NORMAL_MAP_W < NORMAL_MAP_H - 1)
						pt2 = pt + NORMAL_MAP_W;
					else
						continue;
				}
				else {
					// left
					if (pt % NORMAL_MAP_W > 0)
						pt2 = pt - 1;
					else
						continue;
				}
				hasNeigh = false;

				for (pdd2 = 0; pdd2 < skin->pddNMIndex[pt2].size(); pdd2 += 2) {
					if (skin->pddNMIndex[pt2][pdd2] == curPDD) {
						hasNeigh = true;
						if (neigh < 2) {		// (only add each pair once)
							curOfs2 = skin->pddNMIndex[pt2][pdd2 + 1];
							for (key = 0; key < skin->pdds[curPDD]->numSamples-1; key++) {
								nr.set(1, pt, curOfs + key, pt2, curOfs2 + key);
								pdnNeighTable.push_back(nr);
							}
						}
					}
				}

				if (!hasNeigh) {
					// force to zero
					for (key = 0; key < skin->pdds[curPDD]->numSamples-1; key++) {
						nr.set(1, pt, curOfs + key, -1, -1);
						pdnNeighTable.push_back(nr);
					}
				}
			}
		}
	}
	cout << "pdn neighbor table size: " << (int)pdnNeighTable.size() << endl;

	delete []tmNeigh;
}

void insertCoeffs(vector<int> &vVerts, vector<double> &vCoeffs, int v, double c) {
	int i;
	for (i=0; i < vVerts.size(); i++) {
		if (vVerts[i] == v) {
			vCoeffs[i] += c;
			return;
		}
	}
	vVerts.push_back(v);
	vCoeffs.push_back(c);
}

void USolver::initWeightStencil() {
	baAssert(uMesh->numPts() == skin->numPts, "point size mismatch", true);

	vector<int> *initWPV;// = new vector<int>[skin->numPts * skin->maxInf];
	vector<double> *initWPC;// = new vector<double>[skin->numPts * skin->maxInf];
	double *edgeSums = new double[skin->numPts];
	int pt, edge, ni, i;

	// initialize stencil
	weightPriorVerts = new vector<int>[skin->numPts * skin->maxInf];
	weightPriorCoeffs = new vector<double>[skin->numPts * skin->maxInf];

	initWPV = weightPriorVerts;
	initWPC = weightPriorCoeffs;

	// calculate edge sums
	memset(edgeSums, 0, sizeof(double) * skin->numPts);
	for (pt=0; pt < uMesh->numPts(); pt++) {
		for (edge=0; edge < neighbors[pt].size(); edge++) {
			int pt2 = neighbors[pt][edge];
			// only add each edge once
			if (pt2 < pt)
				continue;

			double edgeLen = (uMesh->getPt(pt) - uMesh->getPt(pt2)).length();
			edgeSums[pt] += edgeLen;
			edgeSums[pt2] += edgeLen;
		}
	}

	// calculate gradient stencil
	for (ni = 0; ni < neighTable.size(); ni++) {
		NeighborRelation &nr = neighTable[ni];

		int v0Ofs = nr.v0 * skin->maxInf + nr.ind0;

		insertCoeffs(initWPV[v0Ofs], initWPC[v0Ofs], v0Ofs, -1.0 / nr.dist);

		if (nr.v1 > -1) {
			int v1Ofs = nr.v1 * skin->maxInf + nr.ind1;

			insertCoeffs(initWPV[v0Ofs], initWPC[v0Ofs], v1Ofs, 1.0 / nr.dist);
			insertCoeffs(initWPV[v1Ofs], initWPC[v1Ofs], v1Ofs, -1.0 / nr.dist);
			insertCoeffs(initWPV[v1Ofs], initWPC[v1Ofs], v0Ofs, 1.0 / nr.dist);
		}
	}
	// divide by edgeSums
	for (ni = 0; ni < skin->numPts * skin->maxInf; ni++) {
		for (i = 0; i < initWPC[ni].size(); i++)
			initWPC[ni][i] *= edgeSums[ni / skin->maxInf];
	}
/*
	// calculate gradient of gradient stencil
	for (ni = 0; ni < neighTable.size(); ni++) {
		NeighborRelation &nr = neighTable[ni];

		int v0Ofs = nr.v0 * skin->maxInf + nr.ind0;

		for (i = 0; i < initWPV[v0Ofs].size(); i++) {
			insertCoeffs(weightPriorVerts[v0Ofs], weightPriorCoeffs[v0Ofs], 
				initWPV[v0Ofs][i], -initWPC[v0Ofs][i] / nr.dist);
		}

		if (nr.v1 > -1) {
			int v1Ofs = nr.v1 * skin->maxInf + nr.ind1;

			for (i = 0; i < initWPV[v0Ofs].size(); i++) {
				insertCoeffs(weightPriorVerts[v1Ofs], weightPriorCoeffs[v1Ofs], 
					initWPV[v0Ofs][i], initWPC[v0Ofs][i] / nr.dist);
			}
			for (i = 0; i < initWPV[v1Ofs].size(); i++) {
				insertCoeffs(weightPriorVerts[v0Ofs], weightPriorCoeffs[v0Ofs], 
					initWPV[v1Ofs][i], initWPC[v1Ofs][i] / nr.dist);
				insertCoeffs(weightPriorVerts[v1Ofs], weightPriorCoeffs[v1Ofs], 
					initWPV[v1Ofs][i], -initWPC[v1Ofs][i] / nr.dist);
			}
		}
	}
	// divide by edgeSums
	for (ni = 0; ni < skin->numPts * skin->maxInf; ni++) {
		for (i = 0; i < weightPriorCoeffs[ni].size(); i++)
			weightPriorCoeffs[ni][i] *= edgeSums[ni / skin->maxInf];
	}

	delete []initWPV;
	delete []initWPC;*/
}

void USolver::dumpDofPos() {
	cout << "Solver variables: " << endl;
	if (dpDress >= 0)
		cout << "dress: " << dpDress << endl;
	if (dpInt >= 0)
		cout << "intrinsics: " << dpInt << endl;
	if (dpPose >= 0)
		cout << "pose: " << dpPose << endl;
	if (dpWeight >= 0)
		cout << "weights: " << dpWeight << endl;
	if (dpPDD >= 0)
		cout << "PDDs: " << dpPDD << endl;
	if (dpX >= 0)
		cout << "PCA: " << dpX << endl;
	if (dpNM >= 0)
		cout << "normal map: " << dpNM << endl;
	if (dpNMP >= 0)
		cout << "pose-dependent normal map: " << dpNMP << endl;
	cout << "total: " << curVars.size() << endl;
}

void USolver::save(FILE *f) {
	int varInfo[2];

	varInfo[0] = 1;
	varInfo[1] = vDress.size();
	fwrite(varInfo, sizeof(int), 2, f);
	fwrite(vDress.n, sizeof(double), vDress.size(), f);

	varInfo[0] = 2;
	varInfo[1] = vInt.size();
	fwrite(varInfo, sizeof(int), 2, f);
	fwrite(vInt.n, sizeof(double), vInt.size(), f);

	varInfo[0] = 3;
	varInfo[1] = vPose.size();
	fwrite(varInfo, sizeof(int), 2, f);
	fwrite(vPose.n, sizeof(double), vPose.size(), f);

	varInfo[0] = 4;
	varInfo[1] = vWeight.size();
	fwrite(varInfo, sizeof(int), 2, f);
	fwrite(vWeight.n, sizeof(double), vWeight.size(), f);

	varInfo[0] = 5;
	varInfo[1] = vPDD.size();
	fwrite(varInfo, sizeof(int), 2, f);
	fwrite(vPDD.n, sizeof(double), vPDD.size(), f);

	varInfo[0] = 6;
	varInfo[1] = dataSet->numCharacters * (numComponents-1);
	fwrite(varInfo, sizeof(int), 2, f);
	int ch;
	for (ch=0; ch < dataSet->numCharacters; ch++) {
		fwrite(&dataSet->charMu[ch][1], sizeof(double), numComponents-1, f);
	}

	if (uMesh) {
		varInfo[0] = 7;
		varInfo[1] = NORMAL_MAP_SIZE * 3 * numComponents / 2;
		fwrite(varInfo, sizeof(int), 2, f);
		fwrite(vNM, sizeof(Vec3d), NORMAL_MAP_SIZE * numComponents / 2, f);
	}

	varInfo[0] = 8;
	varInfo[1] = vNMP.size();
	fwrite(varInfo, sizeof(int), 2, f);
	fwrite(vNMP.n, sizeof(double), vNMP.size(), f);

	varInfo[0] = 0;
	varInfo[1] = 0;
	fwrite(varInfo, sizeof(int), 2, f);
}

void USolver::load(FILE *f) {
	int varInfo[2];
	int numOrigComps;

	while (1) {
		varInfo[0] = 0;
		fread(varInfo, sizeof(int), 2, f);
		if (varInfo[0] == 0)
			break;

		switch (varInfo[0]) {
			case 1:
				cout << "found " << varInfo[1] << " dress vars" << endl;
//				vDress.zeroElements();
				fread(vDress.n, sizeof(double), min(vDress.size(), varInfo[1]), f);
				if (varInfo[1] > vDress.size())
					fseek(f, sizeof(double)*(varInfo[1] - vDress.size()), SEEK_CUR);
				break;
			case 2:
				cout << "found " << varInfo[1] << " intrinsic vars" << endl;
//				vInt.zeroElements();
				fread(vInt.n, sizeof(double), min(vInt.size(), varInfo[1]), f);
				if (varInfo[1] > vInt.size())
					fseek(f, sizeof(double)*(varInfo[1] - vInt.size()), SEEK_CUR);
				break;
			case 3:
				cout << "found " << varInfo[1] << " pose vars" << endl;
//				vPose.zeroElements();
                fread(vPose.n, sizeof(double), min(vPose.size(), varInfo[1]), f);
				if (varInfo[1] > vPose.size())
					fseek(f, sizeof(double)*(varInfo[1] - vPose.size()), SEEK_CUR);
				break;
			case 4:
				cout << "found " << varInfo[1] << " weight vars; expected " << vWeight.size() << endl;
//				vWeight.zeroElements();
				fread(vWeight.n, sizeof(double), min(vWeight.size(), varInfo[1]), f);
				if (varInfo[1] > vWeight.size())
					fseek(f, sizeof(double)*(varInfo[1] - vWeight.size()), SEEK_CUR);
				break;
			case 5:
				cout << "found " << varInfo[1] << " PDD vars" << endl;
//				vPDD.zeroElements();
				fread(vPDD.n, sizeof(double), min(vPDD.size(), varInfo[1]), f);
				if (varInfo[1] > vPDD.size())
					fseek(f, sizeof(double)*(varInfo[1] - vPDD.size()), SEEK_CUR);
				break;
			case 6:
				if (varInfo[1] % dataSet->numCharacters != 0) {
					cout << "found " << varInfo[1] << " X vars -- incompatible" << endl;
				}
				else {
					cout << "found " << varInfo[1] << " X vars" << endl;
					numOrigComps = varInfo[1] / dataSet->numCharacters;
					cout << numOrigComps << " " << numComponents << endl;
					int ch;
					for (ch=0; ch < dataSet->numCharacters; ch++) {
						fread(&dataSet->charMu[ch][1], sizeof(double), 
							min(numOrigComps, numComponents-1), f);
						if (numOrigComps > numComponents-1) {
							fseek(f, sizeof(double) * (numOrigComps - (numComponents-1)), SEEK_CUR);
						}
						else if (numOrigComps < numComponents-1) {
							int c;
							for (c = numOrigComps+1; c < numComponents; c++)
								dataSet->charMu[ch][c] = boundedRand(-0.5, 0.5);
						}
					}
				}
				break;
			case 7:
				if (uMesh) {
					cout << "found " << varInfo[1] << " normal map vars" << endl;
					fread(vNM, sizeof(double), min(NORMAL_MAP_SIZE * 3 * numComponents / 2, varInfo[1]), f);
					if (varInfo[1] > vPDD.size())
						fseek(f, sizeof(double)*(varInfo[1] - (NORMAL_MAP_SIZE * 3 * numComponents / 2)), SEEK_CUR);
				}
				else
					fseek(f, varInfo[1] * sizeof(double), SEEK_CUR);
				break;
			case 8:
				cout << "found " << varInfo[1] << " NMP vars" << endl;
//				vNMP.zeroElements();
				fread(vNMP.n, sizeof(double), min(vNMP.size(), varInfo[1]), f);
				if (varInfo[1] > vNMP.size())
					fseek(f, sizeof(double)*(varInfo[1] - vNMP.size()), SEEK_CUR);
				break;
			default:
				cout << "unknown variable type: " << varInfo[0] << endl;
				fseek(f, varInfo[1] * sizeof(double), SEEK_CUR);
		}
	}
}

bool USolver::saveDress(char *fname) {
	FILE *f;
	if (!openFile(&f, fname, "wb", "dress file"))
		return false;
	fwrite(vDress.n, sizeof(double), vDress.size(), f);
	fclose(f);
	return true;
}

bool USolver::loadDress(char *fname, char *fname2) {
	FILE *f;
	if (!openFile(&f, fname, "rb", "dress file"))
		return false;

	if (!fname2 || fname2[0] == 0) {
		fread(vDress.n, sizeof(double), vDress.size(), f);
	}
	else {
		int numSmallPts, numBigPts;
		int *mapPts;
		float *mapWeights;

		if (!loadMapping(fname2, numSmallPts, numBigPts, numPts, mapPts, mapWeights))
			return false;

		double *origVals = new double[numSmallPts * numComponents * 3];
		fread(origVals, sizeof(double), numSmallPts * numComponents * 3, f);
		vDress.zeroElements();

		int comp, pt, ind, inf;
		for (comp = 0; comp < numComponents; comp++) {
			for (pt = 0; pt < numOrigPts; pt++) {
				for (inf = 0; inf < MAX_MAPPING_SIZE; inf++) {
					int mapPt = mapPts[pt * MAX_MAPPING_SIZE + inf];
					if (mapPt >= 0) {
						for (ind = 0; ind < 3; ind++) {
							vDress[comp * numOrigPts * 3 + ind] += 
								mapWeights[MAX_MAPPING_SIZE + inf] *
								origVals[mapPt * 3 + ind];
						}
					}
				}
			}
		}

		delete []origVals;
	}

	fclose(f);
	return true;
}

bool USolver::savePoses(char *fname, int mode) {
	int ch = 1;
	FILE *f;
	if (!openFile(&f, fname, "wb", "poses file"))
		return false;

	if (mode == 0 || mode == 1) {
		fwrite(&numComponents, sizeof(int), 1, f);
		fwrite(vInt.n, sizeof(double), vInt.size(), f);
	}
	else {
		// save this character's intrinsics into the mean
		int comp, dof;
		int ip = 0;
		comp = 1;
		fwrite(&comp, sizeof(int), 1, f);	// only one component saved
		double *intrins = new double[skin->skel->numIntrinsicDofs];
		memset(intrins, 0, sizeof(double) * skin->skel->numIntrinsicDofs);
		for (comp=0; comp < numComponents; comp++) {
			double w = dataSet->charMu[ch][comp];
			for (dof = 0; dof < skin->skel->numIntrinsicDofs; dof++) {
				intrins[dof] += w * vInt[ip++];
			}
		}
		fwrite(intrins, sizeof(double), skin->skel->numIntrinsicDofs, f);
		delete []intrins;
	}

	if (mode == 0) {
		fwrite(&dataSet->numExamples, sizeof(int), 1, f);
		fwrite(vPose.n, sizeof(double), vPose.size(), f);
	}
	else if (mode == 2) {
		// just the poses for this character
		int ex, dof;
		int first = dataSet->charIndex[ch];
		int last = first;
		while (last < dataSet->numExamples && dataSet->examples[last].character == ch)
			last++;
		ex = last - first;
		fwrite(&ex, sizeof(int), 1, f);

		int pp = first*skin->skel->numPoseDofs;
		for (ex=first; ex < last; ex++) {
			fwrite(vPose.n + ex*skin->skel->numPoseDofs, sizeof(double), skin->skel->numPoseDofs, f);
		}
	}

	fclose(f);
	return true;
}

bool USolver::loadPoses(char *fname, bool intOnly) {
	int n;
	FILE *f;
	if (!openFile(&f, fname, "rb", "poses file"))
		return false;

	fread(&n, sizeof(int), 1, f);
	cout << "loading " << min(n, numComponents) << " components" << endl;
	fread(vInt.n, sizeof(double), skin->skel->numIntrinsicDofs * min(n, numComponents), f);
	if (numComponents < n)
		fseek(f, sizeof(double) * skin->skel->numIntrinsicDofs * (n - numComponents), SEEK_CUR);

	if (!intOnly) {
		fread(&n, sizeof(int), 1, f);
		int ex, dof;
		cout << "loading " << min(n, dataSet->numExamples) << " examples" << endl;
		for (ex=0; ex < min(n, dataSet->numExamples); ex++) {
			fread(vPose.n + ex*skin->skel->numPoseDofs, sizeof(double), skin->skel->numPoseDofs, f);
		}
	}
	fclose(f);
	return true;
}

void USolver::weightsToVars() {
	int i, j;
	int wp = 0;
	for (i=0; i < numOrigPts; i++) {
		for (j=1; j < skin->maxInf; j++)
			vWeight[wp++] = skin->infWeights[i*skin->maxInf + j] / WEIGHT_DOF_SCALE;
	}
}

void USolver::dressToVars(int comp) {
	int pt;
	int dp = 0;

	if (comp >= 0)
		dp = numOrigPts * 3 * comp;
	else
		vDress.zeroElements();

	for (pt = 0; pt < numOrigPts; pt++) {
		double pWeight = 0.5;
		double mWeight;
		if (skin->dressPts[pt].iszero())
			pWeight = 0;
		mWeight = 1 - pWeight;
		if (mirrorMap[pt] == -1 || skin->dressPts[mirrorMap[pt]].iszero()) {
			mWeight = 0;
			if (pWeight > 0)
				pWeight = 1;
		}
		
		// X
		vDress[dp++] = 
			pWeight * skin->dressPts[pt][0] + 
			mWeight * skin->dressPts[mirrorMap[pt]][0];

		// Y (midline points must have y = 0)
		if (pt == mirrorMap[pt])
			vDress[dp++] = 0;
		else
			vDress[dp++] = 
				pWeight * skin->dressPts[pt][1] - 
				mWeight * skin->dressPts[mirrorMap[pt]][1];

		// Z
		vDress[dp++] = 
			pWeight * skin->dressPts[pt][2] + 
			mWeight * skin->dressPts[mirrorMap[pt]][2];
	}
}

void USolver::updateSkel(double *comps, int numComps, int poseInd) {
	// if poseInd is -1, show dress pose
	// if poseInd is >= 0, show that example's pose
	// otherwise, leave the pose alone
	int i, j;

	// reconstruct intrinsic DOFs from components
	static double *recon = NULL;
	if (!recon)
		recon = new double[skin->skel->numIntrinsicDofs];
	memset(recon, 0, sizeof(double) * skin->skel->numIntrinsicDofs);
	for (i=0; i < min(numComps, numComponents); i++) {
		for (j=0; j < skin->skel->numIntrinsicDofs; j++) {
			recon[j] += comps[i] * vInt[i * skin->skel->numIntrinsicDofs + j];
		}
	}

	// load DOFs into skeleton
	int ip = 0;
	int pp = -1;
	if (poseInd >= 0)
		pp = skin->skel->numPoseDofs * poseInd;
	for (i=0; i < skin->skel->transforms.size(); i++) {
		SkelTransform *curTrans = skin->skel->transforms.getT(i);
		if (curTrans->isIntrinsic) {
			curTrans->loadDofs(recon + ip);
			ip += curTrans->numDofs();
		}
		else if (poseInd == -1) {
			curTrans->zero();
		}
		else if (pp > -1) {
			curTrans->loadDofs(vPose.n + pp);
			pp += curTrans->numDofs();
		}
	}
	skin->skel->updateCoords();
}

void USolver::updateSkel(int ex) {
	updateSkel(&(dataSet->charMu[dataSet->examples[ex].character][0]), 
		numComponents, ex);
}

void USolver::updateWeights() {
	int i, j, comp;
	double sum;

	int wp = 0;
	for (i=0; i < numOrigPts; i++) {
		sum = 0;
		for (j=1; j < skin->maxInf; j++) {
			if (skin->infJoints[i*skin->maxInf + j] >= 0) {
				sum += vWeight[wp];
				skin->infWeights[i*skin->maxInf + j] = vWeight[wp];
				if (numOrigPts < numPts && mirrorMap[i] != i) {
					skin->infWeights[mirrorMap[i]*skin->maxInf + j] = vWeight[wp];
				}
			}
			wp++;
		}
		skin->infWeights[i*skin->maxInf] = 1.0 - sum;
		if (numOrigPts < numPts && mirrorMap[i] != i)
			skin->infWeights[mirrorMap[i]*skin->maxInf] = 1.0 - sum;
	}
}

void USolver::updatePoints(int ch) {
	updatePoints(dataSet->charMu[ch].Ref(), numComponents);
}

void USolver::updatePoints(double *comps, int numComps, bool pddOnly) {
	int comp, pt, var, numVars;

	numComps = min(numComps, numComponents);

	if (uMesh) {
		// reconstruct normal map
		for (pt=0; pt < NORMAL_MAP_SIZE / 2; pt++) {
			skin->baseNormalMap[pt] = vNM[pt];
			for (comp = 1; comp < numComps; comp++ ){
				double w = comps[comp];
				skin->baseNormalMap[pt] += w * vNM[comp * NORMAL_MAP_SIZE / 2 + pt];
			}

			// mirrored part
			int x = pt % NORMAL_MAP_W;
			int y = pt / NORMAL_MAP_W;
			skin->baseNormalMap[(NORMAL_MAP_H - 1 - y) * NORMAL_MAP_W + x] = skin->baseNormalMap[pt];
			skin->baseNormalMap[(NORMAL_MAP_H - 1 - y) * NORMAL_MAP_W + x][1] *= -1;
		}

		// reconstruct normal map PDDs
		memcpy(skin->pddNMKeys, vNMP.n, sizeof(Vec3d) * skin->numPddNMKeys);
		int index = skin->numPddNMKeys * 3;
		for (comp = 1; comp < min(numComps, numPDComponents); comp++ ){
			double w = comps[comp];
			if (w != 0) {
				for (pt=0; pt < skin->numPddNMKeys; pt++) {
					skin->pddNMKeys[pt] += w * Vec3d(vNMP[index+0], vNMP[index+1], vNMP[index+2]);
					index += 3;
				}
			}
		}
	}

	if (!pddOnly) {
		// reconstruct dress points
		memset(skin->dressPts, 0, sizeof(Vec3d) * numPts);
		for (comp = 0; comp < numComps; comp++) {
			double w = comps[comp];
			if (w != 0) {
				for (pt=0; pt < numPts; pt++) {
					skin->dressPts[pt] += w * varDressPt(pt, comp);
				}
			}
		}
	}

	// reconstruct PDDs
	if (vPDD.size() > 0) {
		memcpy(skin->pddPtKeys, vPDD.n, sizeof(Vec3d) * skin->numPddPtKeys);
		int ofs = skin->numPddPtKeys * 3;
		for (comp = 1; comp < min(numComps, numPDComponents); comp++) {
			for (var = 0; var < skin->numPddPtKeys; var++) {
				skin->pddPtKeys[var] += comps[comp] *
					Vec3d(vPDD[ofs+0], vPDD[ofs+1], vPDD[ofs+2]);
				ofs += 3;
			}
		}
	}
}

void USolver::rayPoint(Vec3d pt, Vec3d norm, UExample *ex, double &curDistance, Vec3d &curClosestPt, double &curSurfaceWeight) {
	curClosestPt = Vec3d();
	curDistance = 0;
	curSurfaceWeight = 0;

	ex->mesh->closestRestrictNormal = true;
	ex->mesh->closestNormalRestriction = norm;

	if (ex->mesh->calcClosestPoint(pt, 1.10)) {
		curClosestPt = ex->mesh->closestPt;

		Vec3d delta = ex->mesh->closestPt - pt;
		delta.normalize();

		// if closest point is a vertex, check if we need to flip the sign
		if (ex->mesh->closestTri[1] == -1) {
			curSurfaceWeight = ex->mesh->getPtConf(ex->mesh->closestTri[0]);

			double dotProd = delta * ex->mesh->getPtNormal(ex->mesh->closestTri[0]);
			if (dotProd < 0)
				ex->mesh->closestDist *= -1;

			// check if this vertex adjoins a hole
			if (ex->vertList[ex->mesh->closestTri[0]]) {
				curSurfaceWeight = 0;
				return;
			}

			// check normal match
			if (ex->mesh->closestNormalRestriction * ex->mesh->getPtNormal(ex->mesh->closestTri[0]) < NORMAL_TOL) {
				curSurfaceWeight = 0;
				return;
			}
		}
		// if closest point is an edge, check if we need to flip the sign
		else if (ex->mesh->closestTri[2] == - 1) {
			curSurfaceWeight = ex->mesh->closestBary[0] * ex->mesh->getPtConf(ex->mesh->closestTri[0]) + ex->mesh->closestBary[1] * ex->mesh->getPtConf(ex->mesh->closestTri[1]);

			EdgeInfo *ei = ex->edgeList->findEdge(ex->mesh->closestTri[0], ex->mesh->closestTri[1]);
			if (!ei) {
				cout << "error: unknown edge!!" << endl;
				return;
			}
			double dotProd = delta * ei->normal;
			if (dotProd < 0)
				ex->mesh->closestDist *= -1;

			// check if this edge adjoins a hole
			if (ei->count == 1) {
				curSurfaceWeight = 0;
				return;
			}

			// check normal match
			if (ex->mesh->closestNormalRestriction * ei->normal < NORMAL_TOL) {
				curSurfaceWeight = 0;
				return;
			}
		}
		else {
			curSurfaceWeight = 
				ex->mesh->closestBary[0] * ex->mesh->getPtConf(ex->mesh->closestTri[0]) + 
				ex->mesh->closestBary[1] * ex->mesh->getPtConf(ex->mesh->closestTri[1]) + 
				ex->mesh->closestBary[2] * ex->mesh->getPtConf(ex->mesh->closestTri[2]);

			Vec3d verts[3];
			verts[0] = ex->mesh->getPt(ex->mesh->closestTri[0]);
			verts[1] = ex->mesh->getPt(ex->mesh->closestTri[1]);
			verts[2] = ex->mesh->getPt(ex->mesh->closestTri[2]);
			Vec3d norm = -(verts[1] - verts[0]) ^ (verts[2] - verts[0]);
			norm.normalize();

			// check normal match
			if (ex->mesh->closestNormalRestriction * norm < NORMAL_TOL) {
				curSurfaceWeight = 0;
				return;
			}
		}

		curDistance = ex->mesh->closestDist;
	}
}

void USolver::adjustWeights() {
	double weights[10];	// must be bigger than max inf
	int i, j, comp;
	double sum, adjSum;

	// scale each weight by the inverse of the distance between the 
	// point and the joint, then renormalize

	int wp = 0;
	for (i=0; i < numOrigPts; i++) {
		sum = 0;
		adjSum = 0;

		weights[0] = 1;
		for (j=1; j < skin->maxInf; j++) {
			if (skin->infJoints[i*skin->maxInf + j] >= 0) {
				sum += vWeight[wp];
//				double w = vWeight[wp] / (uMesh->getPt(i) - skin->curJointPos[skin->infJoints[i*skin->maxInf + j]]).length();
				double w = sqr(vWeight[wp]);
				weights[j] = w;
				adjSum += w;
			}
			wp++;
		}
		weights[0] = sqr(1.0 - sum);
//		weights[0] = (1.0 - sum) / (uMesh->getPt(i) - skin->curJointPos[skin->infJoints[i*skin->maxInf + 0]]).length();
		adjSum += weights[0];

		for (j=1; j < skin->maxInf; j++) {
			if (skin->infJoints[i*skin->maxInf + j] >= 0) {
				vWeight[wp - skin->maxInf + j] = weights[j] / adjSum;
			}
		}
	}
}

void USolver::initVars(bool verbose) {
	int numVars = 0;

	if (optDress) {
		dpDress = numVars;
		numVars += vDress.size();
	}
	else
		dpDress = -1;
	
	if (optInt) {
		dpInt = numVars;
		numVars += vInt.size();
	}
	else
		dpInt = -1;

	if (optPose) {
		dpPose = numVars;
		numVars += vPose.size();
	}
	else
		dpPose = -1;

	if (optWeight) {
		dpWeight = numVars;
		numVars += vWeight.size();
	}
	else
		dpWeight = -1;

	if (optPDD) {
		dpPDD = numVars;
		numVars += vPDD.size();
	}
	else
		dpPDD = -1;

	if (optX) {
		dpX = numVars;
		numVars += dataSet->numCharacters * (numComponents-1);
	}
	else
		dpX = -1;

	if (optNM) {
		dpNM = numVars;
		numVars += 3 * NORMAL_MAP_SIZE / 2 * numComponents;
	}
	else
		dpNM = -1;

	if (optNMP) {
		dpNMP = numVars;
		numVars += vNMP.size();
	}
	else
		dpNMP = -1;

	curVars.resize(numVars, true);
	grad.resize(numVars);
	if (verbose)
		dumpDofPos();

	toVars();
}

void USolver::toVars() {
	int i, comp;
	if (optDress)
		memcpy(curVars.n + dpDress, vDress.n, sizeof(double) * vDress.size());
	if (optInt)
		memcpy(curVars.n + dpInt, vInt.n, sizeof(double) * vInt.size());
	if (optPose) {
		for (i=0; i < vPose.size(); i++)
			curVars[dpPose + i] = vPose[i] / POSE_DOF_SCALE;
	}
	if (optWeight) {
		for (i=0; i < vWeight.size(); i++)
			curVars[dpWeight + i] = vWeight[i] / WEIGHT_DOF_SCALE;
	}
	if (optPDD)
		memcpy(curVars.n + dpPDD, vPDD.n, sizeof(double) * vPDD.size());
	if (optX) {
		int ind = 0;
		for (i=0; i < dataSet->numCharacters; i++) {
			for (comp=1; comp < numComponents; comp++) {
				curVars[dpX + ind] = dataSet->charMu[i][comp];
				ind++;
			}
		}
	}
	if (optNM && uMesh) {
		memcpy(curVars.n + dpNM, vNM, sizeof(Vec3d) * NORMAL_MAP_SIZE / 2 * numComponents);
	}
	if (optNMP) {
		memcpy(curVars.n + dpNMP, vNMP.n, sizeof(double) * vNMP.size());
	}
}

void USolver::fromVars(Vecd &vars) {
	int i, comp;
	if (optDress)
		memcpy(vDress.n, vars.n + dpDress, sizeof(double) * vDress.size());
	if (optInt)
		memcpy(vInt.n, vars.n + dpInt, sizeof(double) * vInt.size());
	if (optPose) {
		for (i=0; i < vPose.size(); i++)
			vPose[i] = vars[dpPose + i] * POSE_DOF_SCALE;
	}
	if (optWeight) {
		for (i=0; i < vWeight.size(); i++)
			vWeight[i] = vars[dpWeight + i] * WEIGHT_DOF_SCALE;
	}
	if (optPDD)
		memcpy(vPDD.n, vars.n + dpPDD, sizeof(double) * vPDD.size());
	if (optX) {
		int ind = 0;
		for (i=0; i < dataSet->numCharacters; i++) {
			for (comp=1; comp < numComponents; comp++) {
				dataSet->charMu[i][comp] = vars[dpX + ind];
				ind++;
			}
		}
	}
	if (optNM && uMesh) {
		memcpy(vNM, curVars.n + dpNM, sizeof(Vec3d) * NORMAL_MAP_SIZE / 2 * numComponents);
	}
	if (optNMP) {
		memcpy(vNMP.n, curVars.n + dpNMP, sizeof(double) * vNMP.size());
	}
}

void USolver::applyDfdv(int ex) {
	static Vec3d *weightGrad = NULL;
	if (!weightGrad)
		weightGrad = new Vec3d[skin->maxInf];
	Vec3d dressGrad[3];
	static Vec3d *pddGrad = new Vec3d[skin->numPddPtKeys * 3];
	static Vec3d *intGrad = new Vec3d[skin->skel->numIntrinsicDofs];
	static Vec3d *poseGrad = new Vec3d[skin->skel->numPoseDofs];

	int pt, comp, var, ofs, inf;

	for (pt = 0; pt < skin->numPts; pt++) {
		if (dfdv[pt].iszero())
			continue;

		if (optDress || optX)
			memset(dressGrad, 0, sizeof(Vec3d)*3);
		if (optInt || optX)
			memset(intGrad, 0, sizeof(Vec3d)*skin->skel->numIntrinsicDofs);
		if (optPose)
			memset(poseGrad, 0, sizeof(Vec3d)*skin->skel->numPoseDofs);
		if (optWeight)
			memset(weightGrad, 0, sizeof(Vec3d)*skin->maxInf);
// since pddGrad is so big, we'll clear the entries as they're used
//		if (optPDD || optX)
//			memset(pddGrad, 0, sizeof(Vec3d)*skin->numPddPtKeys*3);

		skin->calcVecGrad(pt, 
			optDress || optX ? dressGrad   : NULL, 
			optInt || optX   ? intGrad     : NULL, 
			optPose          ? poseGrad : NULL, 
			optWeight        ? weightGrad  : NULL, 
			optPDD || optX   ? pddGrad     : NULL);

		if (optPose) {
			ofs = dpPose + skin->skel->numPoseDofs * ex;
			for (var = 0; var < skin->skel->numPoseDofs; var++) {
				grad[ofs + var] += dfdv[pt] * poseGrad[var];
			}
		}

		for (comp = minComp; comp < maxComp; comp++) {
			double w = dataSet->charMu[dataSet->examples[ex].character][comp];

			if (w != 0) {
				if (optDress) {
					// fill in dress shape gradients
					if (pt >= numOrigPts) {
						ofs = dpDress + (comp * numOrigPts + mirrorMap[pt]) * 3;
						grad[ofs + 0] += w * dfdv[pt] * dressGrad[0];
						grad[ofs + 1] -= w * dfdv[pt] * dressGrad[1];
						grad[ofs + 2] += w * dfdv[pt] * dressGrad[2];
					}
					else {
						ofs = dpDress + (comp * numOrigPts + pt) * 3;
						grad[ofs + 0] += w * dfdv[pt] * dressGrad[0];
						if (mirrorMap[pt] != pt)
							grad[ofs + 1] += w * dfdv[pt] * dressGrad[1];
						grad[ofs + 2] += w * dfdv[pt] * dressGrad[2];
					}
				}

				if (optInt) {
					// fill in intrinsic DOF gradients
					ofs = dpInt + comp * skin->skel->numIntrinsicDofs;
					for (var = 0; var < skin->skel->numIntrinsicDofs; var++)
						grad[ofs + var] += w * dfdv[pt] * intGrad[var];
				}
			}

			if (optX && comp > 0) {
				// fill in PCA reconstruction weight gradients
				ofs = dpX + dataSet->examples[ex].character * (numComponents - 1);
				grad[ofs + comp - 1] += varDressPt(pt, comp) * 
					Vec3d(dfdv[pt] * dressGrad[0], dfdv[pt] * dressGrad[1], dfdv[pt] * dressGrad[2]);
				for (var = 0; var < skin->skel->numIntrinsicDofs; var++)
					grad[ofs + comp - 1] += vInt[comp * skin->skel->numIntrinsicDofs + var] * 
						dfdv[pt] * intGrad[var];

				if (comp < this->numPDComponents) {
					for (var = 0; var < skin->pddPtIndex[pt].size(); var += 2) {
						UPoseDepDef *pdd = skin->pdds[skin->pddPtIndex[pt][var]];
						int pddOfs = (skin->pddPtIndex[pt][var+1]) * 3;
						for (inf = 0; inf < (pdd->rbf->mN - 1) * 3; inf++) {
							grad[ofs + comp - 1] += 
								vPDD[(comp * skin->numPddPtKeys * 3) + pddOfs] * 
								dfdv[pt] * pddGrad[pddOfs];
							pddOfs++;
						}
					}
				}
			}
		}

		if (optWeight) {
			// fill in weight gradients
			if (pt >= numOrigPts)
				ofs = dpWeight + mirrorMap[pt]*(skin->maxInf - 1);
			else
				ofs = dpWeight + pt*(skin->maxInf - 1);
			for (inf = 1; inf < skin->maxInf; inf++) {
				if (skin->infJoints[pt*skin->maxInf + inf] >= 0) {
					grad[ofs + inf - 1] += dfdv[pt] * weightGrad[inf];
					grad[ofs + inf - 1] -= dfdv[pt] * weightGrad[0];
				}
			}
		}

		if (optPDD) {
			// fill in PDD gradients
			for (comp = minComp; comp < maxPComp; comp++) {
				double w = dataSet->charMu[dataSet->examples[ex].character][comp];
				if (w != 0) {
					for (var = 0; var < skin->pddPtIndex[pt].size(); var += 2) {
						UPoseDepDef *pdd = skin->pdds[skin->pddPtIndex[pt][var]];
						int ofs = (skin->pddPtIndex[pt][var+1]) * 3;
						for (inf = 0; inf < (pdd->rbf->mN - 1) * 3; inf++) {
							grad[dpPDD + (comp * skin->numPddPtKeys * 3) + ofs] += 
								w * dfdv[pt] * pddGrad[ofs];
							ofs++;
						}
					}
				}
			}
		}

		if (optPDD || optX) {
			// blank out PDD values
			for (var = 0; var < skin->pddPtIndex[pt].size(); var += 2) {
				UPoseDepDef *pdd = skin->pdds[skin->pddPtIndex[pt][var]];
				int ofs = (skin->pddPtIndex[pt][var+1]) * 3;
				for (inf = 0; inf < (pdd->rbf->mN - 1) * 3; inf++) {
					pddGrad[ofs] = Vec3d();
					ofs++;
				}
			}
		}
	}
}

void USolver::calcNMErr(Vec3d &curVal, Vec3d &targetVal, int ex, int pt, int gradOfs, double yMult) {
	int comp, inf, w;
	
	// Fisher distribution
	double cLen = curVal.length();
	double dLen = targetVal.length();
	double delta = (curVal * targetVal) / (cLen * dLen);
	nmErr += 1.0 - delta;
	numNM += 1;
	Vec3d deriv(
		(cLen * dLen) * targetVal[0] - delta * dLen * dLen * curVal[0],
		(cLen * dLen) * targetVal[1] - delta * dLen * dLen * curVal[1],
		(cLen * dLen) * targetVal[2] - delta * dLen * dLen * curVal[2]);
	deriv = deriv / (cLen * cLen * dLen * dLen);
	if (optNM) {
		for (comp = 0; comp < maxPComp; comp++) {
			double w = dataSet->charMu[dataSet->examples[ex].character][comp];
			grad[comp * NORMAL_MAP_SIZE * 3 / 2 + gradOfs + 0] -= w * nmSigma * deriv[0];
			grad[comp * NORMAL_MAP_SIZE * 3 / 2 + gradOfs + 1] -= w * yMult * nmSigma * deriv[1];
			grad[comp * NORMAL_MAP_SIZE * 3 / 2 + gradOfs + 2] -= w * nmSigma * deriv[2];
		}
	}
	if (optNMP) {
		for (inf=0; inf < skin->pddNMIndex[pt].size(); inf += 2) {
			RBF *rbf = skin->pdds[skin->pddNMIndex[pt][inf]]->rbf;
			int ofs = skin->pddNMIndex[pt][inf+1];
			double yMult = 1.0;
			if (skin->pdds[skin->pddNMIndex[pt][inf]]->isMirror)
				yMult = -1.0;
			for (w=1; w < rbf->mN; w++) {
				for (comp = 0; comp < maxPComp; comp++) {
					double cw = dataSet->charMu[dataSet->examples[ex].character][comp];
					grad[dpNMP + (comp * skin->numPddNMKeys + ofs)*3 + 0] -= 
						cw * nmSigma * rbf->curWeights[w] * deriv[0];
					grad[dpNMP + (comp * skin->numPddNMKeys + ofs)*3 + 1] -= 
						cw * yMult * nmSigma * rbf->curWeights[w] * deriv[1];
					grad[dpNMP + (comp * skin->numPddNMKeys + ofs)*3 + 2] -= 
						cw * nmSigma * rbf->curWeights[w] * deriv[2];
				}
				ofs++;
			}
		}
	}

	// Gaussian distribution
/*	Vec3d delta = curVal - targetVal;
	nmErr += delta.length2();
	numNM += 3;
	if (optNM) {
		for (comp = 0; comp < numComponents; comp++) {
			double w = dataSet->charMu[dataSet->examples[ex].character][comp];
			grad[comp * NORMAL_MAP_SIZE * 3 / 2 + gradOfs + 0] += w * nmSigma * 2.0 * delta[0];
			grad[comp * NORMAL_MAP_SIZE * 3 / 2 + gradOfs + 1] += w * yMult * nmSigma * 2.0 * delta[1];
			grad[comp * NORMAL_MAP_SIZE * 3 / 2 + gradOfs + 2] += w * nmSigma * 2.0 * delta[2];
		}
	}
	if (optNMP) {
		for (inf=0; inf < skin->pddNMIndex[pt].size(); inf += 2) {
			RBF *rbf = skin->pdds[skin->pddNMIndex[pt][inf]]->rbf;
			int ofs = skin->pddNMIndex[pt][inf+1];
			double yMult = 1.0;
			if (skin->pdds[skin->pddNMIndex[pt][inf]]->isMirror)
				yMult = -1.0;
			for (w=1; w < rbf->mN; w++) {
				for (comp = 0; comp < numPDComponents; comp++) {
					double cw = dataSet->charMu[dataSet->examples[ex].character][comp];
					grad[dpNMP + (comp * skin->numPddNMKeys + ofs)*3 + 0] += cw * nmSigma * 2.0 * delta[0] * rbf->curWeights[w];
					grad[dpNMP + (comp * skin->numPddNMKeys + ofs)*3 + 1] += cw * yMult * nmSigma * 2.0 * delta[1] * rbf->curWeights[w];
					grad[dpNMP + (comp * skin->numPddNMKeys + ofs)*3 + 2] += cw * nmSigma * 2.0 * delta[2] * rbf->curWeights[w];
				}
				ofs++;
			}
		}
	}*/
}

double USolver::evaluateFunction(Vecd& variables) {
	int ex, pt, tr, ni, inf, samp, ch, dof, ofs, var, tri, w;
	int comp, comp2;
	double mult, ret = 0;
	Vec3d delta;

	numVert = 0; vertErr = 0;
	numNorm = 0; normErr = 0;
	numWeightNeigh = 0; weightNeighErr = 0;
	numNMNeigh = 0; nmNeighErr = 0;
	numPddNeigh = 0; pddNeighErr = 0;
	numPdnNeigh = 0; pdnNeighErr = 0;
	numDressMuNeigh = 0; dressMuNeighErr = 0;
	numDressWNeigh = 0; dressWNeighErr = 0;
	numMarkerSpring = 0; markerSpringErr = 0;
	numNM = 0; nmErr = 0;

	if (singleComp >= 0) {
		minComp = singleComp;
		maxComp = singleComp + 1;
		maxPComp = min(numPDComponents, singleComp+1);
	}
	else {
		minComp = 0;
		maxComp = numComponents;
		maxPComp = numPDComponents;
	}

	// copy current values
	fromVars(variables);
	// zero gradient
	grad.zeroElements();

	updateWeights();

	// point-distance terms
	for (ex=0; ex < dataSet->numExamples; ex++) {
		memset(dfdv, 0, sizeof(Vec3d) * numPts);

		updateSkel(ex);
		updatePoints(dataSet->examples[ex].character);
		skin->skel->updateDerivs();
		skin->skel->updateGlobalDerivs();
		skin->updateMats();
		skin->updateJoints();	// could check if this is 1st or not
		skin->updatePts(0, -1, optNM || optNMP);

		if (dataSet->examples[ex].lookup) {
			for (pt=0; pt < uMesh->numPts(); pt++)
				uMesh->getPt(pt) = skin->curPts[pt];
			uMesh->calcNormals();
			for (pt=0; pt < dataSet->examples[ex].luSize; pt++) {
				int ind = dataSet->examples[ex].lookup[pt];
				if (ind >= 0) {
					dataSet->examples[ex].luOffsets[ind] = uMesh->getPtNormal(pt) * 0.01;
				}
			}
		}

		if (optDress || optInt || optPose || optWeight || optPDD || optX) {
			if (!ignorePoints) {
				// point error
				for (pt = 0; pt < numPts; pt++) {
					if (markerAssignments) {
						int ma;
						for (ma = 0; ma < numMarkers; ma++) {
							if (dataSet->examples[ex].ptsConf[ma] == 0)
								continue;

							Vec3d exPt = dataSet->examples[ex].pts[ma];
							double mult = markerAssignments[ma*numPts + pt];

							delta = skin->curPts[pt] - exPt;
							vertErr += mult * delta.length2();
							numVert += 3 * mult;

							dfdv[pt] += 2.0 * vertSigma * mult * delta;
						}
					}
					else { //if (ex == 16) {
						// don't go past the number of points in this set
						if (!dataSet->examples[ex].lookup && pt > dataSet->examples[ex].numPts)
							break;

						Vec3d exPt;
						dataSet->examples[ex].getPt(pt, &exPt, &mult);

						delta = skin->curPts[pt] - exPt;
						vertErr += mult * delta.length2();
						numVert += 3 * mult;

//						if (ex == 16) {
//							uMesh->getPtColor(pt) = Vec3d(fabs(delta[0]), fabs(delta[1]), fabs(delta[2])) * 10;
//						}

						dfdv[pt] += 2.0 * vertSigma * mult * delta;
					}
				}

				if (uMesh) {
#ifdef DO_NORMAL_ERR
					// normal error
					for (tri = 0; tri < uMesh->numTris(); tri++) {
						int verts[3];
						double d;
						Vec3d ep[3], mp[3], eNorm, mNorm;
						bool skipTri = false;

						mult = 0;
						for (comp = 0; comp < 3; comp++) {
							verts[comp] = uMesh->getTri(tri, comp);
							dataSet->examples[ex].getPt(verts[comp], &ep[comp], &d); 
							if (ep[comp].iszero())
								skipTri = true;
							mult += d;
							mp[comp] = skin->curPts[verts[comp]];
						}
						if (skipTri)
							continue;
						mult /= 3.0;

						eNorm = (ep[1] - ep[0]) ^ (ep[2] - ep[0]);
						mNorm = (mp[1] - mp[0]) ^ (mp[2] - mp[0]);
						delta = mNorm - eNorm;
						normErr += mult * delta.length2();
						numNorm += 3 * mult;
			/*
						dfdv[verts[0]][0] += 2.0 * mult * normSigma * delta * ((mp[1] - mp[0]) ^ Vec3d(-1, 0, 0) + Vec3d(-1, 0, 0) ^ (mp[2] - mp[0]));
						dfdv[verts[0]][1] += 2.0 * mult * normSigma * delta * ((mp[1] - mp[0]) ^ Vec3d(0, -1, 0) + Vec3d(0, -1, 0) ^ (mp[2] - mp[0]));
						dfdv[verts[0]][2] += 2.0 * mult * normSigma * delta * ((mp[1] - mp[0]) ^ Vec3d(0, 0, -1) + Vec3d(0, 0, -1) ^ (mp[2] - mp[0]));
						dfdv[verts[1]][0] += 2.0 * mult * normSigma * delta * (Vec3d(1, 0, 0) ^ (mp[2] - mp[0]));
						dfdv[verts[1]][1] += 2.0 * mult * normSigma * delta * (Vec3d(0, 1, 0) ^ (mp[2] - mp[0]));
						dfdv[verts[1]][2] += 2.0 * mult * normSigma * delta * (Vec3d(0, 0, 1) ^ (mp[2] - mp[0]));
						dfdv[verts[2]][0] += 2.0 * mult * normSigma * delta * ((mp[1] - mp[0]) ^ Vec3d(1, 0, 0));
						dfdv[verts[2]][1] += 2.0 * mult * normSigma * delta * ((mp[1] - mp[0]) ^ Vec3d(0, 1, 0));
						dfdv[verts[2]][2] += 2.0 * mult * normSigma * delta * ((mp[1] - mp[0]) ^ Vec3d(0, 0, 1));
			*/			
						dfdv[verts[0]] += 2.0 * mult * normSigma * Vec3d(
							delta * Vec3d(0, -mp[1][2] + mp[2][2], mp[1][1] - mp[2][1]),
							delta * Vec3d(mp[1][2] - mp[2][2], 0, -mp[1][0] + mp[2][0]),
							delta * Vec3d(-mp[1][1] + mp[2][1], mp[1][0] - mp[2][0], 0));
						dfdv[verts[1]] += 2.0 * mult * normSigma * Vec3d(
							delta * Vec3d(0, mp[0][2] - mp[2][2], mp[2][1] - mp[0][1]),
							delta * Vec3d(mp[2][2] - mp[0][2], 0, mp[0][0] - mp[2][0]),
							delta * Vec3d(mp[0][1] - mp[2][1], mp[2][0] - mp[0][0], 0));
						dfdv[verts[2]] += 2.0 * mult * normSigma * Vec3d(
							delta * Vec3d(0, -mp[0][2] + mp[1][2], mp[0][1] - mp[1][1]),
							delta * Vec3d(mp[0][2] - mp[1][2], 0, -mp[0][0] + mp[1][0]),
							delta * Vec3d(-mp[0][1] + mp[1][1], mp[0][0] - mp[1][0], 0));
					}
#endif
				}
			}

			// marker springs
			int ind;
			for (ind = 0; ind < mSpringWeights.size(); ind++) {
				int joint = mSpringTrans[ind];
				SkelTransform *curTrans = skin->skel->transforms.getT(joint);
				double w = mSpringWeights[ind];
				if (w == 0)
					continue;
				Vec3d mPt;
				if (mSpringVerts.size() > 0)
					mPt = skin->curPts[mSpringVerts[ind]];
				else {
					// if there is spring information just for the 'a' poses, then 
					// we need to index things differently (and skip the non-'a' poses)
					if (mSpringNumEx == dataSet->numCharacters) {
						if (dataSet->charIndex[dataSet->examples[ex].character] == ex)
							mPt = mSpringVerts2[dataSet->examples[ex].character * mSpringWeights.size() + ind];
						else
							break;
					}
					else {
						mPt = mSpringVerts2[ex * mSpringWeights.size() + ind];
					}
				}
				Vec3d bPt = curTrans->globalCoord.v;

				if (!mPt.iszero()) {
					delta = mPt - bPt;
					markerSpringErr += w * delta.length2();
					numMarkerSpring += 3;

					w *= markerSpringSigma;
					if (mSpringVerts.size() > 0) {
						dfdv[mSpringVerts[ind]] += 2.0 * w * delta;
					}

					int dof = 0;
					int ip = 0;
					int pp = dpPose + ex * skin->skel->numPoseDofs;
					for (dof = 0; dof < skin->skel->numDofs; dof++) {
						if (skin->skel->dofIntrins[dof]) {
							for (comp=minComp; comp < maxComp; comp++) {
								if (optInt) {
									grad[dpInt + ip + comp * skin->skel->numIntrinsicDofs] -= 2.0 * w * 
										dataSet->charMu[dataSet->examples[ex].character][comp] * 
										delta * (curTrans->globalDerivs[dof] * Vec3d() 
	//									- vec4to3(skin->curMats[joint] * skin->curJointDerivs[joint*skin->skel->numDofs + dof])
										);
								}
								if (optX && comp > 0) {
									ofs = dpX + dataSet->examples[ex].character * (numComponents - 1);
									grad[ofs + comp - 1] -= 2.0 * w * 
										vInt[comp * skin->skel->numIntrinsicDofs + ip] * 
										delta * (curTrans->globalDerivs[dof] * Vec3d()
	//									- vec4to3(skin->curMats[joint] * skin->curJointDerivs[joint*skin->skel->numDofs + dof])
										);
								}
							}
							ip++;
						}
						else {
							if (optPose) {
								grad[pp++] -= 2.0 * w * delta * 
									(curTrans->globalDerivs[dof] * Vec3d());
							}
						}
					}
				}
			}

			applyDfdv(ex);
		}

		if (optNM || optNMP) { // || optX) {
			// normal map error
			for (pt=0; pt < NORMAL_MAP_SIZE; pt++) {
				int x = pt % NORMAL_MAP_W;
				int y = pt / NORMAL_MAP_W;
				int gradOfs = dpNM + pt * 3;
				int yMult = 1;
				if (pt / NORMAL_MAP_W >= NORMAL_MAP_H / 2) {
					gradOfs = dpNM + ((NORMAL_MAP_H - 1 - y) * NORMAL_MAP_W + x) * 3;
					yMult = -1;
				}
				
				if (!dataSet->examples[ex].normals[pt].iszero()) {
					calcNMErr(normalMap[pt], dataSet->examples[ex].normals[pt], ex, pt, gradOfs, yMult);
					
//					if (optX) {
//						for (comp = 1; comp < numComponents; comp++) {
//							grad[dpX + dataSet->examples[ex].character * (numComponents - 1) + comp-1] +=
//								nmSigma * 2.0 * delta * 
//						}
//					}
				}
			}
		}

		uiWait();
	} // examples loop


	evaluateWeightPrior();
	evaluateNMPrior();
	evaluatePDDPrior();
	evaluatePDNPrior();

/*
	if (optPDD) {
		// encourage PDDs to be small
		for (var = 0; var < vPDD.size(); var++) {
			// cauchy term
			ret += pddSigma * log(1.0 + sqr(100.0*vPDD[var]));
			grad[dpPDD + var] += pddSigma * (1.0/(1.0+100.0*vPDD[var])) * 2.0 * 100.0 * vPDD[var];

//			ret += pddSigma * sqr(vPDD[var]);
//			grad[dpPDD + var] += 2.0 * vPDD[var];
		}
	}*/

/*	if (optPose) {
		// small hack: encourage all joints to be zero in minPoses ("d" poses)
		for (ex = 0; ex < dataSet->numExamples; ex++) {
			if (!dataSet->examples[ex].minPose)
				continue;

			ofs = skin->skel->numPoseDofs * ex;
			for (tr=0; tr < skin->skel->transforms.size(); tr++) {
				SkelTransform *curTransform = skin->skel->transforms.getT(tr);
				if (!curTransform->isIntrinsic) {
					if (strcmp(curTransform->className, "SkelEulerRotation") == 0) {
						if (strncmp(curTransform->name+1, "Elbow", 5) == 0) {
							ret += 100*minPoseSigma * sqr(vPose[ofs]);
							grad[dpPose + ofs] += 2.0 * 100*minPoseSigma * vPose[ofs];
						}
						else {
							ret += minPoseSigma * sqr(vPose[ofs]);
							grad[dpPose + ofs] += 2.0 * minPoseSigma * vPose[ofs];
						}
					}
					else if (strcmp(curTransform->className, "SkelQuatRotation") == 0) {
						ret += minPoseSigma * sqr(vPose[ofs+0]);
						grad[dpPose + ofs+0] += 2.0 * minPoseSigma * vPose[ofs+0];
						ret += minPoseSigma * sqr(vPose[ofs+1]);
						grad[dpPose + ofs+1] += 2.0 * minPoseSigma * vPose[ofs+1];
						ret += minPoseSigma * sqr(vPose[ofs+2]);
						grad[dpPose + ofs+2] += 2.0 * minPoseSigma * vPose[ofs+2];
						ret += minPoseSigma * sqr(vPose[ofs+3] - 1.0);
						grad[dpPose + ofs+3] += 2.0 * minPoseSigma * (vPose[ofs+3] - 1.0);
					}
					ofs += curTransform->numDofs();
				}
			}
		}
	}*/

	// prior on X
	if (optX) {
		int ind = dpX;
		for (ch = 0; ch < dataSet->numCharacters; ch++) {
			for (comp = max(1, minComp); comp < maxComp; comp++) {
				ret += sqr(variables[ind]);
				grad[ind] += 2.0 * variables[ind];
				ind++;
			}
		}
	}

	// for marker files, make symmetry a soft constraint
	if (!uMesh && optDress) {
		for (comp=minComp; comp < maxComp; comp++) {
			double symSigma = 10000;
			for (pt=0; pt < numPts; pt++) {
				if (mirrorMap[pt] == -1)
					continue;
				Vec3d v0 = varDressPt(pt, comp);
				Vec3d v1 = varDressPt(mirrorMap[pt], comp);
				v1[1] = -v1[1];
				ret += symSigma * (v0 - v1).length2();
				if (pt != mirrorMap[pt]) {
					grad[dpDress+(comp*numPts + pt)*3 + 0] += 2.0 * symSigma * (v0[0] - v1[0]);
					grad[dpDress+(comp*numPts + mirrorMap[pt])*3 + 0] -= 2.0 * symSigma * (v0[0] - v1[0]);
					grad[dpDress+(comp*numPts + pt)*3 + 2] += 2.0 * symSigma * (v0[2] - v1[2]);
					grad[dpDress+(comp*numPts + mirrorMap[pt])*3 + 2] -= 2.0 * symSigma * (v0[2] - v1[2]);
				}
				grad[dpDress+(comp*numPts + pt)*3 + 1] += 2.0 * symSigma * (v0[1] - v1[1]);
				grad[dpDress+(comp*numPts + mirrorMap[pt])*3 + 1] += 2.0 * symSigma * (v0[1] - v1[1]);
			}
		}
	}

	// scaling
	if (optPose)
		for (var = 0; var < vPose.size(); var++)
			grad[dpPose + var] *= POSE_DOF_SCALE;
	if (optWeight)
		for (var = 0; var < vWeight.size(); var++)
			grad[dpWeight + var] *= WEIGHT_DOF_SCALE;

	ret += vertSigma * vertErr;
	ret += normSigma * normErr;
	ret += weightNeighSigma * weightNeighErr;
	ret += nmNeighSigma * nmNeighErr;
	ret += pddNeighSigma * pddNeighErr;
	ret += pdnNeighSigma * pdnNeighErr;
	ret += dressMuNeighSigma * dressMuNeighErr;
	ret += dressWNeighSigma * dressWNeighErr;
	ret += markerSpringSigma * markerSpringErr;
	ret += nmSigma * nmErr;

	lastErr = ret;

//	cout << lastErr << " (" << markerSpringErr << ", " << pddNeighErr << ")" << endl;

	uShow();

	redrawV();
	uiWait();

	return ret;
}

void USolver::evaluateGradient(Vecd& variables, Vecd& gradient) {
	if (optPose && globalPoseOnly) {
		int i, j;
		for (i = 0; i < dataSet->numExamples; i++) {
			for (j = 7; j < skin->skel->numPoseDofs; j++)
			grad[dpPose + skin->skel->numPoseDofs*i + j] = 0;
		}
	}
	if (optPose && lockPoseZero) {
		int ex, j;
		for (ex = 0; ex < dataSet->numExamples; ex++) {
			if (dataSet->charIndex[dataSet->examples[ex].character] == ex) {
				for (j = 0; j < skin->skel->numPoseDofs; j++) {
					grad[dpPose + ex*skin->skel->numPoseDofs + j] = 0;
				}
			}
		}
	}

	gradient = grad;
}

void USolver::evaluateWeightPrior() {
	int pt, i;

	if (!optWeight || !uMesh)
		return;

	for (pt = 0; pt < skin->numPts * skin->maxInf; pt++) {
		double err = 0;

		for (i=0; i < weightPriorVerts[pt].size(); i++) {
			err += weightPriorCoeffs[pt][i] * skin->infWeights[weightPriorVerts[pt][i]];
		}

		if (err == 0)
			continue;

		weightNeighErr += sqr(err);
		numWeightNeigh++;

		// calculate gradient
		for (i=0; i < weightPriorVerts[pt].size(); i++) {
			int gPt = weightPriorVerts[pt][i] / skin->maxInf;
			if (gPt >= numOrigPts)
				gPt = mirrorMap[gPt];

			int inf = weightPriorVerts[pt][i] % skin->maxInf;
			int ofs = dpWeight + gPt * (skin->maxInf - 1);
			if (inf == 0) {
				for (inf=0; inf < skin->maxInf-1; inf++) {
					if (skin->infJoints[gPt * skin->maxInf + inf + 1] >= 0)
						grad[ofs + inf] -= 2.0 * err * weightNeighSigma * weightPriorCoeffs[pt][i];
				}
			}
			else {
				grad[ofs + inf - 1] += 2.0 * err * weightNeighSigma * weightPriorCoeffs[pt][i];
			}
		}
	}
/*
	// neighbor-based regularization (ni = neighbor index)
	for (ni = 0; ni < neighTable.size(); ni++) {
		double delta;
		int v0Ofs, v1Ofs;
		NeighborRelation &nr = neighTable[ni];

		v0Ofs = dpWeight + nr.v0*(skin->maxInf-1);
		v1Ofs = dpWeight + nr.v1*(skin->maxInf-1);
		if (nr.v0 >= numOrigPts)
			v0Ofs = dpWeight + mirrorMap[nr.v0]*(skin->maxInf-1);
		if (nr.v1 >= numOrigPts)
			v1Ofs = dpWeight + mirrorMap[nr.v1]*(skin->maxInf-1);

		// weights
		delta = skin->infWeights[nr.v0*skin->maxInf + nr.ind0];
		if (nr.v1 > -1) {
			delta -= skin->infWeights[nr.v1*skin->maxInf + nr.ind1];
		}
		weightNeighErr += nr.dist * sqr(delta);
		numWeightNeigh += nr.dist;

		if (!_finite(weightNeighErr)) {
			cout << "bad values for gradient regularization " << nr.v0 << ", " << nr.v1 << endl;
			exit(0);
		}
		// update the gradients (note that index zero must be treated specially,
		// since it's 1 - the sum of the others)
		if (nr.ind0 == 0) {
			for (inf=0; inf < skin->maxInf-1; inf++) {
				if (skin->infJoints[nr.v0*skin->maxInf + inf + 1] >= 0)
					grad[v0Ofs + inf] -= 2.0 * nr.dist * weightNeighSigma * delta;
			}
		}
		else {
			grad[v0Ofs + nr.ind0-1] += 2.0 * nr.dist * weightNeighSigma * delta;
		}
		if (nr.v1 > -1) {
			if (nr.ind1 == 0) {
				for (inf=0; inf < skin->maxInf-1; inf++) {
					if (skin->infJoints[nr.v1*skin->maxInf + inf + 1] >= 0)
						grad[v1Ofs + inf] += 2.0 * nr.dist * weightNeighSigma * delta;
				}
			}
			else {
				grad[v1Ofs + nr.ind1-1] -= 2.0 * nr.dist * weightNeighSigma * delta;
			}
		}
	}*/
}

void USolver::evaluateNMPrior() {
	int comp, nmX, nmY;
	int ind;

	if (!optNM || !uMesh)
		return;

	for (comp = minComp; comp < maxComp; comp++) {
		for (nmY = 0; nmY < NORMAL_MAP_H/2; nmY++) {
			ind = comp * NORMAL_MAP_SIZE / 2 + nmY * NORMAL_MAP_W;
			for (nmX = 0; nmX < NORMAL_MAP_W; nmX++) {
				Vec3d v = vNM[ind + nmX];
				Vec3d delta;
				int ind2;

				if (nmY < NORMAL_MAP_H/2-1) {
					ind2 = ind + nmX + NORMAL_MAP_W;

/*					// Fisher distribution
					double cLen = v.length();
					double dLen = vNM[ind2].length();
					double delta = (v * vNM[ind2]) / (cLen * dLen);
					nmNeighErr += delta;
					numNMNeigh++;
					Vec3d deriv(
						(cLen * dLen) * vNM[ind2][0] - delta * dLen * dLen * v[0],
						(cLen * dLen) * vNM[ind2][1] - delta * dLen * dLen * v[1],
						(cLen * dLen) * vNM[ind2][2] - delta * dLen * dLen * v[2]);
					deriv = deriv / (cLen * cLen * dLen * dLen);
					grad[dpNM + (ind+nmX)*3+0] += nmNeighSigma * deriv[0];
					grad[dpNM + (ind+nmX)*3+1] += nmNeighSigma * deriv[1];
					grad[dpNM + (ind+nmX)*3+2] += nmNeighSigma * deriv[2];
					deriv(
						(cLen * dLen) * v[0] - delta * cLen * cLen * vNM[ind2][0],
						(cLen * dLen) * v[1] - delta * cLen * cLen * vNM[ind2][1],
						(cLen * dLen) * v[2] - delta * cLen * cLen * vNM[ind2][2]);
					deriv = deriv / (cLen * cLen * dLen * dLen);
					grad[dpNM + ind2 * 3 + 0] += nmNeighSigma * deriv[0];
					grad[dpNM + ind2 * 3 + 1] += nmNeighSigma * deriv[1];
					grad[dpNM + ind2 * 3 + 2] += nmNeighSigma * deriv[2];
*/
					
					// Gaussian distribution
					delta = v - vNM[ind2];
					nmNeighErr += delta.length2();
					numNMNeigh += 3;

					grad[dpNM + (ind+nmX)*3+0] += nmNeighSigma * 2.0 * delta[0];
					grad[dpNM + (ind+nmX)*3+1] += nmNeighSigma * 2.0 * delta[1];
					grad[dpNM + (ind+nmX)*3+2] += nmNeighSigma * 2.0 * delta[2];
					grad[dpNM + ind2 * 3 + 0] -= nmNeighSigma * 2.0 * delta[0];
					grad[dpNM + ind2 * 3 + 1] -= nmNeighSigma * 2.0 * delta[1];
					grad[dpNM + ind2 * 3 + 2] -= nmNeighSigma * 2.0 * delta[2];
					
				}
				if (nmX < NORMAL_MAP_W-1) {
					ind2 = ind + nmX + 1;

/*					// Fisher distribution
					double cLen = v.length();
					double dLen = vNM[ind2].length();
					double delta = (v * vNM[ind2]) / (cLen * dLen);
					nmNeighErr += delta;
					numNMNeigh++;
					Vec3d deriv(
						(cLen * dLen) * vNM[ind2][0] - delta * dLen * dLen * v[0],
						(cLen * dLen) * vNM[ind2][1] - delta * dLen * dLen * v[1],
						(cLen * dLen) * vNM[ind2][2] - delta * dLen * dLen * v[2]);
					deriv = deriv / (cLen * cLen * dLen * dLen);
					grad[dpNM + (ind+nmX)*3+0] += nmNeighSigma * deriv[0];
					grad[dpNM + (ind+nmX)*3+1] += nmNeighSigma * deriv[1];
					grad[dpNM + (ind+nmX)*3+2] += nmNeighSigma * deriv[2];
					deriv(
						(cLen * dLen) * v[0] - delta * cLen * cLen * vNM[ind2][0],
						(cLen * dLen) * v[1] - delta * cLen * cLen * vNM[ind2][1],
						(cLen * dLen) * v[2] - delta * cLen * cLen * vNM[ind2][2]);
					deriv = deriv / (cLen * cLen * dLen * dLen);
					grad[dpNM + ind2 * 3 + 0] += nmNeighSigma * deriv[0];
					grad[dpNM + ind2 * 3 + 1] += nmNeighSigma * deriv[1];
					grad[dpNM + ind2 * 3 + 2] += nmNeighSigma * deriv[2];
*/

					// Gaussian distribution
					delta = v - vNM[ind2];
					nmNeighErr += delta.length2();

					grad[dpNM + (ind+nmX)*3+0] += nmNeighSigma * 2.0 * delta[0];
					grad[dpNM + (ind+nmX)*3+1] += nmNeighSigma * 2.0 * delta[1];
					grad[dpNM + (ind+nmX)*3+2] += nmNeighSigma * 2.0 * delta[2];
					grad[dpNM + ind2 * 3 + 0] -= nmNeighSigma * 2.0 * delta[0];
					grad[dpNM + ind2 * 3 + 1] -= nmNeighSigma * 2.0 * delta[1];
					grad[dpNM + ind2 * 3 + 2] -= nmNeighSigma * 2.0 * delta[2];
				}
			}
		}
	}
}

void USolver::evaluatePDDPrior() {
	int ni, comp;

	if (!optPDD)
		return;

	for (ni = 0; ni < pddNeighTable.size(); ni++) {
		Vec3d delta;
		int ofs0, ofs1;
		NeighborRelation &nr = pddNeighTable[ni];
//		int sign0 = (nr.v0 >= numOrigPts)?-1:1;
//		int sign1 = (nr.v1 >= numOrigPts)?-1:1;

		for (comp = minComp; comp < maxPComp; comp++) {
			ofs0 = (comp * skin->numPddPtKeys + nr.ind0) * 3;
			ofs1 = (comp * skin->numPddPtKeys + nr.ind1) * 3;
			if (nr.ind1 < 0)
				ofs1 = -1;
			delta = Vec3d(vPDD[ofs0 + 0], vPDD[ofs0 + 1], vPDD[ofs0 + 2]);
			if (ofs1 >= 0)
				delta -= Vec3d(vPDD[ofs1 + 0], vPDD[ofs1 + 1], vPDD[ofs1 + 2]);

			pddNeighErr += nr.dist * delta.length2();
			numPddNeigh += nr.dist;

			grad[dpPDD + ofs0 + 0] += 2.0 * nr.dist * pddNeighSigma * delta[0];
			grad[dpPDD + ofs0 + 1] += 2.0 * nr.dist * pddNeighSigma * delta[1];
			grad[dpPDD + ofs0 + 2] += 2.0 * nr.dist * pddNeighSigma * delta[2];
			if (ofs1 >= 0) {
				grad[dpPDD + ofs1 + 0] -= 2.0 * nr.dist * pddNeighSigma * delta[0];
				grad[dpPDD + ofs1 + 1] -= 2.0 * nr.dist * pddNeighSigma * delta[1];
				grad[dpPDD + ofs1 + 2] -= 2.0 * nr.dist * pddNeighSigma * delta[2];
			}
		}
	}

/*		for (pt=0; pt < numOrigPts; pt++) {
			// dressPose regularization
			double dressNeighWeight;
			for (ni=0; ni < neighbors[pt].size(); ni++) {
				if (neighbors[pt][ni] >= numOrigPts)
					continue;
				for (comp = 0; comp < numComponents; comp++) {
					int ofs0 = dressDofPos + (comp * numOrigPts + pt) * 3;
					int ofs1 = dressDofPos + (comp * numOrigPts + neighbors[pt][ni]) * 3;
					double delta;
					int i;
					for (i=0; i < 3; i++) {
						delta = vars[ofs0+i] - vars[ofs1+i];
						if (comp == 0) {
	                        dressMuNeighErr += sqr(delta);
							numDressMuNeigh++;
							grad[ofs0+i] += dressMuNeighSigma * 2.0 * delta;
							grad[ofs1+i] -= dressMuNeighSigma * 2.0 * delta;
						}
						else {
	                        dressWNeighErr += sqr(delta);
							numDressWNeigh++;
							grad[ofs0+i] += dressWNeighSigma * 2.0 * delta;
							grad[ofs1+i] -= dressWNeighSigma * 2.0 * delta;
						}
					}
				}
			}
		}*/
}

void USolver::evaluatePDNPrior() {
	int ni, comp;

	if (!optNMP)
		return;

	for (ni = 0; ni < pdnNeighTable.size(); ni++) {
		Vec3d delta;
		int ofs0, ofs1;
		NeighborRelation &nr = pdnNeighTable[ni];

		for (comp = minComp; comp < maxPComp; comp++) {
			ofs0 = (comp * skin->numPddNMKeys + nr.ind0) * 3;
			ofs1 = (comp * skin->numPddNMKeys + nr.ind1) * 3;
			if (nr.ind1 < 0)
				ofs1 = -1;
			delta = Vec3d(vNMP[ofs0 + 0], vNMP[ofs0 + 1], vNMP[ofs0 + 2]);
			if (ofs1 >= 0)
				delta -= Vec3d(vNMP[ofs1 + 0], vNMP[ofs1 + 1], vNMP[ofs1 + 2]);

			pdnNeighErr += delta.length2();
			numPdnNeigh += 3;

			grad[dpNMP + ofs0 + 0] += 2.0 * pdnNeighSigma * delta[0];
			grad[dpNMP + ofs0 + 1] += 2.0 * pdnNeighSigma * delta[1];
			grad[dpNMP + ofs0 + 2] += 2.0 * pdnNeighSigma * delta[2];
			if (ofs1 >= 0) {
				grad[dpNMP + ofs1 + 0] -= 2.0 * pdnNeighSigma * delta[0];
				grad[dpNMP + ofs1 + 1] -= 2.0 * pdnNeighSigma * delta[1];
				grad[dpNMP + ofs1 + 2] -= 2.0 * pdnNeighSigma * delta[2];
			}
		}
	}
}

void USolver::solverStep() {
//	cout << "---" << endl;
	cout << lastErr << " (" << markerSpringErr << ", " << pddNeighErr << ")" << endl;
}

void USolver::renormalizeX() {
	int comp, ch, var, n;
	double *avgX = new double[numComponents];
	memset(avgX, 0, sizeof(double) * numComponents);

	// first, set the average x to be zero
	for (ch = 0; ch < dataSet->numCharacters; ch++) {
		for (comp = 0; comp < numComponents; comp++) {
			avgX[comp] += dataSet->charMu[ch][comp];
		}
	}
	for (comp = 0; comp < numComponents; comp++) {
		avgX[comp] /= dataSet->numCharacters;
		cout << avgX[comp] << " ";
	}
	cout << endl;
	for (ch = 0; ch < dataSet->numCharacters; ch++) {
		for (comp = 1; comp < numComponents; comp++) {
			dataSet->charMu[ch][comp] -= avgX[comp];
		}
	}
	for (comp = 1; comp < numComponents; comp++) {
		// fix dress points
		n = numOrigPts * 3;
		for (var = 0; var < n; var++) {
			vDress[var] += avgX[comp] * vDress[comp * n + var];
		}
		// fix bone intrinsics
		n = skin->skel->numIntrinsicDofs;
		for (var = 0; var < n; var++) {
			vInt[var] += avgX[comp] * vInt[comp * n + var];
		}
		// fix normal maps
		if (uMesh) {
			n = NORMAL_MAP_SIZE / 2;
			for (var = 0; var < n; var++) {
				vNM[var] += avgX[comp] * vNM[comp * n + var];
			}
		}
		if (comp < numPDComponents) {
			// fix PDDs
			n = vPDD.size() / numPDComponents;
			for (var = 0; var < n; var++) {
				vPDD[var] += avgX[comp] * vPDD[comp * n + var];
			}
			// fix PDNs
			n = vNMP.size() / numPDComponents;
			for (var = 0; var < n; var++) {
				vNMP[var] += avgX[comp] * vNMP[comp * n + var];
			}
		}
	}

	// now set the variance to be 1
	for (comp = 1; comp < numComponents; comp++) {
		double sigma = 0;
		for (ch = 0; ch < dataSet->numCharacters; ch++) {
			sigma += sqr(dataSet->charMu[ch][comp]);
		}
		sigma = sqrt(sigma / (dataSet->numCharacters - 1));
		cout << sigma << " ";

		for (ch = 0; ch < dataSet->numCharacters; ch++) {
			dataSet->charMu[ch][comp] /= sigma;
		}
		// fix dress points
		n = numOrigPts * 3;
		for (var = 0; var < n; var++) {
			vDress[comp * n + var] *= sigma;
		}
		// fix bone intrinsics
		n = skin->skel->numIntrinsicDofs;
		for (var = 0; var < n; var++) {
			vInt[comp * n + var] *= sigma;
		}
		// fix normal maps
		if (uMesh) {
			n = NORMAL_MAP_SIZE / 2;
			for (var = 0; var < n; var++) {
				vNM[comp * n + var] *= sigma;
			}
		}
		if (comp < numPDComponents) {
			// fix PDDs
			n = vPDD.size() / numPDComponents;
			for (var = 0; var < n; var++) {
				vPDD[comp * n + var] *= sigma;
			}
			// fix PDNs
			n = vNMP.size() / numPDComponents;
			for (var = 0; var < n; var++) {
				vNMP[comp * n + var] *= sigma;
			}
		}
	}
	cout << endl;

	delete []avgX;
}

void USolver::autoWeight() {
	// update sigmas
	if (numVert > 0) {
		cout << "vert sigma: " << vertSigma << " -> ";
		vertSigma = numVert / vertErr;
		cout << vertSigma << endl;
	}
	if (numNorm > 0) {
		cout << "norm sigma: " << normSigma << " -> ";
		normSigma = numNorm / normErr;
		cout << normSigma << endl;
	}
	if (numWeightNeigh > 0) {
		cout << "weight neigh: " << weightNeighSigma << " -> ";
		weightNeighSigma = /*0.01 */ numWeightNeigh / weightNeighErr;
		cout << weightNeighSigma << endl;
	}
	if (numMarkerSpring > 0) {
		cout << "marker springs: " << markerSpringSigma << " -> ";
		markerSpringSigma = 237 * numMarkerSpring / markerSpringErr;
		cout << markerSpringSigma << endl;
	}
	if (numPddNeigh > 0) {
		cout << "PDD neigh: " << pddNeighSigma << " -> ";
		pddNeighSigma = numPddNeigh / pddNeighErr;
		cout << pddNeighSigma << endl;
	}
/*	if (numNM > 0) {
		cout << "Normal maps: " << nmSigma << " -> ";
		nmSigma = numNM / nmErr;
		cout << nmSigma << endl;
	}*/
	if (numNMNeigh > 0) {
		cout << "Normal map neigh: " << nmNeighSigma << " -> ";
		nmNeighSigma = /*0.01 */ numNMNeigh / nmNeighErr;
		cout << nmNeighSigma << endl;
	}
	if (numPdnNeigh > 0) {
		cout << "PDN neigh: " << pdnNeighSigma << " -> ";
		pdnNeighSigma = numPdnNeigh / pdnNeighErr;
		cout << pdnNeighSigma << endl;
	}
/*	if (uSolver.numDressMuNeigh > 0) {
		uSolver.dressMuNeighSigma = uSolver.numDressMuNeigh / uSolver.dressMuNeighErr;
	}
	if (uSolver.numDressWNeigh > 0) {
		uSolver.dressWNeighSigma = uSolver.numDressWNeigh / uSolver.dressWNeighErr;
	}
	cout << "s_v = " << uSolver.vertSigma << ", ";
	cout << "s_l = " << uSolver.lamdaSigma << ", ";
/*	cout << "s_w = " << uSolver.weightNeighSigma << ", ";
	cout << "s_p = " << uSolver.pddNeighSigma << ", ";
	cout << "s_dm = " << uSolver.dressMuNeighSigma << ", ";
	cout << "s_dw = " << uSolver.dressWNeighSigma << endl;
*/
}

void USolver::eStep() {
	int ex, pt, ma;
	int *maCount = new int[numMarkers * numPts];
	memset(markerAssignments, 0, sizeof(double) * numMarkers * numPts);
	memset(maCount, 0, sizeof(int) * numMarkers * numPts);

	for (ex=0; ex < dataSet->numExamples; ex++) {
		memset(dfdv, 0, sizeof(Vec3d) * numPts);

		updateSkel(ex);
		updatePoints(dataSet->examples[ex].character);
		skin->skel->updateDerivs();
		skin->skel->updateGlobalDerivs();
		skin->updateMats();
		skin->updateJoints();	// could check if this is 1st or not
		skin->updatePts(0, -1, optNM || optNMP);

		for (pt = 0; pt < numPts; pt++) {
			for (ma = 0; ma < numMarkers; ma++) {
				if (dataSet->examples[ex].ptsConf[ma] == 0)
					continue;

				Vec3d exPt = dataSet->examples[ex].pts[ma];
				markerAssignments[ma*numPts + pt] += 
					exp(-(skin->curPts[pt] - exPt).length2() / (2.0 * sqr(0.1)));
			}
		}
	}
	for (ma = 0; ma < numMarkers; ma++) {
		double sum = 0;
		for (pt = 0; pt < numPts; pt++) {
			sum += markerAssignments[ma*numPts + pt];
		}
		for (pt = 0; pt < numPts; pt++) {
			markerAssignments[ma*numPts + pt] /= sum;
//			if (ma == 0 && markerAssignments[ma*numPts + pt] > 0.001)
//				cout << pt << " " << markerAssignments[ma*numPts + pt] << endl;
		}
	}

	delete []maCount;
#ifdef OLD_VERSION
	// ex = example (scan) index
	int ex;
	// pt = point index
	int pt;
	// inf = influence iterator
	int inf;
	// comp = component iterator
	int comp;
	// dof = dof iterator
	int dof;
	// y_ij = vertex with index j from example i
	Vec3d y_ij;
	// ii, jj, kk are used to iterate over vectors and matrices
	int ii, jj, kk;
	// current and last character id
	int charID = -1, lastID = -1;
	Vec3d *pddComps = new Vec3d[numComponents];

	int numPts = 0;
	dataSet->phiLogEntropy = 0;

	// static data
	static bool doneInit = false;
	static VLVecd muAcc;
	static Vec3d *aAcc;

	// initialize static data
	if (!doneInit) {
		muAcc.SetSize(numComponents);
		aAcc = new Vec3d[numComponents];
		doneInit = true;
	}

	charID = dataSet->examples[0].character;
	for (ex = 0; ex < dataSet->numExamples; ex++) {
		if (maxSolveEx >= 0 && ex >= maxSolveEx)
			break;

		// calculate the current reconstruction
		updateSkel(ex);
		updatePoints(dataSet->examples[ex].character);
		skin->updateMats();
		skin->updateJoints();	// could check if this is the first or not
		skin->updatePts();

		if (lastID != charID) {
			// clear the current average and covariance estimates
			dataSet->charMu[charID].MakeZero();
			dataSet->charMu[charID][0] = 1;
//			dataSet->charPhi[charID].MakeDiag(sigma2);
			dataSet->charPhi[charID].MakeZero();
			dataSet->charPhi[charID][0][0] = 1;

			muAcc.MakeZero();

			numPts = 0;
		}

		if (numComponents < 2)
			continue;

		// iterate over all observed points
		double maxPts = numPts;
		if (!dataSet->examples[ex].lookup)
			maxPts = min(dataSet->examples[ex].numPts, maxPts);
		for (pt = 0; pt < maxPts; pt++) {
			double mult;
			Vec3d exPt;
			if (!dataSet->examples[ex].getPt(pt, &exPt, &mult) || mult == 0)
				continue;
			mult *= vertSigma;

			memset(pddComps, 0, sizeof(Vec3d)*numComponents);
/*			if (pddDofPos >= 0) {
				for (inf = 0; inf < skin->maxInf; inf++) {
					double weights[KNN_MAX_SAMPLES];
					int joint = skin->infJoints[pt*skin->maxInf + inf];
					if (joint < 0)
						continue;
					int numSamples = skin->transPddOfs[joint+1]-skin->transPddOfs[joint];
					int ofs = pt*skin->maxInf + inf;
					knnQuatInterp(
							skin->skel->transforms.getT(joint)->curCoord.q, 
							skin->pddQuats + skin->transPddOfs[joint], numSamples,
							skin->pddKeys + skin->ptPddOfs[ofs], weights);
					for (comp = 0; comp < numComponents; comp++) {
						int samp;
						for (samp = 0; samp < numSamples; samp++) {
							double *pdd = &varPdd(comp, pt, inf, samp, 0);
							pddComps[comp] += weights[samp] * Vec3d(pdd[0], pdd[1], pdd[2]);
						}
					}
				}
			}*/

			for (comp=0; comp < numComponents; comp++)
				aAcc[comp].zeroElements();

			Vec3d trans;
			for (inf = 0; inf < skin->maxInf; inf++) {
				if (skin->infJoints[pt*skin->maxInf + inf] >= 0) {
					double w = skin->infWeights[pt*skin->maxInf + inf];
					int joint = skin->infJoints[pt*skin->maxInf + inf];
					
					for (comp = 0; comp < numComponents; comp++) {
						aAcc[comp] += w * vec4to3(skin->curMats[joint] * (vec3to4z(varDressPt(pt, comp) + pddComps[comp])));
					}
					trans += w * (skin->curMats[joint]*Vec3d() - vec4to3(skin->curMats[joint] * vec3to4z(skin->dressJoints[joint])));
				}
			}

			numPts++;

			// phi += A'A
			for (ii=1; ii < numComponents; ii++) {
				for (jj=1; jj < numComponents; jj++) {
					dataSet->charPhi[charID][ii][jj] += mult * aAcc[ii] * aAcc[jj];
					if (!_finite(dataSet->charPhi[charID][ii][jj]))
						cout << "bad";
				}
			}

			// muAcc += (t-b-trans)'A
			for (ii=1; ii < numComponents; ii++) {
				muAcc[ii] += mult * (exPt - aAcc[0] - trans) * aAcc[ii];
			}
		}

		lastID = charID;
		if (ex < dataSet->numExamples-1) {
			charID = dataSet->examples[ex+1].character;
		}
		if ((ex == dataSet->numExamples-1) || (lastID != charID)) {
//			cout << (1.0 * numPts / numPts) << endl;
//			dataSet->charPhi[lastID] *= (1.0 * numPts / numPts);
//			muAcc *= (1.0 * numPts / numPts);
			muAcc[0] = 1;
			dataSet->charPhi[lastID][0][0] = 1;
			for (ii=1; ii < numComponents; ii++) {
				dataSet->charPhi[lastID][ii][ii] += 1;
			}

			// actual phi is the inverse
/*			if (lastID == 0) {
				cout << "saving" << endl;
				ofstream mat("m.txt");
				mat << dataSet->charPhi[lastID];
				mat << endl;
				mat << muAcc;
				mat.close();
			}*/
			dataSet->charPhi[lastID] = inv(dataSet->charPhi[lastID]);
/*			if (lastID == 0) {
				cout << "saving" << endl;
				ofstream mat("m-inv.txt");
				mat << dataSet->charPhi[lastID];
				mat.close();
			}
*/
			dataSet->charMu[lastID] = dataSet->charPhi[lastID] * muAcc;
			dataSet->charMu[lastID][0] = 1;

			// update log entropy
			double det = 1;
			VLMatd U, V; VLVecd diag;
			SVDFactorization(dataSet->charPhi[lastID], U, V, diag);
			for (ii = 0; ii < numComponents; ii++)
				det *= diag[ii];
			dataSet->phiLogEntropy += log(det);

			// add mu'mu into curPhi, and set 1st row & column to mu
			for (ii=0; ii < numComponents; ii++) {
				for (jj=0; jj < numComponents; jj++) {
					if (ii == 0)
						dataSet->charPhi[lastID][ii][jj] = dataSet->charMu[lastID][jj];
					else if (jj == 0)
						dataSet->charPhi[lastID][ii][jj] = dataSet->charMu[lastID][ii];
					else {
						dataSet->charPhi[lastID][ii][jj] = 
							dataSet->charPhi[lastID][ii][jj] + dataSet->charMu[lastID][ii] * dataSet->charMu[lastID][jj];
					}
				}
			}
		}
	}

	delete []pddComps;
#endif
}


// obsolete stuff
/*
double USolver::evaluateFunctionExp() {
	int ex, pt, tr, ni, inf, samp, ch, dof;
	int comp;
	double ret = 0;

	updateWeights();

	// point-distance terms
	for (ex=0; ex < dataSet->numExamples; ex++) {
		updateSkel(ex);
		updatePoints(dataSet->examples[ex].character);
		skin->skel->updateDerivs();
		skin->skel->updateGlobalDerivs();
		skin->updateMats();
		skin->updateJoints();	// could check if this is the first or not
		skin->updatePts();

		for (pt = 0; pt < numPts; pt++) {
			// don't go past the number of points in this set
			if (!dataSet->examples[ex].numPts && pt > dataSet->examples[ex].numPts)
				break;

			Vec3d exPt;
			double curConf;
			if (dataSet->examples[ex].getPt(pt, &exPt, &curConf) && curConf > 0) {
				ret += vertSigma * skin->calcGradExp(pt, exPt, 
					NULL, NULL, NULL, NULL, NULL, 
					dataSet->charPhi[dataSet->examples[ex].character], this);
			}
		}
	}

	return ret;
}

void USolver::calcPtErrAndGrad(int pt, Vec3d &target, int ex, double mult) {
	int comp, inf, var, ofs;

	if (mult == 0)
		return;

	static double *weightGrad = NULL;
	if (!weightGrad)
		weightGrad = new double[skin->maxInf];
	Vec3d dressGrad;
	static int maxPddVars = skin->maxInf * 3 * 16;
	static double *pddGrad = new double[maxPddVars];
	static double *intGrad = new double[skin->skel->numIntrinsicDofs];

	Vec3d delta = skin->curPts[pt] - target;
	vertErr += mult * delta.length2();
	numVert += 3;

	dressGrad.zeroElements();
	memset(weightGrad, 0, sizeof(double)*skin->maxInf);
	if (USE_PDD_DOFS)
		memset(pddGrad, 0, sizeof(double)*maxPddVars);
	memset(intGrad, 0, sizeof(double)*skin->skel->numIntrinsicDofs);
	
	skin->calcGrad(pt, 2.0 * mult * delta, 
		optDress | optX ? dressGrad.n : NULL, 
		optInt | optX   ? intGrad     : NULL, 
		optPose         ? grad.n + dpPose + ex * skin->skel->numPoseDofs : NULL, 
		optWeight       ? weightGrad  : NULL, 
		optPDD | optX   ? pddGrad     : NULL);

	for (comp = 0; comp < numComponents; comp++) {
		double w = dataSet->charMu[dataSet->examples[ex].character][comp];

		if (w != 0) {
			if (optDress) {
				// fill in dress shape gradients
				if (pt >= numOrigPts) {
					ofs = dpDress + (comp * numOrigPts + mirrorMap[pt]) * 3;
					grad.n[ofs + 0] += w * dressGrad[0];
					grad.n[ofs + 1] -= w * dressGrad[1];
					grad.n[ofs + 2] += w * dressGrad[2];
				}
				else {
					ofs = dpDress + (comp * numOrigPts + pt) * 3;
					grad.n[ofs + 0] += w * dressGrad[0];
					grad.n[ofs + 1] += w * dressGrad[1];
					grad.n[ofs + 2] += w * dressGrad[2];
				}
			}

			if (optInt) {
				// fill in intrinsic DOF gradients
				ofs = dpInt + comp * skin->skel->numIntrinsicDofs;
				for (var = 0; var < skin->skel->numIntrinsicDofs; var++)
					grad[ofs + var] += w * intGrad[var];
			}
		}

		if (optX && comp > 0) {
			// fill in PCA reconstruction weight gradients
			ofs = dpX + dataSet->examples[ex].character * (numComponents - 1);
			grad[ofs + comp - 1] += varDressPt(pt, comp) * dressGrad;
			for (var = 0; var < skin->skel->numIntrinsicDofs; var++)
				grad[ofs + comp - 1] += vInt[comp * skin->skel->numIntrinsicDofs + var] * 
					intGrad[var];
		}
	}

	if (optWeight) {
		// fill in weight gradients
		if (pt >= numOrigPts)
			ofs = dpWeight + mirrorMap[pt]*(skin->maxInf - 1);
		else
			ofs = dpWeight + pt*(skin->maxInf - 1);
		for (inf = 1; inf < skin->maxInf; inf++) {
			if (skin->infJoints[pt*skin->maxInf + inf] >= 0) {
				grad.n[ofs + inf - 1] += weightGrad[inf];
				grad.n[ofs + inf - 1] -= weightGrad[0];
			}
		}
	}

	if (optPDD) {
		// fill in PDD gradients
		int varOfs;
		int numVars;
		int numSrcPddVars = skin->ptPddOfs[numOrigPts * skin->maxInf] * 3;
		if (pt >= numOrigPts) {
			varOfs = skin->ptPddOfs[mirrorMap[pt]*skin->maxInf];
			numVars  = skin->ptPddOfs[(mirrorMap[pt]+1)*skin->maxInf] - varOfs;
		}
		else {
			varOfs = skin->ptPddOfs[pt*skin->maxInf];
			numVars  = skin->ptPddOfs[(pt+1)*skin->maxInf] - varOfs;
		}
		if (numVars > 0) {
			for (comp = 0; comp < numComponents; comp++) {
				double w = dataSet->charMu[dataSet->examples[ex].character][comp];
				if (w != 0) {
					for (var = 0; var < numVars; var++) {
						Vec3d v;
						int ofs = dpPDD + comp*numSrcPddVars + (varOfs+var)*3;
						grad[ofs + 0] += w * pddGrad[var*3 + 0];
						if (pt >= numOrigPts)
							grad[ofs + 1] -= w * pddGrad[var*3 + 1];
						else
							grad[ofs + 1] += w * pddGrad[var*3 + 1];
						grad[ofs + 2] += w * pddGrad[var*3 + 2];
					}
				}
			}
		}
	}
}


// closest-point mesh error function
		if (dataSet->examples[ex].mesh && uMesh) {
			for (pt=0; pt < uMesh->numPts(); pt++) {
				uMesh->getPt(pt) = skin->curPts[pt];
			}
			uMesh->calcNormals();

			for (pt = 0; pt < numPts; pt++) {
				double curDistance, mult;
				Vec3d curClosestPt;
				rayPoint(skin->curPts[pt], uMesh->getPtNormal(pt), 
					&dataSet->examples[ex], curDistance, curClosestPt, mult);

				calcPtErrAndGrad(pt, curClosestPt, ex, mult * vertSigma);
			}
		}
		else {



void USolver::addPdd() {
	static double *errors = NULL;
	static double *eCount = NULL;
	int ex, pt, inf, tr, comp;
	int i;

	if (errors == NULL) {
		errors = new double[skin->numTrans * dataSet->numExamples];
		eCount = new double[skin->numTrans * dataSet->numExamples];
		memset(errors, 0, sizeof(double) * skin->numTrans * dataSet->numExamples);
		memset(eCount, 0, sizeof(double) * skin->numTrans * dataSet->numExamples);
		// don't use example 0
		for (i=0; i < skin->numTrans; i++)
			eCount[i] = -1;
		// don't allow influences for joint 2
		for (i=0; i < dataSet->numExamples; i++)
			eCount[i*skin->numTrans + 2] = -1;
	}
	else {
		for (i=0; i < skin->numTrans * dataSet->numExamples; i++) {
			if (eCount[i] > 0) eCount[i] = 0;
			errors[i] = -1;
		}
	}

	
//	for (ex=0; ex < dataSet->numExamples; ex++) {
//		updateSkel(ex);
//		updatePoints(dataSet->examples[ex].character);
//		skin->updateMats();
//		if (ex == 0) {
//			skin->updateJoints();
//		}
//		skin->updatePts();
//
//		for (pt = 0; pt < numPts; pt++) {
//			// don't go past the number of points in this set
//			if (pt > dataSet->examples[ex].numPts)
//				break;
//
//			// ignore points in the PDD mask
//			if (pddMask && pddMask[pt] == 0)
//				continue;
//
//			double mult = 1;
//			if (dataSet->examples[i].ptsConf)
//				mult = dataSet->examples[ex].ptsConf[pt];
//			if (mult > 0) {
//				Vec3d delta = skin->curPts[pt] - dataSet->examples[ex].pts[pt];
//				for (inf=0; inf < skin->maxInf; inf++) {
//					int curTrans = skin->infJoints[pt * skin->maxInf + inf];
//					if (curTrans >= 0 && eCount[ex * skin->numTrans + curTrans] >= 0) {
//						errors[ex * skin->numTrans + curTrans] += mult * delta.length2();
////						eCount[ex * skin->numTrans + curTrans]++;
//					}
//				}
//			}
//		}
//	}
	for (ex=0; ex < dataSet->numExamples; ex++) {
		updateSkel(ex);
		for (tr = 3; tr < skin->numTrans; tr++) {
			// hack: skip uncombined and half-joints
			if (tr == 6 || tr == 7 || tr == 10 || tr == 11 || tr == 12 ||
				tr == 24 || tr == 25 || tr == 28 || tr == 29 || tr == 30 ||
				tr == 19 || tr == 36)
				continue;

			QuatNorm q = skin->skel->transforms.getT(tr)->curCoord.q;

			int primaryTr = tr;
			if (mirrorTrans[tr] < tr) {
				primaryTr = mirrorTrans[tr];
				q.y *= -1;
				q.w *= -1;
			}

			double minDist = quatDist(q, QuatNorm(0,0,0,1));

			int numSamples = skin->transPddOfs[primaryTr+1]-skin->transPddOfs[primaryTr];
			for (inf=0; inf < numSamples; inf++) {
				Vec4d &v = skin->pddQuats[skin->transPddOfs[primaryTr]+inf];
				double dist = quatDist(q, QuatNorm(v[0],v[1],v[2],v[3]));
				if (dist < minDist) {
					minDist = dist;
				}
			}

			if (errors[ex * skin->numTrans + primaryTr] >= 0)
				errors[ex * skin->numTrans + primaryTr] = min(errors[ex * skin->numTrans + primaryTr], minDist);
			else
				errors[ex * skin->numTrans + primaryTr] = minDist;
		}
	}

	// find the maximum
	double maxV = 0;
	int maxEx = -1, maxTrans = 0;
	for (ex=0; ex < dataSet->numExamples; ex++) {
		for (tr=0; tr < skin->numTrans; tr++) {
			double v;
			if (eCount[ex * skin->numTrans + tr] >= 0) {
				v = errors[ex * skin->numTrans + tr]; // / eCount[ex * skin->numTrans + tr];
				// downweight based on # of PDDs already assigned to this joint
				//v /= 1.0 + skin->transPddOfs[tr+1] - skin->transPddOfs[tr];
				if (v > maxV) {
					maxV = v;
					maxEx = ex;
					maxTrans = tr;
				}
			}
		}
	}

	if (maxEx >= 0) {
		// save old PDDs
		int oldNumPddVars = skin->ptPddOfs[numOrigPts * skin->maxInf] * 3;
		double *oldPdds = new double[oldNumPddVars * numComponents];
		memcpy(oldPdds, vPDD.n, sizeof(double) * oldNumPddVars * numComponents);

		// add new PDD
		cout << "adding a key for example " << maxEx << ", joint " << skin->skel->transforms.getT(maxTrans)->name << ", error " << maxV << endl;
		updateSkel(maxEx);
		QuatNorm q = skin->skel->transforms.getT(maxTrans)->curCoord.q;
		q.normalize();
		skin->insertPdd(maxTrans, q);
		if (mirrorTrans[maxTrans] != maxTrans)
			skin->insertPdd(mirrorTrans[maxTrans], QuatNorm(q.x, -q.y, q.z, -q.w));
		cout << "quat = " << q << "; value = " << maxV << endl;

		// update number of variables
		vPDD.resize(skin->ptPddOfs[numOrigPts * skin->maxInf] * 3 * numComponents, true);

		// don't use this one again!
		eCount[maxEx * skin->numTrans + maxTrans] = -1;

		delete []oldPdds;
	}
}

*/

