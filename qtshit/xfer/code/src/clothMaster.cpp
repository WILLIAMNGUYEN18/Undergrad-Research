#include "doppel2.h"
#include "cli.h"
#include "clothMaster.h"
#include "cloth.h"
#include "cvMaster.h"
#include "cvgf.h"
#include "markers.h"

extern CharVec charVec;

Cloth *cloth = NULL;
vector<int> markerRefs;
MarkerSet *markers;

void ClothUI::drawGL() {
	if (cloth) {
		renderTriMesh(cloth->tm, dispMode);
	}
}

void clothInit(const char *params) {
	TriMesh *tm, *tm2;
	char fname[80], fname2[80];
	fname[0] = 0; fname2[0] = 0;
	params = extractString(params, fname, 80);
	params = extractString(params, fname2, 80);

	if (fname[0] == 0) {
		tm = charVec.tm;
	}
	else {
		tm = new TriMesh;
		if (!tm->loadFile(fname))
			return;
		tm->calcNormals();
		charVec.vis = false;
	}

	cloth = new Cloth(tm);

	if (fname2[0] != 0) {
		tm2 = new TriMesh();
		if (tm2->loadFile(fname2)) {
			int i;
			for (i=0; i < min(tm->numPts(), tm2->numPts()); i++) {
				tm->getPt(i) = tm2->getPt(i);
			}
			tm->calcNormals();
		}
	}
}

void clothStep(const char *params) {
	int i;
	int numSteps = 1;
	params = extractInt(params, &numSteps);
	
	for (i=0; i < numSteps; i++) {
		cloth->simulate();
		redrawV();
		uiWait();
	}
}

void clothDump(const char *params) {
	int i, min = 0, max = 10;
	params = extractInt(params, &min);
	params = extractInt(params, &max);
	max = min(max, cloth->numVerts);

	for (i=min; i < max; i++) {
		cout << i << ": v=" << cloth->vel->m_block[i] << ", f=" << cloth->totalForces->m_block[i] << endl;
	}
}

void clothMangle(const char *params) {
	int i;
	int mode = 0;
	params = extractInt(params, &mode);

	if (mode == 0) {
		Vec3d scale(1, 1, 1);
		params = extractDouble(params, &scale[0]);
		params = extractDouble(params, &scale[1]);
		params = extractDouble(params, &scale[2]);
		cout << "scaling by " << scale << endl;

		for (i=0; i < cloth->tm->numPts(); i++) {
			cloth->tm->getPt(i)[0] *= scale[0];
			cloth->tm->getPt(i)[1] *= scale[1];
			cloth->tm->getPt(i)[2] *= scale[2];
		}
	}
	else if (mode == 2) {
		for (i=0; i < cloth->tm->numPts(); i++) {
			cloth->tm->getPt(i)[0] += 0.5*cloth->tm->getPt(i)[1];
		}
	}
	else {
		double noise = 0.01;
		params = extractDouble(params, &noise);

		for (i=0; i < cloth->tm->numPts(); i++) {
			cloth->tm->getPt(i) += Vec3d(boundedRand(-noise, noise),boundedRand(-noise, noise),boundedRand(-noise, noise));
		}
	}
	cloth->tm->calcNormals();
	redrawV();
}

void clothMarkers(const char *params) {
	char mrefName[80], markerName[80];
	params = extractString(params, mrefName, 80);
	params = extractString(params, markerName, 80);

	// load marker refs
	ifstream is;
	int i, n;
	if (!openIFStream(&is, mrefName, "marker refs"))
		return;
	is >> n;
	markerRefs.resize(n);
	for (i=0; i < n; i++) {
		is >> markerRefs[i];
	}
	is.close();
	cout << "loaded " << n << " marker refs" << endl;

	// load markers
	markers = new MarkerSet;
	if (markerName[strlen(markerName)-3] == 'm')
		markers->loadFromMkr(markerName);
	else {
		cout << "loading from landmarks..." << endl;
		markers->loadFromLandmarks(markerName);
		// wipe out extraneous markers
		int i;
		for (i=73; i < 85; i++) {
			markers->markers[i].pos = Vec3d();
		}
	}

	// add forces
	for (i=0; i < min(n, markers->numMarkers); i++) {
		if (markerRefs[i] >= 0 && !markers->markers[i].pos.iszero()) {
			cloth->forces.push_back(new CAttract(markerRefs[i], markers->v(i), 1000, 0.01, 0));
		}
	}

//	cloth->forces.push_back(new CAttract(0, Vec3d(0,0,0), 1000, 0.01, 0));
//	cloth->forces.push_back(new CAttract(8, Vec3d(0,6,0), 1000, 0.01, 0));
}


void initClothMaster() {
	registerFunction(clothInit, "clothInit");
	registerFunction(clothStep, "clothStep");
	registerFunction(clothDump, "clothDump");
	registerFunction(clothMangle, "clothMangle");
	registerFunction(clothMarkers, "clothMarkers");
}