/*
	doppel_config.cpp

	A class for storing the current configuration of the doppel
	tool.

	Brett Allen
	May, 2002
*/

#include "doppel_config.h"
#include "trimesh.h"

DoppelConfig config;

DoppelConfig::DoppelConfig() {
	className = "DoppelConfig";
}

void DoppelConfig::copyToGlobal() {
/*	meshUI->copyFromData();
	markerUI->copyFromData();
	skelUI->copyFromData();*/
}

void DoppelConfig::copyFromGlobal() {
}

void DoppelConfig::save(char *fname) {
	ofstream os(fname);

	if (!os.good()) {
		cout << "can't open '" << fname << "'" << endl;
		return;
	}

	saveLoad.saveInstance(os, this);
	os.close();
}

bool DoppelConfig::load(char *fname) {
	ifstream is(fname);

	if (!is.good()) {
		cout << "can't open '" << fname << "'" << endl;
		return false;
	}

	saveLoad.loadInstance(is, this);
	is.close();
	return true;
}

void DoppelConfig::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();
/*
	prop = new SLProperty();
	prop->globalPtr = &meshNT;
	prop->type = SL_NAMETABLE_T;
	meshNT.shallowSave = true;
	classInfo->members.addT("meshes", prop);

	prop = new SLProperty();
	prop->globalPtr = &markerNT;
	prop->type = SL_NAMETABLE_T;
	markerNT.shallowSave = true;
	classInfo->members.addT("markers", prop);

	prop = new SLProperty();
	prop->globalPtr = &skelNT;
	prop->type = SL_NAMETABLE_T;
	skelNT.shallowSave = true;
	classInfo->members.addT("skels", prop);

	prop = new SLProperty();
	prop->globalPtr = &mainWin->viewer->camera;
	prop->type = SL_INSTANCE;
	classInfo->members.addT("camera", prop);
*/
	classInfo->size = sizeof(DoppelConfig);
	classInfo->newInstance = DoppelConfig::newInstance;
	sl->classes.addT("DoppelConfig", classInfo);
}

SLInterface *DoppelConfig::newInstance(int i) {
	return new DoppelConfig[i];
}