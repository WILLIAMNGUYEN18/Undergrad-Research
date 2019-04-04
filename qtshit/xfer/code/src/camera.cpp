#include <iostream>
#include <fstream>
#include "camera.h"
#include "ba.h"

using namespace std;

Camera::Camera() {
	className = "Camera";
}

bool Camera::load(char *fname) {
	ifstream in;
	if (!openIFStream(&in, fname, "camera"))
		return false;
	saveLoad.loadInstance(in, this);
	in.close();	
	return true;
}

void Camera::drawGL() {
	trans.glTranslate();

	double mat[16];
	rot.toMatrixD().getGLMatrix(mat);
	glMultMatrixd(mat);
}

void Camera::drawRIB(ostream &os) {
	os << "Translate " << trans[0] << " " << trans[1] << " " << trans[2] << endl;
	rot.saveRIB(os);
}

Camera Camera::interp(Camera &cam0, Camera &cam1, double w) {
	Camera c;

	c.rot = slerp(cam0.rot, cam1.rot, w);
	c.lightRot = slerp(cam0.lightRot, cam1.lightRot, w);
	c.trans = (1.0 - w) * cam0.trans + w * cam1.trans;

	return c;
}

SLInterface* Camera::newInstance(int i) {
	return new Camera[i];
}

void Camera::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();

	prop = new SLProperty();
	prop->offset = offsetof(Camera, rot);
	prop->type = SL_QUATNORM;
	classInfo->members.addT("rot", prop);

	prop = new SLProperty();
	prop->offset = offsetof(Camera, lightRot);
	prop->type = SL_QUATNORM;
	classInfo->members.addT("lightRot", prop);

	prop = new SLProperty();
	prop->offset = offsetof(Camera, trans);
	prop->type = SL_VEC3D;
	classInfo->members.addT("trans", prop);

	classInfo->size = sizeof(Camera);
	classInfo->newInstance = Camera::newInstance;
	sl->classes.addT("Camera", classInfo);
}
