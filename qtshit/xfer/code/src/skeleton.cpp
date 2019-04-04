#include <fstream>
#include "skeleton.h"
#include "ba.h"

// Skeleton ==================================================

Skeleton::Skeleton() {
	className = "Skeleton";
}

void Skeleton::init() {
	int i, j;
	SkelTransform *curT;

	// set names and indices
	for (i=0; i < transforms.size(); i++) {
		curT = transforms.getT(i);

		strcpy(curT->name, transforms.getName(i));
		curT->index = i;
	}

	// set children
	for (i=0; i < transforms.size(); i++) {
		curT = transforms.getT(i);

		j = transforms.lookupName(curT->parent);
		if (j < 0)
			curT->parentPtr = NULL;
		else {
			transforms.getT(j)->addChild(curT);
		}
	}

	// set symmetric references
	for (i=0; i < transforms.size(); i++) {
		transforms.getT(i)->initRefs(this);
	}

	// make intrinsic map
	numIntrinsicDofs = 0;
	numPoseDofs = 0;
	for (i=0; i < transforms.size(); i++) {
		curT = transforms.getT(i);
		if (curT->isIntrinsic)
			numIntrinsicDofs += curT->numDofs();
		else
			numPoseDofs += curT->numDofs();
	}
	numDofs = numIntrinsicDofs + numPoseDofs;
	int dof = 0;
	dofIntrins = new bool[numDofs];
	for (i=0; i < transforms.size(); i++) {
		curT = transforms.getT(i);
		if (curT->isIntrinsic)
			for (j=0; j < curT->numDofs(); j++)
				dofIntrins[dof++] = true;
		else
			for (j=0; j < curT->numDofs(); j++)
				dofIntrins[dof++] = false;
	}

	updateCoords();
}

void Skeleton::drawGL(SkelTransform *tr, double alpha) {
	int i;

	if (tr == NULL)
		tr = transforms.getT(0);

	glPushMatrix();
		tr->drawGL(alpha);
		for (i=0; i < tr->children.size(); i++) {
			drawGL(tr->children[i], alpha);
		}
	glPopMatrix();
}

void Skeleton::copyVals(Skeleton *sk, int mode) {
	int i;
	
	for (i=0; i < sk->transforms.size(); i++) {
		SkelTransform *st = transforms.getT(sk->transforms.getName(i));

		if (st) {
			if (st->isIntrinsic && (mode & COPY_INT))
				st->copyVal(sk->transforms.getT(i));
			else if (!st->isIntrinsic && (mode & COPY_POSE))
				st->copyVal(sk->transforms.getT(i));
		}
	}
}

void Skeleton::interpVals(Skeleton *sk0, Skeleton *sk1, double interp, bool poseOnly) {
	int i;
	
	for (i=0; i < transforms.size(); i++) {
		SkelTransform *st0 = sk0->transforms.getT(transforms.getName(i));
		SkelTransform *st1 = sk1->transforms.getT(transforms.getName(i));

		if (st0 && st1) {
			if ((!poseOnly) || (!st0->isIntrinsic))
				transforms.getT(i)->interpVal(st0, st1, interp);
		}
	}
}

void Skeleton::updateCoords() {
	int i;
	SkelTransform *curT;

	for (i=0; i < transforms.size(); i++) {
		curT = transforms.getT(i);

		curT->updateCoord();
		if (curT->parentPtr)
			curT->globalCoord = curT->parentPtr->globalCoord * curT->curCoord;
		else
			curT->globalCoord = curT->curCoord;
	}
}

void Skeleton::updateDerivs() {
	int i;
	for (i=0; i < transforms.size(); i++) {
		transforms.getT(i)->updateDerivs();
	}
}

void Skeleton::updateGlobalDerivs() {
	int dof = 0;
	int i, j;

	for (i=0; i < transforms.size(); i++) {
		SkelTransform *curTrans = transforms.getT(i);
		curTrans->dofInd = dof;
		curTrans->updateGlobalDerivs(numDofs);
		dof += curTrans->numDofs();
	}
}

void Skeleton::zero() {
	int i;
	for (i=0; i < transforms.size(); i++) {
		transforms.getT(i)->zero();
	}
}

Skeleton *Skeleton::load(char *fname) {
	Skeleton *sk;

	ifstream is(fname);
	if (!is.good()) {
		cout << "can't open: '" << fname << "'" << endl;
		return NULL;
	}
	cout << "loading skel: " << fname << endl;

	sk = (Skeleton*)saveLoad.load(is);
	sk->init();

	return sk;
}

void Skeleton::loadPose(istream &in) {
	char partName[256];

	while (1) {
		in >> partName;
		if (!in.good() || strlen(partName) < 1)
			break;

		SkelTransform *t = transforms.getT(partName);

		if (t == NULL)
			cout << "unknown part: '" << partName << "'" << endl;
		else {
			t->loadPose(in);
		}
	}

	updateCoords();
}

void Skeleton::savePose(ostream &out) {
	int i;
	for (i=0; i < transforms.size(); i++) {
		SkelTransform *t = transforms.getT(i);

		if (t->numDofs() == 0)
			continue;

		t->savePose(out);
	}
}

void Skeleton::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();

	prop = new SLProperty();
	prop->offset = offsetof(Skeleton, transforms);
	prop->type = SL_NAMETABLE_T;
	classInfo->members.addT("transforms", prop);

	classInfo->size = sizeof(Skeleton);
	classInfo->newInstance = Skeleton::newInstance;
	sl->classes.addT("Skeleton", classInfo);
}

SLInterface *Skeleton::newInstance(int count) {
	return new Skeleton[count];
}

void Skeleton::copyFrom(Skeleton *skel) {
	int i;
	for (i=0; i < skel->transforms.size(); i++)
		transforms.addT(skel->transforms.getName(i), skel->transforms.getT(i)->clone());
	init();
}

void Skeleton::mirrorLR() {
	int i;
	char curName[256];
	SkelTransform *curLeft, *curRight;

	for (i=0; i < transforms.size(); i++) {
		strcpy(curName, transforms.getName(i));
		if ((curName[0] == 'l') && ((curName[strlen(curName)-1] == 'Q') || (curName[strlen(curName)-1] == 'A'))) {
			curLeft = transforms.getT(i);
			curName[0] = 'r';
			curRight = transforms.getT(curName);
			if (curLeft && curRight) {
				curRight->mirrorVal(curLeft);
			}
			else
				cout << "warning: can't find node: " << curName << endl;
		}
	}
}

void Skeleton::allocDerivs(int numDofs) {
	int i;
	for (i=0; i < transforms.size(); i++) {
		if (transforms.getT(i)->globalDerivs) {
			delete [](transforms.getT(i)->globalDerivs);
		}
		transforms.getT(i)->globalDerivs = new Mat4d[numDofs];
	}
}

// SkelTransform =============================================

SkelTransform::SkelTransform() {
	className = "SkelTransform";

	isIntrinsic = false;

	index = 0;
	name[0] = 0;
	color = Vec3d(0, 0, 0);

	parent[0] = 0;

	globalDerivs = NULL;
}

SkelTransform::~SkelTransform() {
	if (globalDerivs)
		delete []globalDerivs;
}

void SkelTransform::addChild(SkelTransform *st) {
	children.push_back(st);
	st->parentPtr = this;
	strcpy(st->parent, name);
}

int SkelTransform::numDofs() {
	return 0;
}

void SkelTransform::loadDofs(double *v, double scale) {
}

void SkelTransform::unloadDofs(double *v, double scale) {
}

double &SkelTransform::getDofAddr(int dofI) {
	static double foo;
	return foo;
}

void SkelTransform::updateGlobalDerivs(int maxDof) {
	int dof;

	if (!parentPtr) {
		// zero all derivatives
		memset(globalDerivs, 0, sizeof(Mat4d) * maxDof);
	}
	else {
		// add current transformation to previous derivative
		for (dof = 0; dof < maxDof; dof++) {
			globalDerivs[dof] = parentPtr->globalDerivs[dof] * curCoord.mat;
		}
	}

	// now handle my own dofs
	for (dof = 0; dof < numDofs(); dof++) {
		if (parentPtr)
			globalDerivs[dof + dofInd] += parentPtr->globalCoord.mat * (*curCoord.deriv[dof]);
		else
			globalDerivs[dof + dofInd] += (*curCoord.deriv[dof]);
	}
}

void SkelTransform::loadPose(istream &in) {
	int i;
	for (i=0; i < numDofs(); i++) {
		in >> getDofAddr(i);
	}
}

void SkelTransform::savePose(ostream &os) {
	if (numDofs() == 0)
		return;

	os << name;
	int i;
	for (i=0; i < numDofs(); i++) {
		os << " " << getDofAddr(i);
	}
	os << endl;
}

void SkelTransform::renderSkelRIB(bool showMarkers, ostream &rib) {
	if (showMarkers) {
/*		int i;
		for (i=0; i < model->numMarkers; i++) {
			if (model->markers[i].transform == index) {
				if (!fancy) {
					(color * 0.5).glColor();
					glDisable(GL_LIGHTING);
					glLineWidth(1);
					glBegin(GL_LINES);
						glVertex3d(0, 0, 0);
						model->markers[i].curVal.glVertex();
					glEnd();
					glEnable(GL_LIGHTING);
				}

				color.glColor();
				glbSphere(model->markers[i].curVal, 0.005);
			}
		}*/
	}
}

void SkelTransform::registerProps(SLClassInfo *classInfo) {
	SLProperty *prop;

	prop = new SLProperty();
	prop->offset = offsetof(SkelTransform, parent);
	prop->type = SL_CHAR256;
	classInfo->members.addT("parent", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelTransform, color);
	prop->type = SL_VEC3D;
	classInfo->members.addT("color", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelTransform, isIntrinsic);
	prop->type = SL_BOOL;
	classInfo->members.addT("isIntrinsic", prop);
}

void SkelTransform::copyFrom(SkelTransform *st) {
	index = st->index;
	strcpy(name, st->name);
	color = st->color;
	curCoord = st->curCoord;
	globalCoord = st->globalCoord;
	isIntrinsic = st->isIntrinsic;
	strcpy(parent, st->parent);
}

// SkelTranslation ===========================================

SkelTranslation::SkelTranslation() {
	className = "SkelTranslation";
	xyz = 4 + 2 + 1;
	box = Vec3d();
}

int SkelTranslation::numDofs() {
	int nd = 0;
	if (xyz & 1) nd++;
	if (xyz & 2) nd++;
	if (xyz & 4) nd++;
	return nd;
}

void SkelTranslation::loadDofs(double *v, double scale) {
	int ind = 0;
	if (xyz & 1)
		curVal[0] = v[ind++] * scale;
	if (xyz & 2)
		curVal[1] = v[ind++] * scale;
	if (xyz & 4)
		curVal[2] = v[ind++] * scale;
}

void SkelTranslation::unloadDofs(double *v, double scale) {
	int ind = 0;
	if (xyz & 1)
		v[ind++] = curVal[0] / scale;
	if (xyz & 2)
		v[ind++] = curVal[1] / scale;
	if (xyz & 4)
		v[ind++] = curVal[2] / scale;
}

double &SkelTranslation::getDofAddr(int dofI) {
	int index = dofI;
	if ((xyz & 1) == 0)
		index++;
	if (((xyz & 2) == 0) && (index >= 1)) {
		index++;
	}
	return curVal[index];
}

void SkelTranslation::updateCoord() {
	curCoord.mat = Mat4d();
	curCoord.mat[0][3] = curVal[0];
	curCoord.mat[1][3] = curVal[1];
	curCoord.mat[2][3] = curVal[2];

	curCoord.q = QuatNorm();
	curCoord.v = curVal;
}

void SkelTranslation::updateDerivs() {
	Mat4d m;
	int i;

	curCoord.initDeriv(numDofs());

	int dof = 0;
	if (xyz & 1) {
		m *= 0;
		m[0][3] = 1;
		curCoord.setDeriv(dof, m);
		dof++;
	}
	if (xyz & 2) {
		m *= 0;
		m[1][3] = 1;
		curCoord.setDeriv(dof, m);
		dof++;
	}
	if (xyz & 4) {
		m *= 0;
		m[2][3] = 1;
		curCoord.setDeriv(dof, m);
		dof++;
	}
}

void SkelTranslation::copyVal(SkelTransform *k) {
	curVal = ((SkelTranslation*)k)->curVal;
}

void SkelTranslation::interpVal(SkelTransform *k0, SkelTransform *k1, double interp) {
	curVal = (1.0 - interp) * ((SkelTranslation*)k0)->curVal +
			interp * ((SkelTranslation*)k1)->curVal;
}

void SkelTranslation::zero() {
	curVal = Vec3d();
}

void SkelTranslation::loadPose(istream &in) {
	if (xyz & 1)
		in >> curVal[0];
	if (xyz & 2)
		in >> curVal[1];
	if (xyz & 4)
		in >> curVal[2];
}

double SkelTranslation::distance(SkelTransform *k) {
	return (((SkelTranslation*)k)->curVal - curVal).length2();
}

void SkelTranslation::drawGL(double alpha) {
	if (!color.iszero()) {
		color.glColor(alpha);
		if (!box.iszero()) {
			glbBox(curVal * 0.5, box);
		}
		else {
			glbDirectedCyl(curVal, curVal.length(), 0.005, 0.005);
		}
	}
	curVal.glTranslate();
}

void SkelTranslation::renderSkelRIB(bool showMarkers, ostream &rib)  {
	if (!color.iszero()) {
		rib << "Color [ " << color << " ]" << endl;
//		glbDirectedCyl(curVal, curVal.length(), 0.005, 0.005);
	}
	rib << "Translate " << curVal << endl;
	SkelTransform::renderSkelRIB(showMarkers, rib);
}

void SkelTranslation::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();

	SkelTransform::registerProps(classInfo);

	prop = new SLProperty();
	prop->offset = offsetof(SkelTranslation, curVal);
	prop->type = SL_VEC3D;
	classInfo->members.addT("curVal", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelTranslation, box);
	prop->type = SL_VEC3D;
	classInfo->members.addT("box", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelTranslation, xyz);
	prop->type = SL_INT;
	classInfo->members.addT("xyz", prop);

	classInfo->size = sizeof(SkelTranslation);
	classInfo->newInstance = SkelTranslation::newInstance;
	sl->classes.addT("SkelTranslation", classInfo);
}

SLInterface *SkelTranslation::newInstance(int count) {
	return new SkelTranslation[count];
}

void SkelTranslation::copyFrom(SkelTransform *st) {
	SkelTransform::copyFrom(st);
	curVal = ((SkelTranslation*)st)->curVal;
	xyz = ((SkelTranslation*)st)->xyz;
	box = ((SkelTranslation*)st)->box;
}

SkelTransform *SkelTranslation::clone() {
	SkelTransform *ret = new SkelTranslation();
	ret->copyFrom(this);
	return ret;
}

// SkelQuatRotation ==========================================

SkelQuatRotation::SkelQuatRotation() {
	className = "SkelQuatRotation";
}

int SkelQuatRotation::numDofs() {
	return 4;
}

void SkelQuatRotation::loadDofs(double *v, double scale) {
	curQuat.x = v[0] * scale;
	curQuat.y = v[1] * scale;
	curQuat.z = v[2] * scale;
	curQuat.w = v[3] * scale;
}

void SkelQuatRotation::unloadDofs(double *v, double scale) {
	v[0] = curQuat.x / scale;
	v[1] = curQuat.y / scale;
	v[2] = curQuat.z / scale;
	v[3] = curQuat.w / scale;
}

double &SkelQuatRotation::getDofAddr(int dofI) {
	return curQuat[dofI];
}

void SkelQuatRotation::updateCoord() {
	curCoord.mat = curQuat.toMatrixD();

	curCoord.q = curQuat;
	curCoord.q.normalize();
	curCoord.v = Vec3d();
}

void SkelQuatRotation::updateDerivs() {
	curCoord.initDeriv(numDofs());

	Mat4d val;
	Mat4d curDerivs[4];
	curQuat.getMatrices(val, curDerivs);

	int i;
	for (i=0; i < 4; i++)
		curCoord.setDeriv(i, curDerivs[i]); //  * scaleDeriv()
}

void SkelQuatRotation::copyVal(SkelTransform *k) {
	curQuat = ((SkelQuatRotation*)k)->curQuat;
}

void SkelQuatRotation::interpVal(SkelTransform *k0, SkelTransform *k1, double interp) {
	curQuat = slerp(((SkelQuatRotation*)k0)->curQuat, ((SkelQuatRotation*)k1)->curQuat, interp);
}

void SkelQuatRotation::mirrorVal(SkelTransform *k) {
	curQuat = ((SkelQuatRotation*)k)->curQuat;
	curQuat.x *= -1;
	curQuat.w *= -1;
}

void SkelQuatRotation::zero() {
	curQuat = QuatNorm();
}

void SkelQuatRotation::loadPose(istream &in) {
	in >> curQuat.x >> curQuat.y >> curQuat.z >> curQuat.w;
}

void SkelQuatRotation::normalize() {
	curQuat.normalize();
}

double SkelQuatRotation::distance(SkelTransform *k) {
	QuatNorm q0 = curQuat;
	QuatNorm q1 = ((SkelQuatRotation*)k)->curQuat;
	q0.normalize();
	q1.normalize();
//	return sqr(q0.x - q1.x) + sqr(q0.y - q1.y) + sqr(q0.z - q1.z) + sqr(q0.w - q1.w);
	return sqr(quatDist(q0, q1));
}

void SkelQuatRotation::drawGL(double alpha) {
	if (!color.iszero()) {
		color.glColor(alpha);
		glbSphere(Vec3d(0, 0, 0), 0.02);
	}

	double mat[16];
	curQuat.toMatrixD().getGLMatrix(mat);
	glMultMatrixd(mat);
}

void SkelQuatRotation::renderSkelRIB(bool showMarkers, ostream &rib)  {
	rib << "Color 0.5 0.5 0.5" << endl;
	rib << "Sphere 0.02 -0.02 0.02 360" << endl;

	curQuat.saveRIB(rib);

	SkelTransform::renderSkelRIB(showMarkers, rib);
}

void SkelQuatRotation::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();

	SkelTransform::registerProps(classInfo);

	prop = new SLProperty();
	prop->offset = offsetof(SkelQuatRotation, curQuat);
	prop->type = SL_QUATNORM;
	classInfo->members.addT("curQuat", prop);

	classInfo->size = sizeof(SkelQuatRotation);
	classInfo->newInstance = SkelQuatRotation::newInstance;
	sl->classes.addT("SkelQuatRotation", classInfo);
}

SLInterface *SkelQuatRotation::newInstance(int count) {
	return new SkelQuatRotation[count];
}

void SkelQuatRotation::copyFrom(SkelTransform *st) {
	SkelTransform::copyFrom(st);
	curQuat = ((SkelQuatRotation*)st)->curQuat;
}

SkelTransform *SkelQuatRotation::clone() {
	SkelTransform *ret = new SkelQuatRotation();
	ret->copyFrom(this);
	return ret;
}

// SkelEulerRotation =========================================

SkelEulerRotation::SkelEulerRotation() {
	className = "SkelEulerRotation";

	axis = 0;
	curAngle = 0;
	
	min = -PI;
	max = PI*3.0/2.0;
}

int SkelEulerRotation::numDofs() {
	return 1;
}

void SkelEulerRotation::loadDofs(double *v, double scale) {
	curAngle = v[0] * scale;
}

void SkelEulerRotation::unloadDofs(double *v, double scale) {
	v[0] = curAngle / scale;
}

double &SkelEulerRotation::getDofAddr(int) {
	return curAngle;
}

void SkelEulerRotation::updateCoord() {
	switch (axis) {
	case 0:
		curQuat = QuatNorm(sin(curAngle / 2.0), 0, 0, -cos(curAngle / 2.0));
		break;
	case 1:
		curQuat = QuatNorm(0, sin(curAngle / 2.0), 0, -cos(curAngle / 2.0));
		break;
	case 2:
		curQuat = QuatNorm(0, 0, sin(curAngle / 2.0), -cos(curAngle / 2.0));
		break;
	}

	curCoord.mat = curQuat.toMatrixD();
	curCoord.q = curQuat;
	curCoord.v = Vec3d();
}

void SkelEulerRotation::updateDerivs() {
	Mat4d derivMat;
	derivMat *= 0;

	curCoord.initDeriv(numDofs());

	double cosTheta = cos(curAngle);
	double sinTheta = sin(curAngle);

	switch( axis ) {
	case 0:
		derivMat[1][1] = -sinTheta;
		derivMat[1][2] = -cosTheta;
		derivMat[2][1] = cosTheta;
		derivMat[2][2] = -sinTheta;
		break;
	case 1:
		derivMat[0][0] = -sinTheta;
		derivMat[0][2] = cosTheta;
		derivMat[2][0] = -cosTheta;
		derivMat[2][2] = -sinTheta;
		break;
	case 2:
		derivMat[0][0] = -sinTheta;
		derivMat[0][1] = -cosTheta;
		derivMat[1][0] = cosTheta;
		derivMat[1][1] = -sinTheta;
		break;
	}

	curCoord.setDeriv(0, derivMat); // * scaleDeriv());
}

void SkelEulerRotation::copyVal(SkelTransform *k) {
	curAngle = ((SkelEulerRotation*)k)->curAngle;
}

void SkelEulerRotation::interpVal(SkelTransform *k0, SkelTransform *k1, double interp) {
	curAngle = (1.0 - interp) * ((SkelEulerRotation*)k0)->curAngle +
		interp * ((SkelEulerRotation*)k1)->curAngle;
}

void SkelEulerRotation::mirrorVal(SkelTransform *k) {
	curAngle = -((SkelEulerRotation*)k)->curAngle;
}

void SkelEulerRotation::zero() {
	curAngle = 0;
}

void SkelEulerRotation::loadPose(istream &in) {
	in >> curAngle;
}

void SkelEulerRotation::normalize() {
	while (curAngle < min)
		curAngle += 2.0*PI;
	while (curAngle > max)
		curAngle -= 2.0*PI;
}

double SkelEulerRotation::distance(SkelTransform *k) {
	return sqr(curAngle - ((SkelEulerRotation*)k)->curAngle);
}

void SkelEulerRotation::drawGL(double alpha) {
	if (!color.iszero()) {
		color.glColor(alpha);
		switch (axis) {
		case 0:
			glbDirectedCyl(Vec3d(1, 0, 0), 0.04, 0.02, 0);
			break;
		case 1:
			glbDirectedCyl(Vec3d(0, 1, 0), 0.04, 0.02, 0);
			break;
		case 2:
			glbDirectedCyl(Vec3d(0, 0, 1), 0.04, 0.02, 0);
			break;
		}
	}

	switch (axis) {
	case 0:
		glRotated(RAD_TO_DEG * curAngle, 1, 0, 0);
		break;
	case 1:
		glRotated(RAD_TO_DEG * curAngle, 0, 1, 0);
		break;
	case 2:
		glRotated(RAD_TO_DEG * curAngle, 0, 0, 1);
		break;
	}
}

void SkelEulerRotation::renderSkelRIB(bool showMarkers, ostream &rib)  {
	rib << "Color 0.5 0.5 0.5" << endl;
	/*
	switch (axis) {
	case 0:
		glbDirectedCyl(Vec3d(1, 0, 0), 0.04, 0.02, 0);
		break;
	case 1:
		glbDirectedCyl(Vec3d(0, 1, 0), 0.04, 0.02, 0);
		break;
	case 2:
		glbDirectedCyl(Vec3d(0, 0, 1), 0.04, 0.02, 0);
		break;
	}
*/
	switch (axis) {
	case 0:
		rib << "Rotate " << (RAD_TO_DEG * curAngle) << " 1 0 0" << endl;
		break;
	case 1:
		rib << "Rotate " << (RAD_TO_DEG * curAngle) << " 0 1 0" << endl;
		break;
	case 2:
		rib << "Rotate " << (RAD_TO_DEG * curAngle) << " 0 0 1" << endl;
		break;
	}

	SkelTransform::renderSkelRIB(showMarkers, rib);
}

void SkelEulerRotation::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();

	SkelTransform::registerProps(classInfo);

	prop = new SLProperty();
	prop->offset = offsetof(SkelEulerRotation, axis);
	prop->type = SL_INT;
	classInfo->members.addT("axis", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelEulerRotation, curAngle);
	prop->type = SL_DOUBLE;
	classInfo->members.addT("curAngle", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelEulerRotation, min);
	prop->type = SL_DOUBLE;
	classInfo->members.addT("min", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelEulerRotation, max);
	prop->type = SL_DOUBLE;
	classInfo->members.addT("max", prop);

	classInfo->size = sizeof(SkelEulerRotation);
	classInfo->newInstance = SkelEulerRotation::newInstance;
	sl->classes.addT("SkelEulerRotation", classInfo);
}

SLInterface *SkelEulerRotation::newInstance(int count) {
	return new SkelEulerRotation[count];
}

void SkelEulerRotation::copyFrom(SkelTransform *st) {
	SkelTransform::copyFrom(st);
	axis = ((SkelEulerRotation*)st)->axis;
	curAngle = ((SkelEulerRotation*)st)->curAngle;
	min = ((SkelEulerRotation*)st)->min;
	max = ((SkelEulerRotation*)st)->max;
	curQuat = ((SkelEulerRotation*)st)->curQuat;
}

SkelTransform *SkelEulerRotation::clone() {
	SkelTransform *ret = new SkelEulerRotation();
	ret->copyFrom(this);
	return ret;
}

// SkelPolarAxisRotation =========================================

SkelPolarAxisRotation::SkelPolarAxisRotation() {
	className = "SkelPolarAxisRotation";

	curPhi = 0;
	curPsi = 0;
	curTheta = 0;

	min = -PI;
	max = PI*3.0/2.0;
}

int SkelPolarAxisRotation::numDofs() {
	return 3;
}

void SkelPolarAxisRotation::loadDofs(double *v, double scale) {
	curPhi = v[0] * scale;
	curPsi = v[1] * scale;
	curTheta = v[2] * scale;
}

void SkelPolarAxisRotation::unloadDofs(double *v, double scale) {
	v[0] = curPhi / scale;
	v[1] = curPsi / scale;
	v[2] = curTheta / scale;
}

double &SkelPolarAxisRotation::getDofAddr(int dofI) {
	switch (dofI) {
	case 0:
		return curPhi;
	case 1:
		return curPsi;
	default:
		return curTheta;
	}
}

void SkelPolarAxisRotation::updateCoord() {
	curAxis = Vec3d(sin(curPhi) * cos(curPsi), sin(curPhi) * sin(curPsi), cos(curPhi));
	double sa = sin(curTheta / 2.0);
	curQuat = QuatNorm(sa * curAxis[0], sa * curAxis[1], sa * curAxis[2], -cos(curTheta / 2.0));

	curCoord.mat = curQuat.toMatrixD();
	curCoord.q = curQuat;
	curCoord.v = Vec3d();
}

void SkelPolarAxisRotation::updateDerivs() {
	Mat4d derivMat;
	derivMat *= 0;

	curCoord.initDeriv(numDofs());

	// Maple-generated code
	double t1 = cos(curTheta);
	double t2 = 1.0-t1;
	double t3 = sin(curPhi);
	double t4 = t3*t3;
	double t5 = t2*t4;
	double t6 = cos(curPsi);
	double t7 = t6*t6;
	double t8 = t5*t7;
	double t10 = sin(curPsi);
	double t11 = t6*t10;
	double t12 = t5*t11;
	double t13 = cos(curPhi);
	double t14 = sin(curTheta);
	double t15 = t13*t14;
	double t17 = t2*t3;
	double t18 = t6*t13;
	double t19 = t17*t18;
	double t20 = t3*t10;
	double t21 = t20*t14;
	double t22 = t19+t21;
	double t24 = t10*t10;
	double t25 = t5*t24;
	double t27 = t10*t13;
	double t28 = t17*t27;
	double t29 = t3*t6;
	double t30 = t29*t14;
	double t31 = t28-t30;
	double t32 = t19-t21;
	double t33 = t28+t30;
	double t34 = t13*t13;
	double t35 = t2*t34;
	double t37 = t14*t4;
	double t40 = t37*t11;
	double t41 = t13*t1;
	double t43 = t3*t14;
	double t44 = t43*t18;
	double t45 = t20*t1;
	double t50 = t43*t27;
	double t51 = t29*t1;
	double t57 = 2.0*t12;
	double t58 = -t25+t8;
	double t64 = 2.0*t17*t11*t13;
	double t66 = t35*t6;
	double t67 = t5*t6;
	double t68 = t27*t14;
	double t74 = t35*t10;
	double t75 = t5*t10;
	double t76 = t18*t14;
	// hopefully this matrix is the same as the quaternion-generated one...
	Mat4d M;
	M[0][0] = t1+t8;
	M[0][1] = t12-t15;
	M[0][2] = t22;
	M[0][3] = 0.0;
	M[1][0] = t12+t15;
	M[1][1] = t1+t25;
	M[1][2] = t31;
	M[1][3] = 0.0;
	M[2][0] = t32;
	M[2][1] = t33;
	M[2][2] = t1+t35;
	M[2][3] = 0.0;
	M[3][0] = 0.0;
	M[3][1] = 0.0;
	M[3][2] = 0.0;
	M[3][3] = 1.0;
	
	Mat4d curDeriv;
	// theta
	curDeriv[0][0] = -t14+t37*t7;
	curDeriv[0][1] = t40-t41;
	curDeriv[0][2] = t44+t45;
	curDeriv[0][3] = 0.0;
	curDeriv[1][0] = t40+t41;
	curDeriv[1][1] = -t14+t37*t24;
	curDeriv[1][2] = t50-t51;
	curDeriv[1][3] = 0.0;
	curDeriv[2][0] = t44-t45;
	curDeriv[2][1] = t50+t51;
	curDeriv[2][2] = -t14+t14*t34;
	curDeriv[2][3] = 0.0;
	curDeriv[3][0] = 0.0;
	curDeriv[3][1] = 0.0;
	curDeriv[3][2] = 0.0;
	curDeriv[3][3] = 0.0;
	curCoord.setDeriv(2, curDeriv); // * scaleDeriv());
	// psi
	curDeriv[0][0] = -t57;
	curDeriv[0][1] = t58;
	curDeriv[0][2] = -t31;
	curDeriv[0][3] = 0.0;
	curDeriv[1][0] = t58;
	curDeriv[1][1] = t57;
	curDeriv[1][2] = t22;
	curDeriv[1][3] = 0.0;
	curDeriv[2][0] = -t33;
	curDeriv[2][1] = t32;
	curDeriv[2][2] = 0.0;
	curDeriv[2][3] = 0.0;
	curDeriv[3][0] = 0.0;
	curDeriv[3][1] = 0.0;
	curDeriv[3][2] = 0.0;
	curDeriv[3][3] = 0.0;
	curCoord.setDeriv(1, curDeriv); // * scaleDeriv());
	// phi
	curDeriv[0][0] = 2.0*t17*t7*t13;
	curDeriv[0][1] = t64+t43;
	curDeriv[0][2] = t66-t67+t68;
	curDeriv[0][3] = 0.0;
	curDeriv[1][0] = t64-t43;
	curDeriv[1][1] = 2.0*t17*t24*t13;
	curDeriv[1][2] = t74-t75-t76;
	curDeriv[1][3] = 0.0;
	curDeriv[2][0] = t66-t67-t68;
	curDeriv[2][1] = t74-t75+t76;
	curDeriv[2][2] = -2.0*t2*t13*t3;
	curDeriv[2][3] = 0.0;
	curDeriv[3][0] = 0.0;
	curDeriv[3][1] = 0.0;
	curDeriv[3][2] = 0.0;
	curDeriv[3][3] = 0.0;
	curCoord.setDeriv(0, curDeriv); // * scaleDeriv());
}

void SkelPolarAxisRotation::copyVal(SkelTransform *k) {
	curPhi = ((SkelPolarAxisRotation*)k)->curPhi;
	curPsi = ((SkelPolarAxisRotation*)k)->curPsi;
	curTheta = ((SkelPolarAxisRotation*)k)->curTheta;
}

void SkelPolarAxisRotation::interpVal(SkelTransform *k0, SkelTransform *k1, double interp) {
	curPhi = (1.0 - interp) * ((SkelPolarAxisRotation*)k0)->curPhi + interp * ((SkelPolarAxisRotation*)k1)->curPhi;
	curPsi = (1.0 - interp) * ((SkelPolarAxisRotation*)k0)->curPsi + interp * ((SkelPolarAxisRotation*)k1)->curPsi;
	curTheta = (1.0 - interp) * ((SkelPolarAxisRotation*)k0)->curTheta + interp * ((SkelPolarAxisRotation*)k1)->curTheta;
}

void SkelPolarAxisRotation::mirrorVal(SkelTransform *k) {
	curPhi = ((SkelPolarAxisRotation*)k)->curPhi;
	curPsi = PI - ((SkelPolarAxisRotation*)k)->curPsi;
	curTheta = -((SkelPolarAxisRotation*)k)->curTheta;
}

void SkelPolarAxisRotation::zero() {
	curTheta = 0;
}

void SkelPolarAxisRotation::loadPose(istream &in) {
	in >> curPhi >> curPsi >> curTheta;
}

void SkelPolarAxisRotation::normalize() {
	while (curTheta < min)
		curTheta += 2.0*PI;
	while (curTheta > max)
		curTheta -= 2.0*PI;
}

double SkelPolarAxisRotation::distance(SkelTransform *k) {
	return sqr(curTheta - ((SkelPolarAxisRotation*)k)->curTheta);
}

void SkelPolarAxisRotation::drawGL(double alpha)  {
	if (!color.iszero()) {
		color.glColor(alpha);
		glbDirectedCyl(curAxis, 0.06, 0.02, 0);
	}

	glRotated(RAD_TO_DEG * curTheta, curAxis[0], curAxis[1], curAxis[2]);
}

void SkelPolarAxisRotation::renderSkelRIB(bool showMarkers, ostream &rib)  {
	rib << "Color 0.5 0.5 0.5" << endl;
//		glbDirectedCyl(curAxis, 0.06, 0.02, 0);

	rib << "Rotate " << (RAD_TO_DEG * curTheta) << " " << curAxis << endl;

	SkelTransform::renderSkelRIB(showMarkers, rib);
}

void SkelPolarAxisRotation::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();

	SkelTransform::registerProps(classInfo);

	prop = new SLProperty();
	prop->offset = offsetof(SkelPolarAxisRotation, curPhi);
	prop->type = SL_DOUBLE;
	classInfo->members.addT("curPhi", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelPolarAxisRotation, curPsi);
	prop->type = SL_DOUBLE;
	classInfo->members.addT("curPsi", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelPolarAxisRotation, curTheta);
	prop->type = SL_DOUBLE;
	classInfo->members.addT("curTheta", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelPolarAxisRotation, min);
	prop->type = SL_DOUBLE;
	classInfo->members.addT("min", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelPolarAxisRotation, max);
	prop->type = SL_DOUBLE;
	classInfo->members.addT("max", prop);

	classInfo->size = sizeof(SkelPolarAxisRotation);
	classInfo->newInstance = SkelPolarAxisRotation::newInstance;
	sl->classes.addT("SkelPolarAxisRotation", classInfo);
}

SLInterface *SkelPolarAxisRotation::newInstance(int count) {
	return new SkelPolarAxisRotation[count];
}

void SkelPolarAxisRotation::copyFrom(SkelTransform *st) {
	SkelTransform::copyFrom(st);
	curPhi = ((SkelPolarAxisRotation*)st)->curPhi;
	curPsi = ((SkelPolarAxisRotation*)st)->curPsi;
	curTheta = ((SkelPolarAxisRotation*)st)->curTheta;
}

SkelTransform *SkelPolarAxisRotation::clone() {
	SkelTransform *ret = new SkelPolarAxisRotation();
	ret->copyFrom(this);
	return ret;
}


// SkelSymmetricTranslation ==================================

SkelSymmetricTranslation::SkelSymmetricTranslation() {
	className = "SkelSymmetricTranslation";

	origTranslation[0] = 0;
	orig = NULL;
	axisMult[0] = 1;
	axisMult[1] = 1;
	axisMult[2] = 1;

	box = Vec3d();
}

void SkelSymmetricTranslation::initRefs(Skeleton *skel) {
	int origInd = skel->transforms.lookupName(origTranslation);
	if (origInd < 0) {
		cout << "warning: " << name << "'s original translation, " << 
			origTranslation << ", is unknown." << endl;
		orig = NULL;
	}
	else {
		orig = (SkelTranslation*)skel->transforms.getT(origInd);
	}
}

void SkelSymmetricTranslation::updateCoord() {
	curVal = orig->curVal;

	curVal[0] *= axisMult[0];
	curVal[1] *= axisMult[1];
	curVal[2] *= axisMult[2];

	curCoord.mat = Mat4d();
	curCoord.mat[0][3] = curVal[0];
	curCoord.mat[1][3] = curVal[1];
	curCoord.mat[2][3] = curVal[2];

	curCoord.q = QuatNorm();
	curCoord.v = curVal;
}

void SkelSymmetricTranslation::updateGlobalDerivs(int maxDof) {
	SkelTransform::updateGlobalDerivs(maxDof);

	if (orig) {
		int dof;
		for (dof = 0; dof < orig->numDofs(); dof++) {
			Mat4d deriv = *orig->curCoord.deriv[dof];
			int row, col;
			for (row=0; row < 3; row++)
				for (col=0; col < 4; col++)
					deriv[row][col] *= axisMult[row];

			if (parentPtr)
				globalDerivs[dof + orig->dofInd] += parentPtr->globalCoord.mat * deriv;
			else
				globalDerivs[dof + orig->dofInd] += deriv;
		}
	}
}

void SkelSymmetricTranslation::loadPose(istream &in) {
	Vec3d junk;
	in >> junk[0] >> junk[1] >> junk[2];
}

double SkelSymmetricTranslation::distance(SkelTransform *k) {
	return (curVal - ((SkelSymmetricTranslation*)k)->curVal).length2();
}

void SkelSymmetricTranslation::drawGL(double alpha)  {
	if (!color.iszero()) {
		color.glColor(alpha);
		if (!box.iszero()) {
			glbBox(curVal * 0.5, box);
		}
		else {
			glbDirectedCyl(curVal, curVal.length(), 0.005, 0.005);
		}
	}
	curVal.glTranslate();
}

void SkelSymmetricTranslation::renderSkelRIB(bool showMarkers, ostream &rib)  {
	if (!color.iszero()) {
		rib << "Color " << color << endl;
//		glbDirectedCyl(curVal, curVal.length(), 0.0075, 0.0075);
	}
	rib << "Translate " << curVal << endl;
	curVal.glTranslate();

	SkelTransform::renderSkelRIB(showMarkers, rib);
}

void SkelSymmetricTranslation::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();

	SkelTransform::registerProps(classInfo);

	prop = new SLProperty();
	prop->offset = offsetof(SkelSymmetricTranslation, origTranslation);
	prop->type = SL_CHAR256;
	classInfo->members.addT("origTranslation", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelSymmetricTranslation, axisMult);
	prop->type = SL_VEC3D;
	classInfo->members.addT("axisMult", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelSymmetricTranslation, box);
	prop->type = SL_VEC3D;
	classInfo->members.addT("box", prop);

	classInfo->size = sizeof(SkelSymmetricTranslation);
	classInfo->newInstance = SkelSymmetricTranslation::newInstance;
	sl->classes.addT("SkelSymmetricTranslation", classInfo);
}

SLInterface *SkelSymmetricTranslation::newInstance(int count) {
	return new SkelSymmetricTranslation[count];
}

void SkelSymmetricTranslation::copyFrom(SkelTransform *st) {
	SkelTransform::copyFrom(st);
	strcpy(origTranslation, ((SkelSymmetricTranslation*)st)->origTranslation);
	axisMult = ((SkelSymmetricTranslation*)st)->axisMult;
	curVal = ((SkelSymmetricTranslation*)st)->curVal;
	box = ((SkelSymmetricTranslation*)st)->box;
}

SkelTransform *SkelSymmetricTranslation::clone() {
	SkelTransform *ret = new SkelSymmetricTranslation();
	ret->copyFrom(this);
	return ret;
}


// SkelPartialTransform ======================================

SkelPartialTransform::SkelPartialTransform() {
	className = "SkelPartialTransform";

	origTransform[0] = 0;
	orig = NULL;
	factor = 1;
}

void SkelPartialTransform::initRefs(Skeleton *skel) {
	int origInd = skel->transforms.lookupName(origTransform);
	if (origInd < 0) {
		cout << "warning: " << name << "'s original transform, " << 
			origTransform << ", is unknown." << endl;
		orig = NULL;
	}
	else {
		orig = (SkelTranslation*)skel->transforms.getT(origInd);
	}
}

void SkelPartialTransform::updateCoord() {
	if (!orig)
		return;
	else if (strcmp(orig->className, "SkelEulerRotation") == 0) {
		SkelEulerRotation *rot = (SkelEulerRotation*)orig;
		QuatNorm curQuat;
		switch (rot->axis) {
		case 0:
			curQuat = QuatNorm(sin(factor * rot->curAngle / 2.0), 0, 0,
				-cos(factor * rot->curAngle / 2.0));
			break;
		case 1:
			curQuat = QuatNorm(0, sin(factor * rot->curAngle / 2.0), 0, 
				-cos(factor * rot->curAngle / 2.0));
			break;
		case 2:
			curQuat = QuatNorm(0, 0, sin(factor * rot->curAngle / 2.0), 
				-cos(factor * rot->curAngle / 2.0));
			break;
		}

		curCoord.mat = curQuat.toMatrixD();
		curCoord.q = curQuat;
		curCoord.v = Vec3d();
	}
	else {
		curCoord.q = slerp(QuatNorm(), orig->curCoord.q, factor);
		curCoord.v = factor * orig->curCoord.v;

		curCoord.mat = curCoord.q.toMatrixD();
		curCoord.mat[0][3] = curCoord.v[0];
		curCoord.mat[1][3] = curCoord.v[1];
		curCoord.mat[2][3] = curCoord.v[2];
	}
}

void SkelPartialTransform::updateGlobalDerivs(int maxDof) {
	SkelTransform::updateGlobalDerivs(maxDof);

	if (!orig)
		return;
	else if (strcmp(orig->className, "SkelEulerRotation") == 0) {
		SkelEulerRotation *rot = (SkelEulerRotation*)orig;
		Mat4d derivMat;
		derivMat *= 0;

		double cosTheta = factor * cos(factor * rot->curAngle);
		double sinTheta = factor * sin(factor * rot->curAngle);

		switch(rot->axis) {
		case 0:
			derivMat[1][1] = -sinTheta;
			derivMat[1][2] = -cosTheta;
			derivMat[2][1] = cosTheta;
			derivMat[2][2] = -sinTheta;
			break;
		case 1:
			derivMat[0][0] = -sinTheta;
			derivMat[0][2] = cosTheta;
			derivMat[2][0] = -cosTheta;
			derivMat[2][2] = -sinTheta;
			break;
		case 2:
			derivMat[0][0] = -sinTheta;
			derivMat[0][1] = -cosTheta;
			derivMat[1][0] = cosTheta;
			derivMat[1][1] = -sinTheta;
			break;
		}

		if (parentPtr)
			globalDerivs[orig->dofInd] += parentPtr->globalCoord.mat * derivMat;
		else
			globalDerivs[orig->dofInd] += derivMat;
	}
	else
		cerr << "WARNING: unsupported parent type for SkelPartialTransform" << endl;

	/* incorrect code
	if (orig) {
		int dof;
		for (dof = 0; dof < orig->numDofs(); dof++) {
			if (parentPtr)
				globalDerivs[dof + orig->dofInd] += factor * parentPtr->globalCoord.mat * (*orig->curCoord.deriv[dof]);
			else
				globalDerivs[dof + orig->dofInd] += factor * (*orig->curCoord.deriv[dof]);
		}
	}
	*/
}

void SkelPartialTransform::loadPose(istream &in) {
}

double SkelPartialTransform::distance(SkelTransform *k) {
	return 0; //(curVal - ((SkelPartialTransform*)k)->curVal).length2();
}

void SkelPartialTransform::drawGL(double alpha)  {
	/*if (!color.iszero()) {
		color.glColor(alpha);
		glbDirectedCyl(curVal, curVal.length(), 0.0075, 0.0075);
	}
	curVal.glTranslate();*/
}

void SkelPartialTransform::renderSkelRIB(bool showMarkers, ostream &rib)  {
/*	if (!color.iszero()) {
		rib << "Color " << color << endl;
//		glbDirectedCyl(curVal, curVal.length(), 0.0075, 0.0075);
	}
	rib << "Translate " << curVal << endl;
	curVal.glTranslate();

	SkelTransform::renderSkelRIB(showMarkers, rib);*/
}

void SkelPartialTransform::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();

	SkelTransform::registerProps(classInfo);

	prop = new SLProperty();
	prop->offset = offsetof(SkelPartialTransform, origTransform);
	prop->type = SL_CHAR256;
	classInfo->members.addT("origTransform", prop);

	prop = new SLProperty();
	prop->offset = offsetof(SkelPartialTransform, factor);
	prop->type = SL_DOUBLE;
	classInfo->members.addT("factor", prop);

	classInfo->size = sizeof(SkelPartialTransform);
	classInfo->newInstance = SkelPartialTransform::newInstance;
	sl->classes.addT("SkelPartialTransform", classInfo);
}

SLInterface *SkelPartialTransform::newInstance(int count) {
	return new SkelPartialTransform[count];
}

void SkelPartialTransform::copyFrom(SkelTransform *st) {
	SkelTransform::copyFrom(st);
	strcpy(origTransform, ((SkelPartialTransform*)st)->origTransform);
	factor = ((SkelPartialTransform*)st)->factor;
}

SkelTransform *SkelPartialTransform::clone() {
	SkelTransform *ret = new SkelPartialTransform();
	ret->copyFrom(this);
	return ret;
}


// SkelCombinedTransform =====================================

SkelCombinedTransform::SkelCombinedTransform() {
	int i;
	for (i=0; i < MAX_COMBINED; i++) {
		orig[i] = NULL;
		origTransform[i][0] = 0;
	}
}

void SkelCombinedTransform::initRefs(Skeleton *skel) {
	int i;
	for (i=0; i < MAX_COMBINED; i++) {
		if (strlen(origTransform[i]) > 0) {
			int origInd = skel->transforms.lookupName(origTransform[i]);
			if (origInd < 0) {
				cout << "warning: " << name << "'s original transform, " << 
					origTransform[i]<< ", is unknown." << endl;
				orig[i] = NULL;
			}
			else {
				orig[i] = (SkelTranslation*)skel->transforms.getT(origInd);
			}
		}
	}
}

void SkelCombinedTransform::updateCoord() {
	curCoord.mat = Mat4d();
	curCoord.q = QuatNorm();
	curCoord.v = Vec3d();

	int i;
	for (i=0; i < MAX_COMBINED; i++) {
		if (orig[i])
			curCoord = curCoord * orig[i]->curCoord;
	}
}

void SkelCombinedTransform::updateGlobalDerivs(int maxDof) {
	SkelTransform::updateGlobalDerivs(maxDof);

	int i;
	for (i=0; i < MAX_COMBINED; i++) {
		if (orig[i]) {
			int dof;
			for (dof = 0; dof < orig[i]->numDofs(); dof++) {
				if (parentPtr)
					globalDerivs[dof + orig[i]->dofInd] += parentPtr->globalCoord.mat * (*orig[i]->curCoord.deriv[dof]);
				else
					globalDerivs[dof + orig[i]->dofInd] += (*orig[i]->curCoord.deriv[dof]);
			}
		}
	}
}

void SkelCombinedTransform::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();

	SkelTransform::registerProps(classInfo);

	int i;
	for (i=0; i < MAX_COMBINED; i++) {
		char name[256];

		sprintf(name, "origTransform%d", i);
		prop = new SLProperty();
		prop->offset = (int)offsetof(SkelCombinedTransform, origTransform[i]);
		prop->type = SL_CHAR256;
		classInfo->members.addT(name, prop);
	}

	classInfo->size = sizeof(SkelCombinedTransform);
	classInfo->newInstance = SkelCombinedTransform::newInstance;
	sl->classes.addT("SkelCombinedTransform", classInfo);
}

SLInterface *SkelCombinedTransform::newInstance(int count) {
	return new SkelCombinedTransform[count];
}

void SkelCombinedTransform::copyFrom(SkelTransform *s) {
	SkelTransform::copyFrom(s);
	int i;
	for (i=0; i < MAX_COMBINED; i++) {
		strcpy(origTransform[i], ((SkelCombinedTransform*)s)->origTransform[i]);
	}
}

SkelTransform *SkelCombinedTransform::clone() {
	SkelTransform *ret = new SkelCombinedTransform();
	ret->copyFrom(this);
	return ret;
}