#include <fstream>
#include "kinematics.h"

// KCoord =================================================

KCoord::KCoord() {
	deriv = NULL;
	numDofs = -1;
}

KCoord::~KCoord() {
	if (deriv) {
		int i;
		for (i=0; i < numDofs; i++)
			if (deriv[i])
				delete deriv[i];
		delete []deriv;
	}
}

void KCoord::initDeriv(int iNumDofs) {
	if (numDofs == iNumDofs)
		return;

	numDofs = iNumDofs;
	deriv = new Mat4dPtr[numDofs];
	int i;
	for (i=0; i < numDofs; i++)
		deriv[i] = 0;
}

void KCoord::setDeriv(int dof, Mat4d &val) {
	if (!deriv || dof < 0 || dof >= numDofs) {
		cerr << "warning: illegal derivative set -- dof " << dof << endl;
		return;
	}

	if (!deriv[dof])
		deriv[dof] = new Mat4d();
	*deriv[dof] = val;
}

KCoord KCoord::operator*(const KCoord& a) {
	KCoord ret;

	ret.mat = mat * a.mat;
	ret.q = a.q * q;
	ret.v = v + q.toMatrixD() * a.v;
//	ret.v = Vec3d(ret.mat[0][3], ret.mat[1][3], ret.mat[2][3]);

	return ret;
}

KCoord& KCoord::operator=(const KCoord& a) {
	mat = a.mat;
	q = a.q;
	v = a.v;
	return *this;
}













#ifdef OLD_VERSION

const int KTranslation::CONSTANT_X = 1;
const int KTranslation::CONSTANT_Y = 2;
const int KTranslation::CONSTANT_Z = 4;

// KDofSet ================================================

KDofSet::KDofSet() {
	globalSize = 0;
	perFrameSize = 0;
	numFrames = 0;
	numDofs = 0;
	nonMarkerSize = 0;
}

KDofSet::~KDofSet() {
	flush();
}

void KDofSet::flush() {
//	if (data)
//		delete []data;
//	data = NULL;
}

KDofSet *KDofSet::clone() {
	KDofSet *ret = new KDofSet();

//	if (data) {
//		ret->data = new double[numDofs];
//		memcpy(ret->data, data, numDofs * sizeof(double));
//	}
	ret->data.resize(data.size());
	int i;
	for (i=0; i < data.size(); i++)
		ret->data[i] = data[i];

	ret->globalSize = globalSize;
	ret->perFrameSize = perFrameSize;
	ret->numFrames = numFrames;
	ret->numDofs = numDofs;
	ret->nonMarkerSize = nonMarkerSize;
	
	return ret;
}

void KDofSet::copyVals(KDofSet *d) {
	if (data.size() != d->data.size()) {
		cout << "error in KDofSet::copyVals: data size mismatch" << endl;
		return;
	}
	else
		data = d->data;
}

int KDofSet::calcIndex(int frame, int ofs) {
	if (ofs < globalSize)
		return ofs;
	else
		return globalSize + perFrameSize * frame + (ofs - globalSize);
}

double &KDofSet::v(int ofs) {
	return data[ofs];
}

double &KDofSet::v(int frame, int ofs) {
	return data[calcIndex(frame, ofs)];
}

ostream& operator <<(ostream& os, KDofSet &k) {
	int i;
	os << k.globalSize << endl;
	os << k.perFrameSize << endl;
	os << k.numFrames << endl;
	os << k.numDofs << endl;
	os << k.nonMarkerSize << endl;
	for (i=0; i < k.numDofs; i++)
		os << k.data[i] << " ";
	os << endl;
	return os;
}

istream& operator >>(istream& is, KDofSet &k) {
	int i;

	is >> k.globalSize;
	is >> k.perFrameSize;
	is >> k.numFrames;
	is >> k.numDofs;
	is >> k.nonMarkerSize;
	
	k.flush();
	k.data.resize(k.numDofs);
	//k.data = new double[k.numDofs];

	for (i=0; i < k.numDofs; i++)
		is >> k.data[i];

	return is;
}

// KMarkerData ============================================

KMarkerData::KMarkerData() {
	data = NULL;
	colors = NULL;
	numMarkers = 0;
	numFrames = 0;

	variableSize = false;
	colorData = false;
}

KMarkerData::~KMarkerData() {
	flush();
}

void KMarkerData::flush() {
	if (data)
		delete []data;
	data = NULL;

	if (colors)
		delete []colors;
	colors = NULL;
}

void KMarkerData::init(int nMarkers, int nFrames) {
	flush();

	numMarkers = nMarkers;
	numFrames = nFrames;

	data = new Vec3d[numMarkers * numFrames];

	if (colorData)
		colors = new Vec3d[numMarkers * numFrames];
}

Vec3d &KMarkerData::v(int frame, int marker) {
	return data[frame * numMarkers + marker];
}

Vec3d &KMarkerData::c(int frame, int marker) {
	static Vec3d defaultColor(128, 128, 128);

	if (!colors)
		return defaultColor;
	else
		return colors[frame * numMarkers + marker];
}

ostream& operator <<(ostream& os, KMarkerData &k) {
	int i, j;

	os << k.numMarkers << endl;
	os << k.numFrames << endl;

	for (i=0; i < k.numFrames; i++) {
		os << i << endl;

		if (k.variableSize)
			os << k.numMarkers << endl;

		for (j=0; j < k.numMarkers; j++) {
			os << j << " " << k.v(i, j);
			if (k.colorData)
				os << " " << (k.c(i, j) * 255.0);
			os << endl;
		}
	}
	return os;
}

istream& operator >>(istream& is, KMarkerData &k) {
	int i, j, junk, curCount;

	is >> i;
	is >> j;

	k.init(i, j);
	curCount = k.numMarkers;

	for (i=0; i < k.numFrames; i++) {
		is >> junk;
		if (junk != i)
			cout << "problem reading marker data -- expected " << i << ", got " << junk << endl;

		if (k.variableSize)
			is >> curCount;

		for (j=0; j < curCount; j++) {
			is >> junk >> k.v(i, j);
			if (k.colorData) {
				is >> k.c(i, j);
				k.c(i, j) /= 255.0;
			}
		}
	}
	return is;
}

// KMarker ================================================

void KMarker::loadDofs(KDofSet *dofs) {
	curVal[0] = dofs->v(dofIndex+0);
	curVal[1] = dofs->v(dofIndex+1);
	curVal[2] = dofs->v(dofIndex+2);
}

void KMarker::unloadDofs(KDofSet *dofs) {
	dofs->v(dofIndex+0) = curVal[0];
	dofs->v(dofIndex+1) = curVal[1];
	dofs->v(dofIndex+2) = curVal[2];
}

void KMarker::updatePos() {
	curPos = model->transforms[transform]->curCoord.mat * curVal;
}

ostream& operator <<(ostream& os, KMarker &k) {
	os << k.dofIndex << " " << k.markerIndex << " " << k.transform << endl;
	os << k.curVal << endl;
	return os;
}

istream& operator >>(istream& is, KMarker &k) {
	is >> k.dofIndex >> k.markerIndex >> k.transform;
	is >> k.curVal;
	return is;
}



// KTransform =============================================

KTransform::KTransform() {
	model = NULL;
	parent = -1;
	perFrame = false;
	hi = 1;
	lo = 0;
	importance = 1;
	isConst = true;
	isCopy = false;
}

KTransform::~KTransform() {
}

void KTransform::addChild(int i) {
	children.push_back(i);
	model->transforms[i]->parent = index;
}

int KTransform::numTotalDofs() {
	return numGDofs() + numFDofs();
}

void KTransform::loadAllDofs(KDofSet *dofs, int frame) {
	loadDofs(dofs, frame);
	int i;
	for (i=0; i < children.size(); i++)
		model->transforms[children[i]]->loadAllDofs(dofs, frame);
}

void KTransform::unloadAllDofs(KDofSet *dofs, int frame) {
	unloadDofs(dofs, frame);
	int i;
	for (i=0; i < children.size(); i++)
		model->transforms[children[i]]->unloadAllDofs(dofs, frame);
}

void KTransform::calcRecursiveTransforms(KCoord *coords) {
	if (parent == -1)
		coords[index] = curCoord;
	else
		coords[index] = coords[parent] * curCoord;

	int i;
	for (i=0; i < children.size(); i++)
		model->transforms[children[i]]->calcRecursiveTransforms(coords);
}

void KTransform::calcRecursiveDerivs(KCoord *coords, int frame) {
	Mat4d prev;
	int i;
	int trueGDofIndex, trueFDofIndex;

	if (parent == -1) {
		coords[index] = curCoord;
	}
	else {
		coords[index] = coords[parent] * curCoord;
		prev = coords[parent].mat;
	}

	trueGDofIndex = model->dofs->calcIndex(frame, dofGIndex);
	trueFDofIndex = model->dofs->calcIndex(frame, dofFIndex);
	for (i=0; i < model->dofs->numDofs; i++) {
		if (i >= trueGDofIndex && i < trueGDofIndex + numGDofs()) {
			coords[index].setDeriv(i, prev * (*curCoord.deriv[i - trueGDofIndex]));
		}
		else if (i >= trueFDofIndex && i < trueFDofIndex + numFDofs()) {
			coords[index].setDeriv(i, prev * (*curCoord.deriv[numGDofs() + i - trueFDofIndex]));
		}
		else {
			if (parent > -1 && coords[parent].deriv[i] != NULL) {
				coords[index].setDeriv(i, (*coords[parent].deriv[i]) * curCoord.mat);
			}
		}
	}

	for (i=0; i < children.size(); i++)
		model->transforms[children[i]]->calcRecursiveDerivs(coords, frame);
}

void KTransform::loadDofs(KDofSet *dofs, int frame) {
	// this standard behavior works for transforms that are all global or 
	// all per-frame, i.e., everything except KPolarAxisRotation
	if (perFrame)
		loadRaw(&dofs->v(frame, dofFIndex));
	else
		loadRaw(&dofs->v(frame, dofGIndex));
}

void KTransform::unloadDofs(KDofSet *dofs, int frame) {
	// this standard behavior works for transforms that are all global or 
	// all per-frame, i.e., everything except KPolarAxisRotation
	if (perFrame)
		unloadRaw(&dofs->v(frame, dofFIndex));
	else
		unloadRaw(&dofs->v(frame, dofGIndex));
}

void KTransform::put(ostream &os) {
	int i;

	os << name << endl;
	os << index << " " << dofFIndex << " " << dofGIndex << " " << perFrame << endl;
	os << color << endl;

	os << parent << endl;
	os << children.size() << " ";
	for (i=0; i < children.size(); i++) {
		os << children[i] << " ";
	}
	os << endl;
	putVal(os);
	os << endl;
}

void KTransform::get(istream &is) {
	int i, count, cur;

	is.getline(name, 40);
	is >> index >> dofFIndex >> dofGIndex >> perFrame;
	is >> color;

	is >> parent;
	is >> count;
	for (i=0; i < count; i++) {
		is >> cur;
		children.push_back(cur);
	}
	getVal(is);
}

void KTransform::putVal(ostream &os) {
}

void KTransform::getVal(istream &is) {
}

ostream& operator <<(ostream& os, KTransform &k) {
	k.put(os);
	return os;
}

istream& operator >>(istream& is, KTransform &k) {
	k.get(is);
	return is;
}

void KTransform::renderSkel(bool showMarkers, bool fancy) {
	if (showMarkers) {
		int i;
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
		}
	}
}

void KTransform::renderSkelRIB(bool showMarkers, ostream &rib) {
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

void KTransform::renderAllSkel(bool showMarkers, bool fancy) {
	glPushMatrix();
		renderSkel(showMarkers, fancy);
		int i;
		for (i=0; i < children.size(); i++)
			model->transforms[children[i]]->renderAllSkel(showMarkers, fancy);
	glPopMatrix();
}

void KTransform::renderAllSkelRIB(bool showMarkers, ostream &rib) {
	rib << "TransformBegin" << endl;
		renderSkelRIB(showMarkers, rib);
		int i;
		for (i=0; i < children.size(); i++)
			model->transforms[children[i]]->renderAllSkelRIB(showMarkers, rib);
	rib << "TransformEnd" << endl;
}

void KTransform::updateUI(Fl_Value_Input *val0, Fl_Value_Input *val1, 
					  Fl_Value_Input *val2, Fl_Value_Input *val3, Fl_Output *out, int frame) {
}

double KTransform::scaleFromDof(double x) {
	return lo + (hi - lo) * x / importance;
}

double KTransform::scaleToDof(double x) {
	return (x - lo) / (hi - lo) * importance;
}

double KTransform::scaleDeriv() {
	return (hi - lo) / importance;

}

// KTranslation ===========================================

KTranslation::KTranslation() {
	constAxes = 0;
	constVal = Vec3d(0, 0, 0);
	curVal = constVal;
}

int KTranslation::type() {
	return KT_TRANSLATION;
}

int KTranslation::numGDofs() {
	if (isConst || perFrame)
		return 0;

	int d = 3;
	if (constAxes & CONSTANT_X)
		d--;
	if (constAxes & CONSTANT_Y)
		d--;
	if (constAxes & CONSTANT_Z)
		d--;
	return d;
}

int KTranslation::numFDofs() {
	if (isConst || !perFrame)
		return 0;

	int d = 3;
	if (constAxes & CONSTANT_X)
		d--;
	if (constAxes & CONSTANT_Y)
		d--;
	if (constAxes & CONSTANT_Z)
		d--;
	return d;
}

void KTranslation::loadRaw(double *raw) {
	int di = 0;

	curVal = constVal;

	if (isConst)
		return;
	if ((constAxes & CONSTANT_X) == 0)
		curVal[0] = scaleFromDof(raw[di++]);
	if ((constAxes & CONSTANT_Y) == 0)
		curVal[1] = scaleFromDof(raw[di++]);
	if ((constAxes & CONSTANT_Z) == 0)
		curVal[2] = scaleFromDof(raw[di++]);
}

void KTranslation::unloadRaw(double *raw) {
	int di = 0;

	if (isConst)
		return;
	if ((constAxes & CONSTANT_X) == 0)
		raw[di++] = scaleToDof(curVal[0]);
	if ((constAxes & CONSTANT_Y) == 0)
		raw[di++] = scaleToDof(curVal[1]);
	if ((constAxes & CONSTANT_Z) == 0)
		raw[di++] = scaleToDof(curVal[2]);
}

void KTranslation::updateTransform() {
	curCoord.mat = Mat4d();
	curCoord.mat[0][3] = curVal[0];
	curCoord.mat[1][3] = curVal[1];
	curCoord.mat[2][3] = curVal[2];

	curCoord.q = QuatNorm();
	curCoord.v = curVal;
}

void KTranslation::updateDerivs() {
	int di = 0;
	Mat4d m;

	curCoord.initDeriv(numTotalDofs());

	if (isConst)
		return;

	if ((constAxes & CONSTANT_X) == 0) {
		m *= 0;
		m[0][3] = scaleDeriv();
		curCoord.setDeriv(di++, m);
	}
	if ((constAxes & CONSTANT_Y) == 0) {
		m *= 0;
		m[1][3] = scaleDeriv();
		curCoord.setDeriv(di++, m);
	}
	if ((constAxes & CONSTANT_Z) == 0) {
		m *= 0;
		m[2][3] = scaleDeriv();
		curCoord.setDeriv(di++, m);
	}
}

void KTranslation::copyVal(KTransform *k) {
	constVal = ((KTranslation*)k)->curVal;
	curVal = ((KTranslation*)k)->curVal;
}

void KTranslation::putVal(ostream &os) {
	os << curVal;
}

void KTranslation::getVal(istream &is) {
	is >> curVal;
	if (isConst)
		constVal = curVal;
}

void KTranslation::renderSkel(bool showMarkers, bool fancy)  {
	if (!color.iszero()) {
		color.glColor();
		if (fancy) {
			glbDirectedCyl(curVal, curVal.length(), 0.005, 0.005);
		}
		else {
			glDisable(GL_LIGHTING);
			glLineWidth(3);
			glBegin(GL_LINES);
				glVertex3d(0, 0, 0);
				curVal.glVertex();
			glEnd();
		}
	}
	curVal.glTranslate();
	KTransform::renderSkel(showMarkers, fancy);
}

void KTranslation::renderSkelRIB(bool showMarkers, ostream &rib)  {
	if (!color.iszero()) {
		rib << "Color [ " << color << " ]" << endl;
//		glbDirectedCyl(curVal, curVal.length(), 0.005, 0.005);
	}
	rib << "Translate " << curVal << endl;
	KTransform::renderSkelRIB(showMarkers, rib);
}

void KTranslation::updateUI(Fl_Value_Input *val0, Fl_Value_Input *val1, 
					  Fl_Value_Input *val2, Fl_Value_Input *val3, Fl_Output *out, int frame) {
	int di;
	if (perFrame)
		di = dofFIndex;
	else
		di = dofGIndex;

	val0->activate();
	val0->label("X:");
	if (isConst || (constAxes & CONSTANT_X))
		val0->value(constVal[0]);
	else
		val0->value(scaleFromDof(model->dofs->v(frame, di++)));
	val1->activate();
	val1->label("Y:");
	if (isConst || (constAxes & CONSTANT_Y))
		val1->value(constVal[1]);
	else
		val1->value(scaleFromDof(model->dofs->v(frame, di++)));
	val2->activate();
	val2->label("Z:");
	if (isConst || (constAxes & CONSTANT_Z))
		val2->value(constVal[2]);
	else
		val2->value(scaleFromDof(model->dofs->v(frame, di++)));
	val3->deactivate();
	out->activate();
	if (isConst)
		out->value("translation (constant)");
	else
		out->value("translation");
}

// KQuatRotation ==========================================

int KQuatRotation::type() {
	return KT_QUATROTATION;
}

int KQuatRotation::numGDofs() {
	if (isConst || perFrame)
		return 0;
	return 4;
}

int KQuatRotation::numFDofs() {
	if (isConst || !perFrame)
		return 0;
	return 4;
}

void KQuatRotation::loadRaw(double *raw) {
	if (isConst) {
		curQuat = constQuat;
		return;
	}
	curQuat.x = scaleFromDof(raw[0]);
	curQuat.y = scaleFromDof(raw[1]);
	curQuat.z = scaleFromDof(raw[2]);
	curQuat.w = scaleFromDof(raw[3]);
}

void KQuatRotation::unloadRaw(double *raw) {
	if (isConst)
		return;
	raw[0] = scaleToDof(curQuat.x);
	raw[1] = scaleToDof(curQuat.y);
	raw[2] = scaleToDof(curQuat.z);
	raw[3] = scaleToDof(curQuat.w);
}

void KQuatRotation::updateTransform() {
	curCoord.mat = curQuat.toMatrixD();

	curCoord.q = curQuat;
	curCoord.v = Vec3d();
}

void KQuatRotation::updateDerivs() {
	curCoord.initDeriv(numTotalDofs());

	if (isConst)
		return;

	Mat4d val;
	Mat4d curDerivs[4];
	curQuat.getMatrices(val, curDerivs);

	int i;
	for (i=0; i < 4; i++)
		curCoord.setDeriv(i, curDerivs[i] * scaleDeriv());
}

void KQuatRotation::copyVal(KTransform *k) {
	constQuat = ((KQuatRotation*)k)->curQuat;
	curQuat = ((KQuatRotation*)k)->curQuat;
}

void KQuatRotation::normalize() {
	curQuat.normalize();
}

void KQuatRotation::putVal(ostream &os) {
	os << curQuat;
}

void KQuatRotation::getVal(istream &is) {
	is >> curQuat;
	if (isConst)
		constQuat = curQuat;
}

void KQuatRotation::renderSkel(bool showMarkers, bool fancy)  {
	if (fancy) {
		glColor3d(0.5, 0.5, 0.5);
		glbSphere(Vec3d(0, 0, 0), 0.02);
	}

	double mat[16];
	curQuat.toMatrixD().getGLMatrix(mat);
	glMultMatrixd(mat);

	KTransform::renderSkel(showMarkers, fancy);
}

void KQuatRotation::renderSkelRIB(bool showMarkers, ostream &rib)  {
	rib << "Color 0.5 0.5 0.5" << endl;
	rib << "Sphere 0.02 -0.02 0.02 360" << endl;

	curQuat.saveRIB(rib);

	KTransform::renderSkelRIB(showMarkers, rib);
}

void KQuatRotation::updateUI(Fl_Value_Input *val0, Fl_Value_Input *val1, 
					  Fl_Value_Input *val2, Fl_Value_Input *val3, Fl_Output *out, int frame) {
	int di;
	if (perFrame)
		di = dofFIndex;
	else
		di = dofGIndex;

	val0->activate();
	val0->label("X:");
	val0->value(scaleFromDof(model->dofs->v(frame, di+0)));
	val1->activate();
	val1->label("Y:");
	val1->value(scaleFromDof(model->dofs->v(frame, di+1)));
	val2->label("Z:");
	val2->activate();
	val2->value(scaleFromDof(model->dofs->v(frame, di+2)));
	val3->label("W:");
	val3->activate();
	val3->value(scaleFromDof(model->dofs->v(frame, di+3)));
	out->activate();
	if (isConst)
		out->value("quaternion (constant)");
	else
		out->value("quaternion");
}

// KEulerRotation =========================================

KEulerRotation::KEulerRotation() {
	axis = 0;
	curAngle = 0;
	constAngle = 0;

	min = -PI;
	max = PI*3.0/2.0;
}

int KEulerRotation::type() {
	return KT_EULERROTATION;
}

int KEulerRotation::numGDofs() {
	if (isConst || perFrame)
		return 0;
	return 1;
}

int KEulerRotation::numFDofs() {
	if (isConst || !perFrame)
		return 0;
	return 1;
}

void KEulerRotation::loadRaw(double *raw) {
	if (isConst)
		curAngle = constAngle;
	else
		curAngle = scaleFromDof(raw[0]);

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
}

void KEulerRotation::unloadRaw(double *raw) {
	if (!isConst)
		raw[0] = scaleToDof(curAngle);
}

void KEulerRotation::updateTransform() {
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

void KEulerRotation::updateDerivs() {
	if (isConst)
		return;

	Mat4d derivMat;
	derivMat *= 0;

	curCoord.initDeriv(numTotalDofs());

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

	curCoord.setDeriv(0, derivMat * scaleDeriv());
}

void KEulerRotation::copyVal(KTransform *k) {
	constAngle = ((KEulerRotation*)k)->curAngle;
	curAngle = ((KEulerRotation*)k)->curAngle;
}

void KEulerRotation::normalize() {
	while (curAngle < min)
		curAngle += 2.0*PI;
	while (curAngle > max)
		curAngle -= 2.0*PI;
}

void KEulerRotation::putVal(ostream &os) {
	os << axis << " " << curAngle;
}

void KEulerRotation::getVal(istream &is) {
	is >> axis >> curAngle;
	if (isConst)
		constAngle = curAngle;
}

void KEulerRotation::renderSkel(bool showMarkers, bool fancy)  {
	if (fancy) {
		glColor3d(0.5, 0, 0.5);
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

	KTransform::renderSkel(showMarkers, fancy);
}

void KEulerRotation::renderSkelRIB(bool showMarkers, ostream &rib)  {
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

	KTransform::renderSkelRIB(showMarkers, rib);
}

void KEulerRotation::updateUI(Fl_Value_Input *val0, Fl_Value_Input *val1, 
					  Fl_Value_Input *val2, Fl_Value_Input *val3, Fl_Output *out, int frame) {
	val0->activate();
	val0->label("axis:");
	val0->value(axis);
	val1->activate();
	val1->label("angle:");
	if (perFrame)
		val1->value(scaleFromDof(model->dofs->v(frame, dofFIndex)));
	else
		val1->value(scaleFromDof(model->dofs->v(frame, dofGIndex)));
	val2->label("");
	val2->deactivate();
	val3->label("");
	val3->deactivate();
	out->activate();
	if (isConst)
		out->value("euler angle (constant)");
	else
		out->value("euler angle");
}

// KPolarAxisRotation =========================================

KPolarAxisRotation::KPolarAxisRotation() {
	curPhi = 0;
	curPsi = 0;
	curTheta = 0;

	constAxis = false;

	min = -PI;
	max = PI*3.0/2.0;
}

int KPolarAxisRotation::type() {
	return KT_POLARAXISROTATION;
}

int KPolarAxisRotation::numGDofs() {
	if (isConst || constAxis)
		return 0;
	if (perFrame)
		return 2;
	return 3;
}

int KPolarAxisRotation::numFDofs() {
	if (isConst)
		return 0;
	if (perFrame)
		return 1;
	return 0;
}

void KPolarAxisRotation::loadRaw(double *raw) {
	cout << "warning: shouldn't call loadRaw on a KPolarAxisRotation!" << endl;
}

void KPolarAxisRotation::loadDofs(KDofSet *dofs, int frame) {
	if (isConst) {
		curPhi = constPhi;
		curPsi = constPsi;
		curTheta = constTheta;
	}
	else {
		if (constAxis) {
			curPhi = constPhi;
			curPsi = constPsi;
		}
		else {
			curPhi = scaleFromDof(dofs->v(frame, dofGIndex + 0));
			curPsi = scaleFromDof(dofs->v(frame, dofGIndex + 1));
		}
		if (perFrame)
			curTheta = scaleFromDof(dofs->v(frame, dofFIndex));
		else
			curTheta = scaleFromDof(dofs->v(frame, dofGIndex + 2));
	}

	curAxis = Vec3d(sin(curPhi) * cos(curPsi), sin(curPhi) * sin(curPsi), cos(curPhi));
	double sa = sin(curTheta / 2.0);
	curQuat = QuatNorm(sa * curAxis[0], sa * curAxis[1], sa * curAxis[2], -cos(curTheta / 2.0));
}

void KPolarAxisRotation::unloadRaw(double *raw) {
	cout << "warning: shouldn't call unloadRaw on a KPolarAxisRotation!" << endl;
}

void KPolarAxisRotation::unloadDofs(KDofSet *dofs, int frame) {
	if (!isConst) {
		if (!constAxis) {
			dofs->v(frame, dofGIndex + 0) = scaleToDof(curPhi);
			dofs->v(frame, dofGIndex + 1) = scaleToDof(curPsi);
		}
		if (perFrame)
			dofs->v(frame, dofFIndex) = scaleToDof(curTheta);
		else
			dofs->v(frame, dofGIndex + 2) = scaleToDof(curTheta);
	}
}

void KPolarAxisRotation::updateTransform() {
	curAxis = Vec3d(sin(curPhi) * cos(curPsi), sin(curPhi) * sin(curPsi), cos(curPhi));
	double sa = sin(curTheta / 2.0);
	curQuat = QuatNorm(sa * curAxis[0], sa * curAxis[1], sa * curAxis[2], -cos(curTheta / 2.0));

	curCoord.mat = curQuat.toMatrixD();
	curCoord.q = curQuat;
	curCoord.v = Vec3d();
}

void KPolarAxisRotation::updateDerivs() {
	if (isConst)
		return;

	Mat4d derivMat;
	derivMat *= 0;

	curCoord.initDeriv(numTotalDofs());

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
	if (constAxis) {
		curCoord.setDeriv(0, curDeriv * scaleDeriv());
		return;
	}
	else
		curCoord.setDeriv(2, curDeriv * scaleDeriv());
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
	// phi
	curCoord.setDeriv(1, curDeriv * scaleDeriv());
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
	curCoord.setDeriv(0, curDeriv * scaleDeriv());
}

void KPolarAxisRotation::copyVal(KTransform *k) {
	constPhi = ((KPolarAxisRotation*)k)->curPhi;
	curPhi = ((KPolarAxisRotation*)k)->curPhi;
	constPsi = ((KPolarAxisRotation*)k)->curPsi;
	curPsi = ((KPolarAxisRotation*)k)->curPsi;
	constTheta = ((KPolarAxisRotation*)k)->curTheta;
	curTheta = ((KPolarAxisRotation*)k)->curTheta;
}

void KPolarAxisRotation::normalize() {
	while (curTheta < min)
		curTheta += 2.0*PI;
	while (curTheta > max)
		curTheta -= 2.0*PI;
}

void KPolarAxisRotation::putVal(ostream &os) {
	os << curPhi << " " << curPsi << " " << curTheta;
}

void KPolarAxisRotation::getVal(istream &is) {
	is >> curPhi >> curPsi >> curTheta;
	if (isConst || constAxis) {
		constPhi = curPhi;
		constPsi = curPsi;
	}
	constTheta = curTheta;
}

void KPolarAxisRotation::renderSkel(bool showMarkers, bool fancy)  {
	if (fancy) {
		glColor3d(0.5, 0.5, 0.5);
		glbDirectedCyl(curAxis, 0.06, 0.02, 0);
	}

	glRotated(RAD_TO_DEG * curTheta, curAxis[0], curAxis[1], curAxis[2]);

	KTransform::renderSkel(showMarkers, fancy);
}

void KPolarAxisRotation::renderSkelRIB(bool showMarkers, ostream &rib)  {
	rib << "Color 0.5 0.5 0.5" << endl;
//		glbDirectedCyl(curAxis, 0.06, 0.02, 0);

	rib << "Rotate " << (RAD_TO_DEG * curTheta) << " " << curAxis << endl;

	KTransform::renderSkelRIB(showMarkers, rib);
}

void KPolarAxisRotation::updateUI(Fl_Value_Input *val0, Fl_Value_Input *val1, 
					  Fl_Value_Input *val2, Fl_Value_Input *val3, Fl_Output *out, int frame) {
	val0->activate();
	val0->label("phi:");
	val0->value(scaleFromDof(model->dofs->v(frame, dofGIndex + 0)));
	val1->activate();
	val1->label("psi:");
	val1->value(scaleFromDof(model->dofs->v(frame, dofGIndex + 1)));
	val2->label("theta:");
	if (perFrame)
		val2->value(scaleFromDof(model->dofs->v(frame, dofFIndex)));
	else
		val2->value(scaleFromDof(model->dofs->v(frame, dofGIndex + 2)));
	val3->label("");
	val3->deactivate();
	out->activate();
	if (isConst)
		out->value("polar axis (constant)");
	else if (constAxis)
		out->value("polar axis (constant axis)");
	else
		out->value("polar axis");
}

// KSymmetricTranslation ==================================

KSymmetricTranslation::KSymmetricTranslation() {
	origTranslation = -1;
	axisMult[0] = 1;
	axisMult[1] = 1;
	axisMult[2] = 1;
	isCopy = true;
}

KTranslation *KSymmetricTranslation::orig() {
	if (origTranslation == -1) {
		cout << "warning: uninitialized symmetric translation" << endl;
		return NULL;
	}
	return (KTranslation*)model->transforms[origTranslation];
}

void KSymmetricTranslation::updateFromOrig() {
	perFrame = orig()->perFrame;
	hi = orig()->hi;
	lo = orig()->lo;
	importance = orig()->importance;
	isConst = orig()->isConst;
	dofFIndex = orig()->dofFIndex;
	dofGIndex = orig()->dofGIndex;
}

int KSymmetricTranslation::type() {
	return KT_SYMMETRICTRANSLATION;
}

int KSymmetricTranslation::numGDofs() {
	return orig()->numGDofs();
}

int KSymmetricTranslation::numFDofs() {
	return orig()->numFDofs();
}

void KSymmetricTranslation::loadRaw(double *raw) {
}

void KSymmetricTranslation::unloadRaw(double *raw) {
}

void KSymmetricTranslation::updateTransform() {
	Vec3d curVal = orig()->curVal;

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

void KSymmetricTranslation::updateDerivs() {
	int di = 0;
	Mat4d m;

	curCoord.initDeriv(numTotalDofs());

	if (isConst)
		return;

	if ((orig()->constAxes & KTranslation::CONSTANT_X) == 0) {
		m *= 0;
		m[0][3] = scaleDeriv() * axisMult[0];
		curCoord.setDeriv(di++, m);
	}
	if ((orig()->constAxes & KTranslation::CONSTANT_Y) == 0) {
		m *= 0;
		m[1][3] = scaleDeriv() * axisMult[1];
		curCoord.setDeriv(di++, m);
	}
	if ((orig()->constAxes & KTranslation::CONSTANT_Z) == 0) {
		m *= 0;
		m[2][3] = scaleDeriv() * axisMult[2];
		curCoord.setDeriv(di++, m);
	}
}

void KSymmetricTranslation::copyVal(KTransform *k) {
}

void KSymmetricTranslation::putVal(ostream &os) {
	os << origTranslation << " " << axisMult[0] << " " << axisMult[1] << " " << axisMult[2];
}

void KSymmetricTranslation::getVal(istream &is) {
	is >> origTranslation >> axisMult[0] >> axisMult[1] >> axisMult[2];
}

void KSymmetricTranslation::renderSkel(bool showMarkers, bool fancy)  {
	Vec3d curVal = orig()->curVal;
	curVal[0] *= axisMult[0];
	curVal[1] *= axisMult[1];
	curVal[2] *= axisMult[2];

	if (!color.iszero()) {
		color.glColor();
		if (fancy) {
			glbDirectedCyl(curVal, curVal.length(), 0.0075, 0.0075);
		}
		else {
			glDisable(GL_LIGHTING);
			glLineWidth(3);
			glBegin(GL_LINES);
				glVertex3d(0, 0, 0);
				curVal.glVertex();
			glEnd();
		}
	}
	curVal.glTranslate();
	KTransform::renderSkel(showMarkers, fancy);
}

void KSymmetricTranslation::renderSkelRIB(bool showMarkers, ostream &rib)  {
	Vec3d curVal = orig()->curVal;
	curVal[0] *= axisMult[0];
	curVal[1] *= axisMult[1];
	curVal[2] *= axisMult[2];

	if (!color.iszero()) {
		rib << "Color " << color << endl;
//		glbDirectedCyl(curVal, curVal.length(), 0.0075, 0.0075);
	}
	rib << "Translate " << curVal << endl;
	curVal.glTranslate();

	KTransform::renderSkelRIB(showMarkers, rib);
}

void KSymmetricTranslation::updateUI(Fl_Value_Input *val0, Fl_Value_Input *val1, 
					  Fl_Value_Input *val2, Fl_Value_Input *val3, Fl_Output *out, int frame) {
	int di;
	if (perFrame)
		di = dofFIndex;
	else
		di = dofGIndex;

	val0->activate();
	val0->label("orig:");
	val0->value(origTranslation);
	val1->deactivate();
	val1->label("");
	val2->deactivate();
	val2->label("");
	val3->deactivate();
	out->activate();
	out->value("symmetric translation");
}

// KinematicModel =========================================

KinematicModel::KinematicModel() {
	coordData = NULL;
}

void KinematicModel::flush() {
	flushCoords();

	if (meshNames) {
		delete []meshNames;
		meshNames = NULL;
	}
}

void KinematicModel::initMT(int nMeshes, int nMarkers, int nTransforms) {
	numMeshes = nMeshes;
	numMarkers = nMarkers;
	numTransforms = nTransforms;

	markers = new KMarker[numMarkers];
	transforms = new KTransformPtr[numTransforms];

	dofs = NULL;
	markerData = NULL;

	curMarkerPos.init(numMarkers, numMeshes);

	meshNames = new Str80[numMeshes];
}

void KinematicModel::calcTransforms() {
	int i, frame;

	for (i=0; i < numMarkers; i++)
		markers[i].loadDofs(dofs);

	for (frame = 0; frame < numMeshes; frame++) {
		for (i=0; i < numTransforms; i++) {
			transforms[i]->loadDofs(dofs, frame);
			transforms[i]->updateTransform();
		}

		transforms[0]->calcRecursiveTransforms(coordData + frame * numTransforms);

		for (i=0; i < numMarkers; i++)
			curMarkerPos.v(frame, i) = trans(markers[i].transform, frame) * markers[i].curVal;
	}
}

void KinematicModel::calcDerivs() {
	int i, frame;

	for (i=0; i < numMarkers; i++)
		markers[i].loadDofs(dofs);

	for (frame = 0; frame < numMeshes; frame++) {
		for (i=0; i < numTransforms; i++) {
			transforms[i]->loadDofs(dofs, frame);
			transforms[i]->updateTransform();
			transforms[i]->updateDerivs();
		}

		transforms[0]->calcRecursiveTransforms(coordData + frame * numTransforms);
		transforms[0]->calcRecursiveDerivs(coordData + frame * numTransforms, frame);

		for (i=0; i < numMarkers; i++) {
			curMarkerPos.v(frame, i) = trans(markers[i].transform, frame) * markers[i].curVal;
		}
	}
}

void KinematicModel::putKin(ostream &os) {
	os << numMeshes << " " << numMarkers << " " << numTransforms << endl;

	int i;
	for (i=0; i < numMarkers; i++)
		os << markers[i] << endl;
	for (i=0; i < numTransforms; i++) {
		os << transforms[i]->type() << endl;
		os << (*transforms[i]) << endl;
	}
}

void KinematicModel::getKin(istream &is) {
	int i, j, k;

	is >> i >> j >> k;

	initMT(i, j, k);

	for (i=0; i < numMarkers; i++)
		is >> markers[i];
	for (i=0; i < numTransforms; i++) {
		is >> j;
		switch (j) {
		case KT_TRANSLATION:
			transforms[i] = new KTranslation();
			is >> (*transforms[i]);
			break;
		case KT_QUATROTATION:
			transforms[i] = new KQuatRotation();
			is >> (*transforms[i]);
			break;
		case KT_EULERROTATION:
			transforms[i] = new KEulerRotation();
			is >> (*transforms[i]);
			break;
		case KT_POLARAXISROTATION:
			transforms[i] = new KPolarAxisRotation();
			is >> (*transforms[i]);
			break;
		}
		transforms[i]->model = this;
	}
}

void KinematicModel::assignDofIndices() {
	int i;

	if (!dofs)
		dofs = new KDofSet();

	dofs->globalSize = 0;
	dofs->perFrameSize = 0;

	for (i=0; i < numTransforms; i++) {
		if (transforms[i]->isCopy)
			continue;

		transforms[i]->dofGIndex = dofs->globalSize;
		dofs->globalSize += transforms[i]->numGDofs();
	}
	dofs->nonMarkerSize = dofs->globalSize;

	for (i=0; i < numMarkers; i++) {
		markers[i].dofIndex = dofs->globalSize;
		dofs->globalSize += 3;
	}

	for (i=0; i < numTransforms; i++) {
		if (transforms[i]->isCopy)
			continue;

		transforms[i]->dofFIndex = dofs->globalSize + dofs->perFrameSize;
		dofs->perFrameSize += transforms[i]->numFDofs();
	}
	
	dofs->numFrames = numMeshes;
	dofs->numDofs = dofs->globalSize + dofs->perFrameSize * numMeshes;
	dofs->data.resize(dofs->numDofs);
	//dofs->data = new double[dofs->numDofs];
}

void KinematicModel::loadDofs(int frame) {
	int i;

	transforms[0]->loadAllDofs(dofs, frame);
	for (i=0; i < numMarkers; i++)
		markers[i].loadDofs(dofs);
}

void KinematicModel::unloadDofs(int frame) {
	int i;

	transforms[0]->unloadAllDofs(dofs, frame);
	for (i=0; i < numMarkers; i++)
		markers[i].unloadDofs(dofs);
}

// --- frame stuff ---

void KinematicModel::initCoords() {
	flushCoords();

	coordData = new KCoord[numTransforms * numMeshes];
	int i;
	for (i=0; i < numTransforms * numMeshes; i++)
		coordData[i].initDeriv(dofs->numDofs);
}

void KinematicModel::flushCoords() {
	if (coordData) {
		delete []coordData;
		coordData = NULL;
	}
}

Mat4d &KinematicModel::trans(int transform, int frame) {
	static Mat4d errMat;
	if (transform < 0 || transform >= numTransforms) {
		cerr << "illegal transform requested: " << transform << endl;
		return errMat;
	}
	if (frame < 0 || frame >= numMeshes) {
		cerr << "illegal frame requested: " << frame << endl;
		return errMat;
	}

	return coordData[frame * numTransforms + transform].mat;
}

Mat4d* &KinematicModel::deriv(int transform, int frame, int dof) {
	static Mat4d *errMat = NULL;
	if (transform < 0 || transform >= numTransforms) {
		cerr << "illegal transform requested: " << transform << endl;
		return errMat;
	}
	if (frame < 0 || frame >= numMeshes) {
		cerr << "illegal frame requested: " << frame << endl;
		return errMat;
	}

	return coordData[frame * numTransforms + transform].deriv[dof];
}

void KinematicModel::setDeriv(int transform, int frame, int dof, Mat4d &m) {
	if (transform < 0 || transform >= numTransforms) {
		cerr << "illegal transform requested: " << transform << endl;
		return;
	}
	if (frame < 0 || frame >= numMeshes) {
		cerr << "illegal frame requested: " << frame << endl;
		return;
	}
	if (dof < 0 || frame >= dofs->numDofs) {
		cerr << "illegal dof requested: " << dof << endl;
		return;
	}

	coordData[frame * numTransforms + transform].setDeriv(dof, m);
}

void KinematicModel::renderMarkers(bool showDeltas, int frame) {
	int i;

	if (!markerData)
		return;

	for (i=0; i < numMarkers; i++) {
		if (markerData->v(frame, markers[i].markerIndex).iszero())
			continue;

		if (showDeltas) {
			glDisable(GL_LIGHTING);
			glColor3f(1, 1, 1);
			glLineWidth(3);
			glBegin(GL_LINES);
				curMarkerPos.v(frame, i).glVertex();
				markerData->v(frame, markers[i].markerIndex).glVertex();
			glEnd();
		}

		glEnable(GL_LIGHTING);
		glColor3f(0.5, 0.5, 0.5);
		glbSphere(markerData->v(frame, markers[i].markerIndex), 0.005);
	}
}

void KinematicModel::loadAll(char *lPrefix) {
	char fname[256];
	int i;

	strcpy(prefix, lPrefix);

	sprintf(fname, "%smarkers.txt", prefix);
	cout << "loading markers from " << fname << endl;
	ifstream fMarkers(fname);
	if (fMarkers.good()) {
		markerData = new KMarkerData();
		fMarkers >> (*markerData);
		fMarkers.close();
	}
	else {
		cout << "can't open " << fname << endl;
	}
	
	sprintf(fname, "%sraw-markers.txt", prefix);
	cout << "loading raw markers from " << fname << endl;
	ifstream rfMarkers(fname);
	if (rfMarkers.good()) {
		rawMarkerData = new KMarkerData();
		rawMarkerData->variableSize = true;
		rawMarkerData->colorData = true;
		rfMarkers >> (*rawMarkerData);
		rfMarkers.close();
	}
	else {
		cout << "can't open " << fname << endl;
	}

// assume kinematics already loaded...
//	initMT(data.kin->markerData->numFrames, data.kin->markerData->numMarkers, data.kin->
//	sprintf(fname, "%s-kin.txt");

	sprintf(fname, "%smeshes.txt", prefix);
	cout << "loading meshes from " << fname << endl;
	ifstream fMeshes(fname);
	if (!fMeshes.good()) {
		cout << "can't open " << fname << endl;
	}
	else {
		fMeshes >> i;
		for (i=0; i < numMeshes; i++)
			fMeshes >> meshNames[i];
	}

	sprintf(fname, "%sdofs.txt", prefix);
	loadDofs(fname);
}

bool KinematicModel::loadKin(char *fname) {
	ifstream in(fname);

	if(!in.good())
		return false;

	getKin(in);
	in.close();

	assignDofIndices();
	loadAll("../data/arm0");
	initCoords();
	calcTransforms();

	return true;
}

bool KinematicModel::saveKin(char *fname) {
	ofstream out(fname);

	if(!out.good())
		return false;

	putKin(out);
	out.close();

	return true;
}

bool KinematicModel::loadDofs(char *fname) {
	ifstream in(fname);

	if (!in.good())
		return false;

	if (dofs)
		delete dofs;

	dofs = new KDofSet;
	loadDofs(in);
	in.close();

	return true;
}

bool KinematicModel::saveDofs(char *fname) {
	ofstream out(fname);

	if (!out.good())
		return false;

	saveDofs(out);
	out.close();

	return true;
}

bool KinematicModel::loadDofs(istream &in) {
	in >> (*dofs);
	return in.good();
}

bool KinematicModel::saveDofs(ostream &out) {
	out << (*dofs);
	return out.good();
}

bool KinematicModel::loadPose(char *fname) {
	int i;
	ifstream in(fname);

	if (!in.good())
		return false;

	char tempStr[256];

	in >> tempStr;
	if (strcmp(tempStr, "pose") != 0) {
		cerr << "warning while loading pose `" << fname << 
			"': expected `pose'; read `" << tempStr << "'." << endl;
	}

	while (1) {
		in >> tempStr;
		if (!in.good() || strcmp(tempStr, "end") == 0)
			break;

		for (i=0; i < numTransforms; i++) {
			if (strcmp(transforms[i]->name, tempStr) == 0)
				break;
		}

		if (i == numTransforms) {
			cerr << "warning while loading pose `" << fname << 
				"': unknown transform `" << tempStr << "'." << endl;

			while (in.good() && strcmp(tempStr, "endp") != 0)
				in >> tempStr;
		}
		else {
			transforms[i]->getVal(in);
			in >> tempStr;
			while (in.good() && strcmp(tempStr, "endp") != 0)
				in >> tempStr;
		}
	}

	in.close();

	return true;
}

bool KinematicModel::savePose(char *fname) {
	ofstream out(fname);

	if (!out.good())
		return false;

	out << "pose" << endl;
	int i;
	for (i=0; i < numTransforms; i++) {
		out << transforms[i]->name << " ";
		transforms[i]->putVal(out);
		out << " endp" << endl;
	}
	out << "end" << endl;
	
	out.close();

	return true;
}

void KinematicModel::fitMarkers(int frame) {
	int i;
	Vec3d curV;
	Mat3d rotM;
	Mat4d curTrans;

	calcTransforms();

	for (i=0; i < numMarkers; i++) {
		// skip missing markers
		if (markerData->v(frame, markers[i].markerIndex).iszero())
			continue;

		curTrans = trans(markers[i].transform, frame);
		curV = vec4to3(curTrans * Vec4d(0, 0, 0, 1));
		curV = markerData->v(frame, markers[i].markerIndex) - curV;
		
		int x, y;
		for (x=0; x < 3; x++)
			for (y=0; y < 3; y++)
				rotM[y][x] = curTrans[x][y];
		curV = rotM * curV;

		markers[i].curVal = curV;
		markers[i].unloadDofs(dofs);
	}

	calcTransforms();
}

void KinematicModel::copyVals(KinematicModel *k) {
	int i;
	for (i=0; i < min(numTransforms, k->numTransforms); i++)
		transforms[i]->copyVal(k->transforms[i]);
}

void KinematicModel::normalizeTransforms() {
	int i, j;
	for (j=0; j < numMeshes; j++) {
		transforms[0]->loadAllDofs(dofs, j);
		for (i=0; i < numTransforms; i++) {
			if (j == 0 || transforms[i]->perFrame)
				transforms[i]->normalize();
		}
		transforms[0]->unloadAllDofs(dofs, j);
	}
}
#endif