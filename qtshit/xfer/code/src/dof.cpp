#include "doppel2.h"
#include "dof.h"
#include "FL/Fl.H"
#include "ba.h"

const int DofSet::LOCAL_TRANS = 0;
const int DofSet::GLOBAL_TRANS = 1;
const int DofSet::MIXED_TRANS = 2;

// SkelTransformDof ---------------------------------------

SkelTransformDof::SkelTransformDof() {
	trans = NULL;
}

void SkelTransformDof::setTransform(SkelTransform *t, int ind) {
	trans = t;
	index = ind;
}

//SkelTransform *SkelTransformDof::getTransform() {
//	return trans;
//}

void SkelTransformDof::initSkel() {
//	isSymmetric = (strcmp(trans->className, "SkelSymmetricTranslation") == 0);
}

bool SkelTransformDof::hasDependency(SkelTransform *t) {
	if (t == trans)
		return true;
	if ((strcmp(t->className, "SkelSymmetricTranslation") == 0) && 
		(((SkelSymmetricTranslation*)t)->orig == trans))
		return true;
	return false;
}

void SkelTransformDof::skelToDof(double &value) {
	value = trans->getDofAddr(index);
}

void SkelTransformDof::dofToSkel(double value) {
	trans->getDofAddr(index) = value;
}

Mat4d SkelTransformDof::deriv(SkelTransform *t) {
	if (t == trans)
		return *trans->curCoord.deriv[index];
	else if ((strcmp(t->className, "SkelSymmetricTranslation") == 0) && 
		(((SkelSymmetricTranslation*)t)->orig == trans))
		return ((SkelSymmetricTranslation*)t)->axisMult[index] * *trans->curCoord.deriv[index];
	else
		return Mat4d() * 0;
}

Vec4d SkelTransformDof::transDeriv(SkelTransform *t) {
	Vec4d ret;
	
	switch (index) {
	case 0:
		ret = Vec4d(1, 0, 0, 0);
		break;
	case 1:
		ret = Vec4d(0, 1, 0, 0);
		break;
	default:
		return Vec4d(0, 0, 1, 0);
	}

	if (t == trans)
		return ret;
	else if ((strcmp(t->className, "SkelSymmetricTranslation") == 0) && 
		(((SkelSymmetricTranslation*)t)->orig == trans))
		return ((SkelSymmetricTranslation*)t)->axisMult[index] * ret;
	else
		return Vec4d(0, 0, 0, 0);
}

// SkelTelescopingDof -------------------------------------

void SkelTelescopingDof::initSkel() {
	if (strcmp(trans->className, "SkelTranslation") != 0) {
		cout << "error: Telescoping DOF added for non-translation!!" << endl;
		return;
	}

	tr = (SkelTranslation*)trans;
	direction = tr->curVal;
	direction.normalize();
}

void SkelTelescopingDof::skelToDof(double &value) {
	direction = tr->curVal;
	value = direction.length();
	direction.normalize();
}

void SkelTelescopingDof::dofToSkel(double value) {
	tr->curVal = direction * value;
}

Mat4d SkelTelescopingDof::deriv(SkelTransform *t) {
	Mat4d ret;
	ret *= 0;
	ret[0][3] = direction[0];
	ret[1][3] = direction[1];
	ret[2][3] = direction[2];

	if (t == trans)
		return ret;
	else if ((strcmp(t->className, "SkelSymmetricTranslation") == 0) && 
		(((SkelSymmetricTranslation*)t)->orig == trans)) {
		ret[0][3] *= ((SkelSymmetricTranslation*)t)->axisMult[0];
		ret[1][3] *= ((SkelSymmetricTranslation*)t)->axisMult[1];
		ret[2][3] *= ((SkelSymmetricTranslation*)t)->axisMult[2];
		return ret;
	}
	else
		return Mat4d() * 0;
}

Vec4d SkelTelescopingDof::transDeriv(SkelTransform *t) {
	if (t == trans)
		return Vec4d(direction[0], direction[1], direction[2], 0);
	else if ((strcmp(t->className, "SkelSymmetricTranslation") == 0) && 
		(((SkelSymmetricTranslation*)t)->orig == trans)) {
		Vec3d axis = ((SkelSymmetricTranslation*)t)->axisMult;
		return Vec4d(axis[0] * direction[0], axis[1] * direction[1], axis[2] * direction[2], 0);
	}
	else
		return Vec4d(0, 0, 0, 0);
}

// DofSet -------------------------------------------------

DofSet::DofSet() {
	className = "DofSet";
}

void DofSet::addTrans(char *name, int kind) {
	SkelTransformDof *dof;

	SkelTransform *trans = skel->transforms.getT(name);
	if (!trans) {
		cout << "error adding transform DOF: unknown transform '" << name << "'" << endl;
	}

	int i;
	for (i=0; i < trans->numDofs(); i++) {
		dof = new SkelTransformDof();
		dof->setTransform(trans, i);

		if (kind == LOCAL_TRANS) {
			localTransDofs.push_back(dof);
		}
		else {
			globalTransDofs.push_back(dof);
		}
		// should support mixed transformations, too...
	}

	capRate = 5;
}

void DofSet::setSkel(Skeleton *sk) {
	origSkel = sk;
	skel = new Skeleton();
	skel->copyFrom(sk);
}

void DofSet::init() {
	numFrames = (int)trueMarkers.size();

	numGlobalDofs = (int)globalTransDofs.size();
	numLocalDofs = (int)localTransDofs.size();
	numDofs = numGlobalDofs + numLocalDofs * numFrames;
	numTransDofs = numGlobalDofs + numLocalDofs;

	skel->allocDerivs(numTransDofs);

	stepCount = 0;
	
	variables.resize(numDofs);
	initSkel();
}

void DofSet::initSkel() {
	int i;
	for (i=0; i < numGlobalDofs; i++)
		globalTransDofs[i]->initSkel();
	for (i=0; i < numLocalDofs; i++)
		localTransDofs[i]->initSkel();
}

void DofSet::skelValsFromOrig(int frame) {
	skel->copyVals(origSkel);

	if (frame == -1) {
		int i;
		for (i=0; i < numFrames; i++)
			skelToDofs(variables, i);
	}
	else 
		skelToDofs(variables, frame);
	initSkel();
}

void DofSet::skelValsToOrig(int frame) {
	dofsToSkel(variables, frame);
	origSkel->copyVals(skel);
	origSkel->updateCoords();
}

void DofSet::interp(NameTableX<double> *weights) {
	int i, j;
	double d;

	for (i=0; i < numLocalDofs; i++) {
		d = 0;

		for (j=0; j < numFrames; j++) {
			d += weights->getX(j) * variables[numGlobalDofs + j*numLocalDofs + i];
		}

		localTransDofs[i]->dofToSkel(d);
	}

	origSkel->copyVals(skel);
	origSkel->updateCoords();
}

void DofSet::loadLandmarks(char *fname) {
	ifstream lmData;
	if (!openIFStream(&lmData, fname, "landmark data"))
		return;
	
	int i, j, count;
	char s[256];
	Marker tempMarker;

	lmData >> count;
	calculatedMarkers->init(count);
	for (i=0; i < count; i++) {
		lmData >> j;
		calcToTrueMap.push_back(j);

		lmData >> s;
		calculatedMarkers->markers[i].trans = skel->transforms.getT(s);

		lmData >> calculatedMarkers->markers[i].boneOfs;
		lmData >> calculatedMarkers->markers[i].pos;

		// for now, we'll eliminate the bone offsets
		//calculatedMarkers->markers[i].pos = calculatedMarkers->markers[i].offset() - ((SkelVec3d*)calculatedMarkers->markers[i].trans)->curVal;
		//calculatedMarkers->markers[i].boneOfs = 1e6;
	}
	lmData.close();
}

void DofSet::skelToDofs(Vecd &values, int frame) {
	// copy from current skeleton values to DOFs
	int i;
	for (i=0; i < numGlobalDofs; i++)
		globalTransDofs[i]->skelToDof(values[i]);
	for (i=0; i < numLocalDofs; i++)
		localTransDofs[i]->skelToDof(values[numGlobalDofs + frame * numLocalDofs + i]);
}

void DofSet::dofsToSkel(Vecd &values, int frame) {
	// copy from DOFs to current skeleton values
	int i;
	for (i=0; i < numGlobalDofs; i++)
		globalTransDofs[i]->dofToSkel(values[i]);
	for (i=0; i < numLocalDofs; i++)
		localTransDofs[i]->dofToSkel(values[numGlobalDofs + frame * numLocalDofs + i]);
	skel->updateCoords();
}

void DofSet::updateTransDerivs() {
	int trans, dof;

	for (dof = 0; dof < numGlobalDofs + numLocalDofs; dof++) {
		SkelTransformDof *curDof;
		if (dof < numGlobalDofs)
			curDof = globalTransDofs[dof];
		else
			curDof = localTransDofs[dof - numGlobalDofs];

		for (trans=0; trans < skel->transforms.size(); trans++) {
			SkelTransform *curTrans = skel->transforms.getT(trans);
			bool isSST = (strcmp(curTrans->className, "SkelSymmetricTranslation") == 0);

			if (curDof->hasDependency(curTrans)) {
				if (!curTrans->parentPtr)
					curTrans->globalDerivs[dof] = curDof->deriv(curTrans);
				else
					curTrans->globalDerivs[dof] = curTrans->parentPtr->globalCoord.mat * curDof->deriv(curTrans);
			}
			else {
				if (!curTrans->parentPtr)
					curTrans->globalDerivs[dof] *= 0;
				else
					curTrans->globalDerivs[dof] = 
						curTrans->parentPtr->globalDerivs[dof] * curTrans->curCoord.mat;
			}
		}
	}
}

void DofSet::evaluateGradient(Vecd &values, Vecd &gradient) {
	int frame, marker, dof, index;
	Vec3d deriv;

	gradient.zeroElements();

	for (frame = 0; frame < numFrames; frame++) {
		dofsToSkel(values, frame);
		skel->updateDerivs();
		updateTransDerivs();

		for (marker = 0; marker < calculatedMarkers->numMarkers; marker++) {
//			if (calculatedMarkers->v(marker).iszero())
//				continue;

			index = marker;
			if (calcToTrueMap.size() > 0)
				index = calcToTrueMap[marker];

			if (trueMarkers[frame]->curPos(index).iszero())
				continue;

			Vec3d dist = calculatedMarkers->curPos(marker) - trueMarkers[frame]->curPos(index);

			for (dof = 0; dof < numTransDofs; dof++) {
				SkelTransformDof *curDof;
				int gradientIndex;
				if (dof < numGlobalDofs) {
					curDof = globalTransDofs[dof];
					gradientIndex = dof;
				}
				else {
					curDof = localTransDofs[dof - numGlobalDofs];
					gradientIndex = dof + (numLocalDofs*frame);
				}

				SkelTransform *curTrans = calculatedMarkers->markers[marker].trans;

				if (calculatedMarkers->markers[marker].boneOfs == 1e6)
					deriv = curTrans->globalDerivs[dof] * calculatedMarkers->offset(marker);
				else {
					if (curDof->hasDependency(curTrans)) {
						Vec4d temp = curTrans->parentPtr->globalCoord.mat * (calculatedMarkers->markers[marker].boneOfs * curDof->transDeriv(curTrans));
						deriv = vec4to3(temp);
					}
					else
						deriv = curTrans->parentPtr->globalDerivs[dof] * calculatedMarkers->offset(marker);
				}
				gradient[gradientIndex] += 2.0 * dist * deriv;
			}
		}
	}
}

double DofSet::evaluateFunction(Vecd &values) {
	int frame, marker, index;
	double ret = 0;
	Vec3d deriv;

	for (frame = 0; frame < numFrames; frame++) {
		dofsToSkel(values, frame);

		for (marker = 0; marker < calculatedMarkers->numMarkers; marker++) {
//			if (calculatedMarkers->v(marker).iszero())
//				continue;

			index = marker;
			if (calcToTrueMap.size() > 0)
				index = calcToTrueMap[marker];

			if (trueMarkers[frame]->curPos(index).iszero())
				continue;

			Vec3d dist = calculatedMarkers->curPos(marker) - trueMarkers[frame]->curPos(index);
			ret += dist.length2();
		}
	}

	lastErr = ret;
	return ret;
}

void DofSet::solverStep() {
/*	if (capRate > 0 && (stepCount % capRate == 0)) {
		cout << "frame " << stepCount << "; error " << lastErr << endl;
		redrawV();
		uiWait();
	}
*/
	redrawVNow();
	char fname[80];
	sprintf(fname, "solve%04d.tga", stepCount); 
	uiScreenshot(fname);
	uiWait();

	stepCount++;
}

void DofSet::drawGL(int f) {
	int marker, index;
	int frame, minFrame, maxFrame;

	if (f < 0) {
		minFrame = 0;
		maxFrame = numFrames;
	}
	else {
		minFrame = f;
		maxFrame = f+1;
	}

	for (frame=minFrame; frame < maxFrame; frame++) {
		dofsToSkel(variables, frame);
		skel->drawGL();

		for (marker = 0; marker < calculatedMarkers->numMarkers; marker++) {
//			if (calculatedMarkers->v(marker).iszero())
//				continue;

			index = marker;
			if (calcToTrueMap.size() > 0)
				index = calcToTrueMap[marker];

			if (trueMarkers[frame]->curPos(index).iszero())
					continue;

			glColor3f(1, 1, 0);
			glbSphere(calculatedMarkers->curPos(marker), 0.01);

			glDisable(GL_LIGHTING);
			glColor3f(0, 0, 0);
			glLineWidth(2.0);
			glBegin(GL_LINES);
			calculatedMarkers->curPos(marker).glVertex();
			trueMarkers[frame]->curPos(index).glVertex();
			glEnd();
			glEnable(GL_LIGHTING);
		}
	}
}

void DofSet::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();

/*	prop = new SLProperty();
	prop->offset = offsetof(Skeleton, transforms);
	prop->type = SL_NAMETABLE_T;
	classInfo->members.addT("transforms", prop);
*/
	classInfo->size = sizeof(DofSet);
	classInfo->newInstance = DofSet::newInstance;
	sl->classes.addT("DofSet", classInfo);
}

SLInterface *DofSet::newInstance(int count) {
	return new DofSet[count];
}
