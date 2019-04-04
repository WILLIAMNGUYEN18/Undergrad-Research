#include "ba.h"
#include "kinematics.h"
#include "markers.h"
#include "skeleton.h"
#include "trimesh.h"

// Marker =================================================

Marker::Marker() {
	className = "Marker";
	name[0] = 0;

	kind = mkrPOS;
	color = Vec3d(1, 0, 0.5);
	skelFrame = -1;
	trans = NULL;
	mesh = NULL;

	boneOfs = 1e6;
}

Vec3d Marker::curPos() {
	if (kind == mkrPT && mesh) {
		pos = mesh->getPt(baryVerts[0]);
		return pos;
	}
	else if (kind == mkrBARY && mesh) {
		pos = baryPos[0] * mesh->getPt(baryVerts[0]);
		if (baryVerts[1] >= 0)
			pos += baryPos[1] * mesh->getPt(baryVerts[1]);
		if (baryVerts[2] >= 0)
			pos += baryPos[2] * mesh->getPt(baryVerts[2]);
		return pos;
	}

	if (trans) {
		if (trans->parentPtr && boneOfs != 1e6)
			return trans->parentPtr->globalCoord.mat * offset();
		else
			return trans->globalCoord.mat * offset();
	}
	return offset();
}

Vec3d Marker::offset() {
	if (boneOfs == 1e6 || !trans)
		return pos;

	return boneOfs * ((SkelVec3d*)trans)->curVal + pos;
}

void Marker::drawGL() {
	color.glColor();
	if (!curPos().iszero())
	glbSphere(curPos(), 0.01);
}

void Marker::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();

	prop = new SLProperty();
	prop->offset = offsetof(Marker, pos);
	prop->type = SL_VEC3D;
	classInfo->members.addT("pos", prop);

	prop = new SLProperty();
	prop->offset = offsetof(Marker, color);
	prop->type = SL_VEC3D;
	classInfo->members.addT("color", prop);

	prop = new SLProperty();
	prop->offset = offsetof(Marker, skelFrame);
	prop->type = SL_INT;
	classInfo->members.addT("skelFrame", prop);

	classInfo->size = sizeof(Marker);
	classInfo->newInstance = Marker::newInstance;
	sl->classes.addT("Marker", classInfo);
}

SLInterface *Marker::newInstance(int count) {
	return new Marker[count];
}

// MarkerSet ==============================================

MarkerSet::MarkerSet() {
	className = "MarkerSet";

	markers = NULL;
	numMarkers = 0;
}

void MarkerSet::init(int nMarkers) {
	numMarkers = nMarkers;
	markers = new Marker[numMarkers];
}

Vec3d &MarkerSet::v(int mIndex) {
	return markers[mIndex].pos;
}

Vec3d &MarkerSet::c(int mIndex) {
	return markers[mIndex].color;
}

Vec3d MarkerSet::offset(int mIndex) {
	return markers[mIndex].offset();
}

Vec3d MarkerSet::curPos(int mIndex) {
	return markers[mIndex].curPos();
}

void MarkerSet::drawGL() {
	int i;
	for (i=0; i < numMarkers; i++)
		markers[i].drawGL();
}

void MarkerSet::loadFromLandmarks(const char *fname) {
	ifstream in;

	init(85); //76);

	if (!openIFStream(&in, fname, "landmark file"))
		return;

	char tempStr[1024];
	int j;
	double tempDouble;

	while (!in.eof()) {
		in.getline(tempStr, 1024);
		if (strncmp(tempStr, "AUX ", 4) == 0) {
			break;
		}
	}

	int index;
	markers[0].pos = Vec3d();
	while (!in.eof()) {
		in >> tempStr;
		if (strncmp(tempStr, "END", 3) == 0)
			break;
		index = atoi(tempStr);
		for (j=0; j < 3; j++)
			in >> tempDouble;
		for (j=0; j < 3; j++) {
			in >> tempDouble;
			if (index > 0 && index <= 73)
				markers[index].pos[j] = tempDouble / 1000.0;
			else if (index >= 74)
				markers[index].pos[j] = tempDouble;
		}
		in.getline(tempStr, 1024);
	}

	// "b" files have butt block marker; gotta remove it
	if (fabs(markers[74].pos[2]) > 10) {
		markers[74].pos = Vec3d();
		markers[75].pos = Vec3d();
	}
}

void MarkerSet::loadText(const char *fname) {
	ifstream in;

	if (!openIFStream(&in, fname, "landmark file"))
		return;

	int i, j;
	Vec3d v;

	in >> j;
	init(j);

	for (i=0; i < j; i++) {
		in >> v;
		markers[i].pos = v;
	}
	in.close();
}

void MarkerSet::loadFromMkr(const char *fname) {
	ifstream in;
	int index;
	char s[256];

	if (!openIFStream(&in, fname, "mkr file"))
		return;
	in.getline(s, 255);

	if (s[0] == '1') {
		int n;
		in >> n;
		init(n);
		for (index = 0; index < n; index++) {
			in >> markers[index].name;
			int mode;
			in >> mode;
			if (mode != 0) {
				markers[index].pos = Vec3d();
				in.getline(s, 255);
			}
			else {
				in >> markers[index].pos[0] >> markers[index].pos[1] >> markers[index].pos[2];
			}
		}
		cout << "loaded " << (index-1) << " markers from " << fname << endl;
	}
	else {
		in.seekg(0);

		init(85);
		for (index=0; index < 85; index++)
			markers[index].pos = Vec3d();


		for (index=1; index < 85; index++) {
			in >> markers[index].name;
			if (!in.good() || markers[index].name[0] == 0)
				break;
			in >> markers[index].pos[0] >> markers[index].pos[1] >> markers[index].pos[2];
		}
		cout << "loaded " << (index-1) << " markers from " << fname << endl;
	}
	in.close();
}

void MarkerSet::saveToMkr(const char *fname) {
	ofstream out;
	int index;

	if (!openOFStream(&out, fname, "mkr file"))
		return;

	for (index=1; index < numMarkers; index++) {
		out << index << " " << markers[index].pos[0] << " " << markers[index].pos[1] << " "<< markers[index].pos[2] << endl;
	}
	out.close();
}

void MarkerSet::save(const char *fname) {
	ofstream out;
	int index;

	if (!openOFStream(&out, fname, "mkr file"))
		return;

	out << "1" << endl;	// version
	out << numMarkers << endl;

	for (index=0; index < numMarkers; index++) {
		if (markers[index].name[0] == 0)
			out << index;
		else
			out << markers[index].name;

		out << " " << markers[index].kind << " ";
		if (markers[index].kind == mkrPOS) {
			out << " " << markers[index].pos[0] << " " << markers[index].pos[1] << " "<< markers[index].pos[2];
		}
		else if (markers[index].kind == mkrPT) {
			out << " " << markers[index].baryVerts[0];
		}
		else {
			out << " " << markers[index].baryVerts[0] << " " << markers[index].baryVerts[1] << " "<< markers[index].baryVerts[2];
			out << " " << markers[index].baryPos[0] << " " << markers[index].baryPos[1] << " "<< markers[index].baryPos[2];
		}

		out << endl;
	}
	out.close();
}

void MarkerSet::load(const char *fname) {
	ifstream in;
	int index;

	if (!openIFStream(&in, fname, "mkr file"))
		return;

	in >> index;
	if (index != 1) {
		cout << "incompatible version (" << index << ") reading " << fname << endl;
		return;
	}

	in >> index;
	init(index);

	for (index=0; index < numMarkers; index++) {
		in >> markers[index].name;
		in >> markers[index].kind;

		if (markers[index].kind == mkrPOS) {
			in >> markers[index].pos[0] >> markers[index].pos[1] >> markers[index].pos[2];
		}
		else if (markers[index].kind == mkrPT) {
			in >> markers[index].baryVerts[0];
		}
		else {
			in >> markers[index].baryVerts[0] >> markers[index].baryVerts[1] >> markers[index].baryVerts[2];
			in >> markers[index].baryPos[0] >> markers[index].baryPos[1] >> markers[index].baryPos[2];
		}
	}
	in.close();
}

void MarkerSet::registerClass(SaveLoad *sl) {
	SLClassInfo *classInfo;
	SLProperty *prop;

	classInfo = new SLClassInfo();

	prop = new SLProperty;
	prop->offset = offsetof(MarkerSet, markers);
	prop->offset2 = offsetof(MarkerSet, numMarkers);
	prop->type = SL_ARRAY;
	classInfo->members.addT("markers", prop);

	classInfo->size = sizeof(MarkerSet);
	classInfo->newInstance = MarkerSet::newInstance;
	sl->classes.addT("MarkerSet", classInfo);
}

SLInterface *MarkerSet::newInstance(int count) {
	return new MarkerSet[count];
}