#include <fstream>
using namespace std;
#include "mocap.h"
#include "skeleton.h"
#include "ba.h"

const char *MAPPING_NAMES[NUM_MAPPINGS] = {
	"baseTQ<a-X>",
	"baseTQ<t-X>",
	"waistQ<a-X>",
	"abdomenQ<a-X>",
	"lClavicleQ<a-X>",
	"lShoulderQ<a-X>",
	"lElbowA<a-X>",
	"lForearmA<a-Y>",
	"lWristA<a-Z>",
	"rClavicalQ<a-X>",
	"rShoulderQ<a-X>",
	"rElbowA<a-X>",
	"rForearmA<a-Y>",
	"rWristA<a-Z>",
	"lHipQ<a-X>",
	"lKneeA<a-Y>",
	"lShinA<a-Z>",
	"lAnkelA<a-Y>",
	"lFootA<a-Y>",
	"rHipQ<a-X>",
	"rKnee<a-Y>",
	"rShin<a-Z>",
	"rAnkelA<a-Y>",
	"rFootA<a-Y>",
	"Neck<a-X>"
};

MocapData::MocapData() {
	data = NULL;
	isText = false;
}

bool MocapData::load(char *fname) {
	int i, j;
	char junk[80];

	if (fname[strlen(fname)-3] == 't') {
		return loadTxt(fname);
	}

	ifstream in(fname);

	cout << "reading " << fname << endl;

	if (!in.good())
		return false;

	in >> junk; // "frames"
	in >> junk; // "="
	in >> frames;
	in >> junk; // "dofs"
	in >> junk; // "="
	in >> frameSize;

	// the last frame always seems to be bad
//	frames--;

//	baseTInd = baseRInd = waistInd = abdomenInd= lClavInd, lShoulderInd, lElbowAInd, lElbowTInd, lWristInd;
//	int rClavInd, rShoulderInd, rElbowAInd, rElbowTInd, rWristInd;

	for (i=0; i < NUM_MAPPINGS; i++)
		mappings[i] = -1;

	for (i=0; i < frameSize; i++) {
		in >> junk; // frame names

		if (strcmp(junk, "end") == 0)
			break;

		for (j=0; j < NUM_MAPPINGS; j++) {
			if (strcmp(junk, MAPPING_NAMES[j]) == 0) {
				mappings[j] = i;
				break;
			}
		}
	}
	frameSize = i;

	if (data)
		delete []data;
	data = new double[frames * frameSize];

	for (j=0; j < frames; j++) {
		for (i=0; i < frameSize; i++) {
			in >> data[j*frameSize + i];
		}
	}

	cout << "loaded " << frames << " frames with " << frameSize << " dofs" << endl;

	for (i=0; i < NUM_MAPPINGS; i++) {
		if (mappings[i] == -1) {
			cout << "warning: " << MAPPING_NAMES[i] << " not found" << endl;
		}
	}

	isText = false;
//	return in.good();
	return true;
}

bool MocapData::loadTxt(char *fname) {
	int i, j;
	char junk[80];

	ifstream in(fname);
	cout << "reading as text " << fname << endl;
	if (!in.good())
		return false;

	in >> frames;

	frameSize = 60;

	if (data)
		delete []data;
	data = new double[frames * frameSize];

	for (j=0; j < frames; j++) {
		for (i=0; i < frameSize; i++) {
			in >> data[j*frameSize + i];
		}
	}

	cout << "loaded " << frames << " frames with " << frameSize << " dofs" << endl;

	isText = true;
	return true;
}

QuatNorm eulerToQN(double z, double x, double y) {
	return QuatNorm(z*DEG_TO_RAD, 0, 0) * QuatNorm(0, -x*DEG_TO_RAD, 0) * QuatNorm(0, 0, y*DEG_TO_RAD);
}

void MocapData::toSkel(int frame, Skeleton *skel) {
	int ofs = frameSize * frame;
	int i;

	Vec3d &baseT = ((SkelTranslation*)skel->transforms.getT("baseT"))->curVal;
	QuatNorm &baseQ = ((SkelQuatRotation*)skel->transforms.getT("baseQ"))->curQuat;
	double &lElbowA = ((SkelEulerRotation*)skel->transforms.getT("lElbowA"))->curAngle;
	double &lForearmA = ((SkelEulerRotation*)skel->transforms.getT("lForearmA"))->curAngle;
	QuatNorm &lShoulderQ = ((SkelQuatRotation*)skel->transforms.getT("lShoulderQ"))->curQuat;

	if (isText) {
//		baseT = 0.01 * Vec3d(data[ofs+0], data[ofs+1], data[ofs+2]);
//		baseQ = QuatNorm(data[ofs+3]*DEG_TO_RAD, data[ofs+4]*DEG_TO_RAD, data[ofs+5]*DEG_TO_RAD);
		lShoulderQ = eulerToQN(data[ofs+33], data[ofs+34], data[ofs+35]);
		lElbowA = -data[ofs+36] * DEG_TO_RAD;
		lForearmA = -data[ofs+37] * DEG_TO_RAD;
	}
	else {
		double xFactor = -1;
		double yFactor = -1;
		double zFactor = -1;

		double clavXFix = 0; //45;
		double clavZFix = 0; //35;
		double rShoulderXFix = 0; //25;
		double rShoulderYFix = 0;
		double shoulderZFix = 0; //-15;

		// base trans
		baseT = Vec3d(
			(data[ofs+mappings[1]+0]+45.0)/1000.0, 
			(data[ofs+mappings[1]+1]+54.0)/1000.0, 
			(data[ofs+mappings[1]+2]-1028.0)/1000.0);
		// base rot
		baseQ = QuatNorm(xFactor*data[ofs+mappings[0]+0]*DEG_TO_RAD, yFactor*data[ofs+mappings[0]+1]*DEG_TO_RAD, zFactor*data[ofs+mappings[0]+2]*DEG_TO_RAD);
		// waist
		QuatNorm &waistQ = ((SkelQuatRotation*)skel->transforms.getT("waistQ"))->curQuat;
		waistQ = QuatNorm(xFactor*data[ofs+mappings[2]+0]*DEG_TO_RAD, yFactor*data[ofs+mappings[2]+1]*DEG_TO_RAD, zFactor*data[ofs+mappings[2]+2]*DEG_TO_RAD);
		// abdomen
		QuatNorm &abdomenQ = ((SkelQuatRotation*)skel->transforms.getT("abdomenQ"))->curQuat;
		abdomenQ = QuatNorm(xFactor*data[ofs+mappings[3]+0]*DEG_TO_RAD, yFactor*data[ofs+mappings[3]+1]*DEG_TO_RAD, zFactor*data[ofs+mappings[3]+2]*DEG_TO_RAD);

		// l clav
		QuatNorm &lClavicleQ = ((SkelQuatRotation*)skel->transforms.getT("lClavicleQ"))->curQuat;
		lClavicleQ = QuatNorm(
			(xFactor*data[ofs+mappings[4]+0]+clavXFix)*DEG_TO_RAD, 
			yFactor*data[ofs+mappings[4]+1]*DEG_TO_RAD, 
			(zFactor*data[ofs+mappings[4]+2]+clavZFix)*DEG_TO_RAD);
		// l shoulder
		lShoulderQ = QuatNorm(
			xFactor*data[ofs+mappings[5]+0]*DEG_TO_RAD, 
			yFactor*data[ofs+mappings[5]+1]*DEG_TO_RAD, 
			(zFactor*data[ofs+mappings[5]+2]+shoulderZFix)*DEG_TO_RAD);
		// l elbow a
		lElbowA = data[ofs+mappings[6]]*DEG_TO_RAD;
		// l elbow t
		lForearmA = data[ofs+mappings[7]]*DEG_TO_RAD;
		// l wrist a
		double &lWristA = ((SkelEulerRotation*)skel->transforms.getT("lWristA"))->curAngle;
		lWristA = data[ofs+mappings[8]]*DEG_TO_RAD;
		
		// r clav
		QuatNorm &rClavicleQ = ((SkelQuatRotation*)skel->transforms.getT("rClavicleQ"))->curQuat;
		rClavicleQ = QuatNorm(
			(xFactor*data[ofs+mappings[9]+0]-clavXFix)*DEG_TO_RAD, 
			yFactor*data[ofs+mappings[9]+1]*DEG_TO_RAD, 
			(zFactor*data[ofs+mappings[9]+2]-clavXFix)*DEG_TO_RAD);
		// r shoulder
		QuatNorm &rShoulderQ = ((SkelQuatRotation*)skel->transforms.getT("rShoulderQ"))->curQuat;
		rShoulderQ = QuatNorm(
			(xFactor*data[ofs+mappings[10]+0]+rShoulderXFix)*DEG_TO_RAD, 
			(yFactor*data[ofs+mappings[10]+1])*DEG_TO_RAD, 
			(zFactor*data[ofs+mappings[10]+2]-shoulderZFix)*DEG_TO_RAD);
		// r elbow a
		double &rElbowA = ((SkelEulerRotation*)skel->transforms.getT("rElbowA"))->curAngle;
		rElbowA = data[ofs+mappings[11]]*DEG_TO_RAD;
		// r elbow t
		double &rForearmA = ((SkelEulerRotation*)skel->transforms.getT("rForearmA"))->curAngle;
		rForearmA = data[ofs+mappings[12]]*DEG_TO_RAD;
		// r wrist a
		double &rWristA = ((SkelEulerRotation*)skel->transforms.getT("rWristA"))->curAngle;
		rWristA = data[ofs+mappings[13]]*DEG_TO_RAD;

		// l hip
		i = 14;
		QuatNorm &lHipQ = ((SkelQuatRotation*)skel->transforms.getT("lHipQ"))->curQuat;
		lHipQ = QuatNorm(xFactor*data[ofs+mappings[i]+0]*DEG_TO_RAD, yFactor*data[ofs+mappings[i]+1]*DEG_TO_RAD, zFactor*data[ofs+mappings[i]+2]*DEG_TO_RAD);
		// l knee
		i++;
		double &lKneeA = ((SkelEulerRotation*)skel->transforms.getT("lKneeA"))->curAngle;
		lKneeA = data[ofs+mappings[i]]*DEG_TO_RAD;
		// l ankle a
		i++;
		double &lShinA = ((SkelEulerRotation*)skel->transforms.getT("lShinA"))->curAngle;
		lShinA = data[ofs+mappings[i]]*DEG_TO_RAD;
		// l ankle t
		i++;
		double &lAnkleA = ((SkelEulerRotation*)skel->transforms.getT("lAnkleA"))->curAngle;
		lAnkleA = data[ofs+mappings[i]]*DEG_TO_RAD;
		// l foot t
		i++;
		double &lFootA = ((SkelEulerRotation*)skel->transforms.getT("lFootA"))->curAngle;
		lAnkleA = data[ofs+mappings[i]]*DEG_TO_RAD;

		// r hip
		i++;
		QuatNorm &rHipQ = ((SkelQuatRotation*)skel->transforms.getT("rHipQ"))->curQuat;
		rHipQ = QuatNorm(xFactor*data[ofs+mappings[i]+0]*DEG_TO_RAD, yFactor*data[ofs+mappings[i]+1]*DEG_TO_RAD, zFactor*data[ofs+mappings[i]+2]*DEG_TO_RAD);
		// r knee
		i++;
		double &rKneeA = ((SkelEulerRotation*)skel->transforms.getT("rKneeA"))->curAngle;
		rKneeA = data[ofs+mappings[i]]*DEG_TO_RAD;
		// r ankle a
		i++;
		double &rShinA = ((SkelEulerRotation*)skel->transforms.getT("rShinA"))->curAngle;
		rShinA = data[ofs+mappings[i]]*DEG_TO_RAD;
		// r ankle t
		i++;
		double &rAnkleA = ((SkelEulerRotation*)skel->transforms.getT("rAnkleA"))->curAngle;
		rAnkleA = data[ofs+mappings[i]]*DEG_TO_RAD;
		// r foot t
		i++;
		double &rFootA = ((SkelEulerRotation*)skel->transforms.getT("rFootA"))->curAngle;
		rAnkleA = data[ofs+mappings[i]]*DEG_TO_RAD;

		// neck
		i++;
		QuatNorm &neckQ = ((SkelQuatRotation*)skel->transforms.getT("neckR"))->curQuat;
		neckQ = QuatNorm(xFactor*data[ofs+mappings[i]+0]*DEG_TO_RAD, yFactor*data[ofs+mappings[i]+1]*DEG_TO_RAD, zFactor*data[ofs+mappings[i]+2]*DEG_TO_RAD);
	}
}
