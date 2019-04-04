#include "ba.h"
#include "uMocap.h"
#include "skeleton.h"

///////////////////////////////////////////////////////////////////////
//  C3D file reader and writer
///////////////////////////////////////////////////////////////////////

#define C3D_REC_SIZE   512

typedef struct c3d_head_t {
	unsigned char	prec_start;
	unsigned char	key;
	short	pnt_cnt;
	short	a_channels;
	short	start_frame;
	short	end_frame;
	short	int_gap;
	float	scale;
	short	rec_start;
	short	a_frames;
	float	freq;
	short	stuff[244];	
} c3d_head;

typedef struct c3d_param_t {
	unsigned char	reserved[2];
	unsigned char	pblocks;
	unsigned char	ftype;
	char stuff[C3D_REC_SIZE-4];	
} c3d_param;

typedef struct c3d_frameSI_t {
	short	x, y, z;
	unsigned char	cam_byte;
	unsigned char	residual;
} c3d_frameSI;

typedef struct c3d_frame_t {
	float	x, y, z;
	float	residual;
} c3d_frame;

float ConvertDecToFloat(char bytes[4]);
void ConvertFloatToDec(float f, char* bytes);


///////////////////////////////////////////////////////////////////////
//  C3D file reader and writer (from Karen)
///////////////////////////////////////////////////////////////////////

float ConvertDecToFloat(char bytes[4]) {
    char p[4]; 
    p[0] = bytes[2]; 
    p[1] = bytes[3]; 
    p[2] = bytes[0]; 
    p[3] = bytes[1]; 
    if (p[0] || p[1] || p[2] || p[3]) 
        --p[3];          // adjust exponent 
    return *(float*)p; 
} 

void ConvertFloatToDec(float f, char* bytes) {
    char* p = (char*)&f; 
    bytes[0] = p[2]; 
    bytes[1] = p[3]; 
    bytes[2] = p[0]; 
    bytes[3] = p[1]; 
    if (bytes[0] || bytes[1] || bytes[2] || bytes[3]) 
        ++bytes[1];      // adjust exponent 
}

// UMocap =================================================

UMocap::UMocap() {
	numFrames = 0;
	data = NULL;
}

void UMocap::init(int iNumFrames, int iNumMarkers) {
	numFrames = iNumFrames;
	numMarkers = iNumMarkers;

	data = new Vec3d[numFrames * numMarkers];
}

bool UMocap::loadC3D(char *fname) {
	char buf[C3D_REC_SIZE];
    FILE *file;
	int i, j;
	Vec3d v;
	c3d_head hdr;
	c3d_param param;
	c3d_frameSI frameSI;
	c3d_frame frame;
	bool bDecFmt = true;

	if (!openFile(&file, fname, "rb", "c3d file"))
		return false;
	 
	//get the header
	if (!fread(buf, C3D_REC_SIZE, 1, file))
		return false;
	memcpy(&hdr, buf, sizeof(hdr));

	//get number format
	if (hdr.rec_start > 2) {
		if (!fread(buf, C3D_REC_SIZE, 1, file))
			return false;
		memcpy(&param, buf, sizeof(param));
		if (param.ftype == 84)
			bDecFmt = false;
	}

	//convert if in dec format
	if (bDecFmt) {
		hdr.freq = ConvertDecToFloat((char*)&hdr.freq);
		hdr.scale = ConvertDecToFloat((char*)&hdr.scale);
	}

	int iNumFrames = hdr.end_frame - hdr.start_frame + 1;
	int iNumMarkers = hdr.pnt_cnt;
	double mC3DScale = hdr.scale;

	// initialize
	init(iNumFrames, iNumMarkers);

	float pntScale = (hdr.scale < 0) ? 1 : hdr.scale;

	//eat parameter records
	for (i = 3; i < hdr.rec_start; i++)
		if (!fread(buf, C3D_REC_SIZE, 1, file))
			return false;
	
	// start retrieving data
	int iRecSize;
	if (mC3DScale < 0)
		iRecSize = sizeof(c3d_frame) + ( hdr.a_channels * hdr.a_frames * sizeof(float));
	else
		iRecSize = sizeof(c3d_frameSI) + ( hdr.a_channels * hdr.a_frames * sizeof(short));

	for (i=0; i<numFrames; i++)
		for(j=0; j<numMarkers; j++)	{
			if (!fread(buf, iRecSize, 1, file))
				return false;
			if (mC3DScale < 0) {
				memcpy(&frame, buf, sizeof(frame));
				v[0] = frame.y / 1000.0;
				v[1] = frame.z / 1000.0;
				v[2] = frame.x / 1000.0;
			}
			else {
				memcpy(&frameSI, buf, sizeof(frameSI));
				v[0] = (float)frameSI.y * pntScale / 1000.0;
				v[1] = (float)frameSI.z * pntScale / 1000.0;
				v[2] = (float)frameSI.x * pntScale / 1000.0;
			}
			data[i * numMarkers + j] = v;
		}

	fclose(file);

//	char *pch = strrchr(fileName, '\\');  //clip leading path
//	if (pch) fileName = pch+1;
	return true;
}

// UMocapPoses ============================================

const int NUM_QUATS = 10;
const char QUAT_NAMES[NUM_QUATS][20] = {
	"baseQ",
	"lHipQ", "lAnkleA", 
	"rHipQ", "rAnkleA", 
	"waistQ", "abdomenQ", "neckR",
	"lShoulderCR", "rShoulderCR"};
const int NUM_ROTS = 12;
const char ROT_NAMES[NUM_ROTS][20] = {
	"lKneeA", "rKneeA", 
	"lElbowA", "lForearmA", "lWristA",
	"rElbowA", "rForearmA", "rWristA",
	"lClavicleZR", "lClavicleXR",
	"rClavicleZR", "rClavicleXR"};
const int NUM_TRANS = 1;
const char TRANS_NAMES[NUM_TRANS][20] = {"baseT"};


UMocapPoses::UMocapPoses() {
	numFrames = 0;
	numVars = NUM_QUATS * 4 + NUM_ROTS + 3 * NUM_TRANS;
	data = NULL;
	fixQ = new QuatNorm[NUM_QUATS];
	fixR = new double[NUM_ROTS];
	fixT = new Vec3d[NUM_TRANS];
	memset(fixR, 0, sizeof(double)*NUM_ROTS);
}

void UMocapPoses::init(int iNumFrames, int iNumVars) {
	numFrames = iNumFrames;
	numVars = iNumVars;

	data = new double[numFrames * numVars];
}

/*
void getToken(ifstream &in, char *s) {
	int sInd = 0;
	char c;

	while (in.good() && sInd < 255) {
		in >> c;
		if (c == ',' || (c < ' ' && sInd > 0))
			break;
		if (c < ' ')
			continue;
		s[sInd] = c;
		sInd++;
	}
	s[sInd] = 0;
}

bool UMocapPoses::loadCSV(char *fname) {
	ifstream in;
	int frame, i, j;
	char s[256];

	numFrames = 0;
	if (!openIFStream(&in, fname, "CSV file"))
		return false;

	getToken(in, s);	// subject name
	getToken(in, s);	// number of frames
	numFrames = atoi(s);
	getToken(in, s);	// "Frame"
	numVars = 0;
	while (1) {
		getToken(in, s);
		if (s[0] == '1' || !in.good())
			break;
		numVars++;
	}
	cout << "loading " << numFrames << " frames and " << numVars << " vars" << endl;
	init(numFrames, numVars);

	for (i=0; i < numFrames; i++) {
		for (j=0; j < numVars; j++) {
			getToken(in, s);
			data[(i*numVars) + j] = atof(s);
		}
		getToken(in, s);	// next frame index
	}

	for (i=0; i < 11; i++)
		fixQ[i] = QuatNorm();
	memset(fixR, 0, sizeof(double)*8);

	return true;
}*/

QuatNorm amcEulerToQuat(Vec3d v) {
	v *= DEG_TO_RAD;
//	return QuatNorm(v[2], 0, 0) * QuatNorm(0, 0, -v[1]) * QuatNorm(0, v[0], 0);
	return QuatNorm(0, v[0], 0) * QuatNorm(0, 0, -v[1]) * QuatNorm(v[2], 0, 0);
}

bool UMocapPoses::load(const char *fname) {
	const char *ext = fname + (strlen(fname) - 3);

	if (stricmp(ext, "amc") == 0) 
		return loadAMC(fname);
	else if (stricmp(ext, "csv") == 0)
		return loadCSV(fname);
	else if (stricmp(ext, "dof") == 0)
		return loadDOF(fname);
	cout << "unknown extension" << endl;
	return false;
}

bool UMocapPoses::loadAMC(const char *fname) {
	ifstream in;
	char s[256];
	int i, frame = 0;
	double d;
	Vec3d v;
	QuatNorm q;
	Mat4d m;

	if (!openIFStream(&in, fname, "AMC file"))
		return false;

	in.getline(s, 255);	// #!OML:ASF
	in.getline(s, 255);	// :FULLY-SPECIFIED
	in.getline(s, 255);	// :DEGREES
	
	vector<double> vars;
	int rotOfs = NUM_QUATS * 4;
	int transOfs = rotOfs + NUM_ROTS;
	int iNumVars = transOfs + NUM_TRANS * 3;

	QuatNorm fixRotT = QuatNorm(PI/2, 0, 0) * QuatNorm(0, -PI/2, 0);
	QuatNorm fixRot = QuatNorm(fixRotT.x, fixRotT.y, fixRotT.z, -fixRotT.w);
//	cout << fixRot.toMatrixD() * Vec3d(1, 0, 0) << endl;
//	cout << fixRot.toMatrixD() * Vec3d(0, 1, 0) << endl;
//	cout << fixRot.toMatrixD() * Vec3d(0, 0, 1) << endl;


/* testing
	for (i=0; i < 10; i++) {
		Vec3d v;
		QuatNorm q(boundedRand(-1,1), boundedRand(-1,1), boundedRand(-1,1), boundedRand(-1,1));
		QuatNorm q2;
		q.normalize();
		cout << q << endl;
//		v = -matToEuler(q.toMatrixD());
//		q2 = (QuatNorm(v[0], 0, 0) * QuatNorm(0, v[1], 0) * QuatNorm(0, 0, v[2]));
//		q2.normalize();
//		cout << q2 << endl;

		v = -matToEulerXZY(q.toMatrixD());
		q2 = (QuatNorm(0, v[2], 0) * QuatNorm(0, 0, v[1]) * QuatNorm(v[0], 0, 0));
		q2.normalize();
		cout << q2 << endl;
		cout << endl;
	}*/

	while (1) {
		in >> i;
		if (!in.good() || frame+1 != i)
			break;

		int ofs = (int)vars.size();
		for (i=0; i < iNumVars; i++) {
			if (i < rotOfs && (i%4 == 3))
				vars.push_back(1);
			else
				vars.push_back(0);
		}

		in >> s; // root
		in >> v;
		v *= 0.0254;
		vars[ofs + transOfs + 0] = -v[2];
		vars[ofs + transOfs + 1] = -v[0];
		vars[ofs + transOfs + 2] = -v[1];
		in >> v;
		q = amcEulerToQuat(v);
		vars[ofs + 0] = q[0];
		vars[ofs + 1] = q[1];
		vars[ofs + 2] = q[2];
		vars[ofs + 3] = q[3];

		in >> s; // lowerback
		in >> v;
		q = amcEulerToQuat(v);
		vars[ofs + 5*4 + 0] = q[0];
		vars[ofs + 5*4 + 1] = q[1];
		vars[ofs + 5*4 + 2] = q[2];
		vars[ofs + 5*4 + 3] = q[3];
		in >> s; // upperback
		in >> v;
		q = amcEulerToQuat(v);
		in >> s; // thorax
		in >> v;
		q = amcEulerToQuat(v) * q;
		vars[ofs + 6*4 + 0] = q[0];
		vars[ofs + 6*4 + 1] = q[1];
		vars[ofs + 6*4 + 2] = q[2];
		vars[ofs + 6*4 + 3] = q[3];
		in >> s; // lowerneck
		in >> v;
		q = amcEulerToQuat(v);
		in >> s; // upperneck
		in >> v;
		q = amcEulerToQuat(v) * q;
		vars[ofs + 7*4 + 0] = q[0];
		vars[ofs + 7*4 + 1] = q[1];
		vars[ofs + 7*4 + 2] = q[2];
		vars[ofs + 7*4 + 3] = q[3];
		in >> s; // head
		in >> v;
		in >> s; // rclavicle
		in >> d;
		vars[ofs + rotOfs + 10] = -d * DEG_TO_RAD;
		in >> d;
		vars[ofs + rotOfs + 11] = d * DEG_TO_RAD;
		in >> s; // rhumerus
		in >> v;
		q = 
			QuatNorm(-v[0] * DEG_TO_RAD, 0, 0) * 
			QuatNorm(0, v[1] * DEG_TO_RAD, 0) *
			QuatNorm(0, 0, v[2] * DEG_TO_RAD) *
			QuatNorm(0, 45 * DEG_TO_RAD, 0);
		q.normalize();
		vars[ofs + 9*4 + 0] = q.x;
		vars[ofs + 9*4 + 1] = q.y;
		vars[ofs + 9*4 + 2] = q.z;
		vars[ofs + 9*4 + 3] = q.w;
		in >> s; // rradius
		in >> d;
		vars[ofs + rotOfs + 5] = d * DEG_TO_RAD;
		in >> s; // rwrist
		in >> d;
		vars[ofs + rotOfs + 6] = d * DEG_TO_RAD;
		in >> s; // rhand
		in >> d;
		vars[ofs + rotOfs + 7] = d * DEG_TO_RAD;
		in >> d;
		in >> s; // rfingers
		in >> d;
		in >> s; // rthumb
		in >> d;
		in >> d;
		in >> s; // lclavicle
		in >> d;
		vars[ofs + rotOfs + 8] = -d * DEG_TO_RAD;
		in >> d;
		vars[ofs + rotOfs + 9] = -d * DEG_TO_RAD;
		in >> s; // lhumerus
		in >> v;
		q = 
			QuatNorm(v[0] * DEG_TO_RAD, 0, 0) * 
			QuatNorm(0, -v[1] * DEG_TO_RAD, 0) *
			QuatNorm(0, 0, v[2] * DEG_TO_RAD) *
			QuatNorm(0, 45 * DEG_TO_RAD, 0);
		q.normalize();
		vars[ofs + 8*4 + 0] = q.x;
		vars[ofs + 8*4 + 1] = q.y;
		vars[ofs + 8*4 + 2] = q.z;
		vars[ofs + 8*4 + 3] = q.w;
//		v = matToEulerXZY(q.toMatrixD());
		in >> s; // lradius
		in >> d;
		vars[ofs + rotOfs + 2] = -d * DEG_TO_RAD;
		in >> s; // lwrist
		in >> d;
		vars[ofs + rotOfs + 3] = -d * DEG_TO_RAD;
		in >> s; // lhand
		in >> d;
		vars[ofs + rotOfs + 4] = -d * DEG_TO_RAD;
		in >> d;
		in >> s; // lfingers
		in >> d;
		in >> s; // lthumb
		in >> d;
		in >> d;
		in >> s; // rfemur
		in >> v;
		q = amcEulerToQuat(v);
		vars[ofs + 4*3 + 0] = q[0];
		vars[ofs + 4*3 + 1] = q[1];
		vars[ofs + 4*3 + 2] = q[2];
		vars[ofs + 4*3 + 3] = q[3];
		in >> s; // rtiba
		in >> d;
		d *= DEG_TO_RAD;
		vars[ofs + rotOfs + 1] = -d;
		in >> s; // rfoot
		in >> v[0];
		v[1] = 0;
		in >> v[2];
		q = amcEulerToQuat(v);
//		vars[ofs + 4*4 + 0] = q[0];
//		vars[ofs + 4*4 + 1] = q[1];
//		vars[ofs + 4*4 + 2] = q[2];
//		vars[ofs + 4*4 + 3] = q[3];
		in >> s; // rtoes
		in >> d;
		in >> s; // lfemur
		in >> v;
		q = amcEulerToQuat(v);
		vars[ofs + 4*1 + 0] = q[0];
		vars[ofs + 4*1 + 1] = q[1];
		vars[ofs + 4*1 + 2] = q[2];
		vars[ofs + 4*1 + 3] = q[3];
		in >> s; // ltiba
		in >> d;
		d *= DEG_TO_RAD;
		vars[ofs + rotOfs + 0] = -d;
		in >> s; // lfoot
		in >> v[0];
		v[1] = 0;
		in >> v[2];
		q = amcEulerToQuat(v);
//		vars[ofs + 4*2 + 0] = q[0];
//		vars[ofs + 4*2 + 1] = q[1];
//		vars[ofs + 4*2 + 2] = q[2];
//		vars[ofs + 4*2 + 3] = q[3];
		in >> s; // ltoes
		in >> d;

		frame++;
	}

	cout << "loaded " << frame << " frames" << endl;
	in.close();

	init(frame, iNumVars);
	for (i=0; i < vars.size(); i++)
		data[i] = vars[i];
	return true;
}

double lud(vector<char*> &dofNames, double *dofs, char *name) {
	int i;
	for (i=0; i < (int)dofNames.size(); i++) {
		if (stricmp(dofNames[i], name) == 0)
			return dofs[i];
	}
	cout << "unknown dof " << name << endl;
	return 0;
}

bool UMocapPoses::loadCSV(const char *fname) {
	int numDofs;
	vector <char*> dofNames;
	double *dofs;
	vector <double> vars;
	int rotOfs = NUM_QUATS * 4;
	int transOfs = rotOfs + NUM_ROTS;
	int iNumVars = transOfs + NUM_TRANS * 3;

	ifstream in;
	char s[4096], *curN;
	int i, frame = 0;
	double d;
	Vec3d v;
	QuatNorm q;
	Mat4d m;

	if (!openIFStream(&in, fname, "CSV file"))
		return false;

	in.getline(s, 4095);	// human
	in.getline(s, 4095);	// dof names
	char *temp = strtok(s, ",\n");
	temp = strtok(NULL, ",\n");
	while (temp && strlen(temp) > 0) {
		dofNames.push_back(strdup(temp));
		temp = strtok(NULL, ",\n");
	}
	numDofs = (int)dofNames.size();
	dofs = new double[numDofs];

	while (1) {
		in.getline(s, 4095);	// dof values
		if (!in.good() || strlen(s) < 2)
			break;

		temp = strtok(s, ",\n");
		temp = strtok(NULL, ",\n");
		for (i=0; i < numDofs; i++) {
			if (temp)
				dofs[i] = atof(temp);
			else
				dofs[i] = 0;
			temp = strtok(NULL, ",\n");
		}

		int ofs = (int)vars.size();
		for (i=0; i < iNumVars; i++) {
			if (i < rotOfs && (i%4 == 3))
				vars.push_back(1);
			else
				vars.push_back(0);
		}

		// base
		vars[ofs + transOfs + 0] = -lud(dofNames, dofs, "pelvis<t-X>") / 1000;
		vars[ofs + transOfs + 1] = -lud(dofNames, dofs, "pelvis<t-Y>") / 1000;
		vars[ofs + transOfs + 2] = lud(dofNames, dofs, "pelvis<t-Z>") / 1000;
		v = Vec3d(lud(dofNames, dofs, "pelvis<a-X>"),
			lud(dofNames, dofs, "pelvis<a-Y>"),
			lud(dofNames, dofs, "pelvis<a-Z>"));
		q = 
			QuatNorm(0, 0, -v[2] * DEG_TO_RAD) *
			QuatNorm(0, v[1] * DEG_TO_RAD, 0) *
			QuatNorm(v[0] * DEG_TO_RAD, 0, 0);
		vars[ofs + 0] = q[0];
		vars[ofs + 1] = q[1];
		vars[ofs + 2] = q[2];
		vars[ofs + 3] = q[3];

		// lhip
		v = Vec3d(lud(dofNames, dofs, "lfemur<a-X>"),
			lud(dofNames, dofs, "lfemur<a-Y>"),
			lud(dofNames, dofs, "lfemur<a-Z>")) * 0.5;
		q = 
			QuatNorm(0, 0, -v[2] * DEG_TO_RAD) *
			QuatNorm(0, v[1] * DEG_TO_RAD, 0) *
			QuatNorm(v[0] * DEG_TO_RAD, 0, 0);
		q.normalize();
		vars[ofs + 1*4 + 0] = q.x;
		vars[ofs + 1*4 + 1] = q.y;
		vars[ofs + 1*4 + 2] = q.z;
		vars[ofs + 1*4 + 3] = q.w;
		// lankle
		v = Vec3d(lud(dofNames, dofs, "lfoot<a-X>"),
			lud(dofNames, dofs, "lfoot<a-Y>"),
			lud(dofNames, dofs, "lfoot<a-Z>")) * 0.5;
		q = 
			QuatNorm(0, 0, -v[2] * DEG_TO_RAD) *
			QuatNorm(0, v[1] * DEG_TO_RAD, 0) *
			QuatNorm(v[0] * DEG_TO_RAD, 0, 0);
		q.normalize();
		vars[ofs + 2*4 + 0] = q.x;
		vars[ofs + 2*4 + 1] = q.y;
		vars[ofs + 2*4 + 2] = q.z;
		vars[ofs + 2*4 + 3] = q.w;
		// rhip
		v = Vec3d(lud(dofNames, dofs, "rfemur<a-X>"),
			lud(dofNames, dofs, "rfemur<a-Y>"),
			lud(dofNames, dofs, "rfemur<a-Z>")) * 0.5;
		q = 
			QuatNorm(0, 0, -v[2] * DEG_TO_RAD) *
			QuatNorm(0, v[1] * DEG_TO_RAD, 0) *
			QuatNorm(v[0] * DEG_TO_RAD, 0, 0);
		q.normalize();
		vars[ofs + 3*4 + 0] = q.x;
		vars[ofs + 3*4 + 1] = q.y;
		vars[ofs + 3*4 + 2] = q.z;
		vars[ofs + 3*4 + 3] = q.w;
		// lankle
		v = Vec3d(lud(dofNames, dofs, "rfoot<a-X>"),
			lud(dofNames, dofs, "rfoot<a-Y>"),
			lud(dofNames, dofs, "rfoot<a-Z>")) * 0.5;
		q = 
			QuatNorm(0, 0, -v[2] * DEG_TO_RAD) *
			QuatNorm(0, v[1] * DEG_TO_RAD, 0) *
			QuatNorm(v[0] * DEG_TO_RAD, 0, 0);
		q.normalize();
		vars[ofs + 4*4 + 0] = q.x;
		vars[ofs + 4*4 + 1] = q.y;
		vars[ofs + 4*4 + 2] = q.z;
		vars[ofs + 4*4 + 3] = q.w;
		// waist
		v = Vec3d(lud(dofNames, dofs, "thorax<a-X>"),
			lud(dofNames, dofs, "thorax<a-Y>"),
			lud(dofNames, dofs, "thorax<a-Z>")) * 0.5;
		q = 
			QuatNorm(0, 0, -v[2] * DEG_TO_RAD) *
			QuatNorm(0, v[1] * DEG_TO_RAD, 0) *
			QuatNorm(v[0] * DEG_TO_RAD, 0, 0);
		QuatNorm q2 = slerp(QuatNorm(), q, 0.3);
		q2.normalize();
		vars[ofs + 5*4 + 0] = q2.x;
		vars[ofs + 5*4 + 1] = q2.y;
		vars[ofs + 5*4 + 2] = q2.z;
		vars[ofs + 5*4 + 3] = q2.w;
		// abdomen
		q2 = slerp(QuatNorm(), q, 0.7);
		q2.normalize();
		vars[ofs + 6*4 + 0] = q2.x;
		vars[ofs + 6*4 + 1] = q2.y;
		vars[ofs + 6*4 + 2] = q2.z;
		vars[ofs + 6*4 + 3] = q2.w;
		// neck
		v = Vec3d(lud(dofNames, dofs, "head<a-X>"),
			lud(dofNames, dofs, "head<a-Y>"),
			lud(dofNames, dofs, "head<a-Z>")) * 0.5;
		q = 
			QuatNorm(0, 0, -v[2] * DEG_TO_RAD) *
			QuatNorm(0, v[1] * DEG_TO_RAD, 0) *
			QuatNorm(v[0] * DEG_TO_RAD, 0, 0);
		q.normalize();
		vars[ofs + 7*4 + 0] = q.x;
		vars[ofs + 7*4 + 1] = q.y;
		vars[ofs + 7*4 + 2] = q.z;
		vars[ofs + 7*4 + 3] = q.w;

		// lShoulder
		v = Vec3d(lud(dofNames, dofs, "lhumerus<a-X>"),
			lud(dofNames, dofs, "lhumerus<a-Y>"),
			lud(dofNames, dofs, "lhumerus<a-Z>"));
		q = 
			QuatNorm(0, 90 * DEG_TO_RAD, 0) *
			QuatNorm(0, 0, -v[2] * DEG_TO_RAD) *
			QuatNorm(0, v[1] * DEG_TO_RAD, 0) *
			QuatNorm(v[0] * DEG_TO_RAD, 0, 0);
		q.normalize();
		q = slerp(q, QuatNorm(0.244482, -0.532536, -0.296743, -0.754041), 0.5);
		vars[ofs + 8*4 + 0] = q.x;
		vars[ofs + 8*4 + 1] = q.y;
		vars[ofs + 8*4 + 2] = q.z;
		vars[ofs + 8*4 + 3] = q.w;
		// rShoulder
		v = Vec3d(lud(dofNames, dofs, "rhumerus<a-X>"),
			lud(dofNames, dofs, "rhumerus<a-Y>"),
			lud(dofNames, dofs, "rhumerus<a-Z>"));
		q = 
			QuatNorm(0, 90 * DEG_TO_RAD, 0) *
			QuatNorm(0, 0, -v[2] * DEG_TO_RAD) *
			QuatNorm(0, v[1] * DEG_TO_RAD, 0) *
			QuatNorm(v[0] * DEG_TO_RAD, 0, 0);
		q.normalize();
		q = slerp(q, QuatNorm(-0.269887, -0.512545, 0.311801, -0.753153), 0.5);
		vars[ofs + 9*4 + 0] = q.x;
		vars[ofs + 9*4 + 1] = q.y;
		vars[ofs + 9*4 + 2] = q.z;
		vars[ofs + 9*4 + 3] = q.w;

		// lKneeA
		vars[ofs + rotOfs + 0] = -lud(dofNames, dofs, "ltibia<a-Y>") * DEG_TO_RAD;
		// rKneeA
		vars[ofs + rotOfs + 1] = -lud(dofNames, dofs, "rtibia<a-Y>") * DEG_TO_RAD;
		// lElbowA
		vars[ofs + rotOfs + 2] = lud(dofNames, dofs, "lradius<a-Z>") * DEG_TO_RAD;
		// lForearmA
		vars[ofs + rotOfs + 3] = -lud(dofNames, dofs, "lradius<a-Y>") * DEG_TO_RAD;
		// lWristA
		vars[ofs + rotOfs + 4] = -lud(dofNames, dofs, "lhand<a-X>") * DEG_TO_RAD;
		// rElbowA
		vars[ofs + rotOfs + 5] = lud(dofNames, dofs, "rradius<a-Z>") * DEG_TO_RAD;
		// rForearmA
		vars[ofs + rotOfs + 6] = -lud(dofNames, dofs, "rradius<a-Y>") * DEG_TO_RAD;
		// lWristA
		vars[ofs + rotOfs + 7] = -lud(dofNames, dofs, "rhand<a-X>") * DEG_TO_RAD;

		// lClavicle
		vars[ofs + rotOfs + 8] = lud(dofNames, dofs, "lclavicle<a-Z>") * DEG_TO_RAD;
		vars[ofs + rotOfs + 9] = -lud(dofNames, dofs, "lclavicle<a-X>") * DEG_TO_RAD;
		// rClavicle
		vars[ofs + rotOfs + 10] = lud(dofNames, dofs, "rclavicle<a-Z>") * DEG_TO_RAD;
		vars[ofs + rotOfs + 11] = -lud(dofNames, dofs, "rclavicle<a-X>") * DEG_TO_RAD;

		frame++;
	}

	cout << "loaded " << frame << " frames" << endl;
	in.close();

	init(frame, iNumVars);
	for (i=0; i < vars.size(); i++)
		data[i] = vars[i];

	delete []dofs;
	return true;
}



bool UMocapPoses::loadDOF(const char *fname) {
	int numDofs;
	vector <char*> dofNames;
	double *dofs;
	vector <double> vars;
	int rotOfs = NUM_QUATS * 4;
	int transOfs = rotOfs + NUM_ROTS;
	int iNumVars = transOfs + NUM_TRANS * 3;

	ifstream in;
	char s[4096], *curN;
	int i, frame = 0;
	double d;
	Vec3d v;
	QuatNorm q;
	Mat4d m;

	if (!openIFStream(&in, fname, "DOF file"))
		return false;

	in.getline(s, 4095);	// frames = # dofs = #
	in.getline(s, 4095);	// dof names
	char *temp = strtok(s, " \t\n");
//	temp = strtok(NULL, " \n");
	while (temp && strlen(temp) > 0) {
		dofNames.push_back(strdup(temp));
		temp = strtok(NULL, " \t\n");
	}
	numDofs = (int)dofNames.size();
	cout << "found " << numDofs << " dofs" << endl;
	dofs = new double[numDofs];

	while (1) {
		in.getline(s, 4095);	// dof values
		if (!in.good() || strlen(s) < 2)
			break;

		temp = strtok(s, " \t\n");
//		temp = strtok(NULL, " \t\n");
		for (i=0; i < numDofs; i++) {
			if (temp)
				dofs[i] = atof(temp);
			else
				dofs[i] = 0;
			temp = strtok(NULL, " \t\n");
		}

		int ofs = (int)vars.size();
		for (i=0; i < iNumVars; i++) {
			if (i < rotOfs && (i%4 == 3))
				vars.push_back(1);
			else
				vars.push_back(0);
		}

		// base
		vars[ofs + transOfs + 0] = lud(dofNames, dofs, "pelvis_trans_x") / 1;
		vars[ofs + transOfs + 1] = -lud(dofNames, dofs, "pelvis_trans_z") / 1; //ok
		vars[ofs + transOfs + 2] = lud(dofNames, dofs, "pelvis_trans_y") / 1;
		v = Vec3d(lud(dofNames, dofs, "pelvis_euler_x"),
			lud(dofNames, dofs, "pelvis_euler_y"),
			lud(dofNames, dofs, "pelvis_euler_z"));
		q = 
			QuatNorm(0, v[2], 0) *
			QuatNorm(0, 0, -v[1]) *
			QuatNorm(-v[0], 0, 0);
		vars[ofs + 0] = q[0];
		vars[ofs + 1] = q[1];
		vars[ofs + 2] = q[2];
		vars[ofs + 3] = q[3];

		// lhip
		v = Vec3d(lud(dofNames, dofs, "l_thigh_euler_x"),
			lud(dofNames, dofs, "l_thigh_euler_y"),
			lud(dofNames, dofs, "l_thigh_euler_z"));
		q = 
			QuatNorm(0, v[2], 0) *
			QuatNorm(0, 0, -v[1]) *
			QuatNorm(-v[0], 0, 0);
		q.normalize();
		vars[ofs + 1*4 + 0] = q.x;
		vars[ofs + 1*4 + 1] = q.y;
		vars[ofs + 1*4 + 2] = q.z;
		vars[ofs + 1*4 + 3] = q.w;
		// lankle
		v = Vec3d(lud(dofNames, dofs, "l_ankle_euler_x"),
			lud(dofNames, dofs, "l_ankle_euler_y"),
			lud(dofNames, dofs, "l_ankle_euler_z"));
		q = 
			QuatNorm(0, v[2], 0) *
			QuatNorm(0, 0, -v[1]) *
			QuatNorm(-v[0], 0, 0);
		q.normalize();
		vars[ofs + 2*4 + 0] = q.x;
		vars[ofs + 2*4 + 1] = q.y;
		vars[ofs + 2*4 + 2] = q.z;
		vars[ofs + 2*4 + 3] = q.w;
		// rhip
		v = Vec3d(lud(dofNames, dofs, "r_thigh_euler_x"),
			lud(dofNames, dofs, "r_thigh_euler_y"),
			lud(dofNames, dofs, "r_thigh_euler_z"));
		q = 
			QuatNorm(0, v[2], 0) *
			QuatNorm(0, 0, -v[1]) *
			QuatNorm(-v[0], 0, 0);
		q.normalize();
		vars[ofs + 3*4 + 0] = q.x;
		vars[ofs + 3*4 + 1] = q.y;
		vars[ofs + 3*4 + 2] = q.z;
		vars[ofs + 3*4 + 3] = q.w;
		// lankle
		v = Vec3d(lud(dofNames, dofs, "r_ankle_euler_x"),
			lud(dofNames, dofs, "r_ankle_euler_y"),
			lud(dofNames, dofs, "r_ankle_euler_z"));
		q = 
			QuatNorm(0, v[2], 0) *
			QuatNorm(0, 0, -v[1]) *
			QuatNorm(-v[0], 0, 0);
		q.normalize();
		vars[ofs + 4*4 + 0] = q.x;
		vars[ofs + 4*4 + 1] = q.y;
		vars[ofs + 4*4 + 2] = q.z;
		vars[ofs + 4*4 + 3] = q.w;
		// waist
		v = Vec3d(lud(dofNames, dofs, "abdomen_euler_x"),
			lud(dofNames, dofs, "abdomen_euler_y"),
			lud(dofNames, dofs, "abdomen_euler_z"));
		q = 
			QuatNorm(0, v[2], 0) *
			QuatNorm(0, 0, -v[1]) *
			QuatNorm(-v[0], 0, 0);
		vars[ofs + 5*4 + 0] = q.x;
		vars[ofs + 5*4 + 1] = q.y;
		vars[ofs + 5*4 + 2] = q.z;
		vars[ofs + 5*4 + 3] = q.w;
		// abdomen
		v = Vec3d(lud(dofNames, dofs, "spine_euler_x"),
			lud(dofNames, dofs, "spine_euler_y"),
			lud(dofNames, dofs, "spine_euler_z"));
		q = 
			QuatNorm(0, v[2], 0) *
			QuatNorm(0, 0, -v[1]) *
			QuatNorm(-v[0], 0, 0);
		vars[ofs + 6*4 + 0] = q.x;
		vars[ofs + 6*4 + 1] = q.y;
		vars[ofs + 6*4 + 2] = q.z;
		vars[ofs + 6*4 + 3] = q.w;
		// neck
		v = Vec3d(lud(dofNames, dofs, "neck_euler_x"),
			lud(dofNames, dofs, "neck_euler_y"),
			lud(dofNames, dofs, "neck_euler_z"));
		q = 
			QuatNorm(0, v[2], 0) *
			QuatNorm(0, 0, -v[1]) *
			QuatNorm(-v[0], 0, 0);
		q.normalize();
		vars[ofs + 7*4 + 0] = q.x;
		vars[ofs + 7*4 + 1] = q.y;
		vars[ofs + 7*4 + 2] = q.z;
		vars[ofs + 7*4 + 3] = q.w;

		// lShoulder
		v = Vec3d(lud(dofNames, dofs, "l_bicep_euler_x"),
			lud(dofNames, dofs, "l_bicep_euler_y"),
			lud(dofNames, dofs, "l_bicep_euler_z"));
		q = 
			QuatNorm(-v[2], 0, 0) *
			QuatNorm(0, -v[1], 0) *
			QuatNorm(0, 0, v[0]) *
			QuatNorm(0, 90 * DEG_TO_RAD, 0) *
			QuatNorm(-2.1, 0, 0);
		q.normalize();
//		q = slerp(q, QuatNorm(0.244482, -0.532536, -0.296743, -0.754041), 0.5);
		vars[ofs + 8*4 + 0] = q.x;
		vars[ofs + 8*4 + 1] = q.y;
		vars[ofs + 8*4 + 2] = q.z;
		vars[ofs + 8*4 + 3] = q.w;
		// rShoulder
		v = Vec3d(lud(dofNames, dofs, "r_bicep_euler_x"),
			lud(dofNames, dofs, "r_bicep_euler_y"),
			lud(dofNames, dofs, "r_bicep_euler_z"));
		q = 
			QuatNorm(v[2], 0, 0) *
			QuatNorm(0, v[1], 0) *
			QuatNorm(0, 0, v[0]) *
			QuatNorm(0, 90 * DEG_TO_RAD, 0) *
			QuatNorm(2.1, 0, 0);
		q.normalize();
//		q = slerp(q, QuatNorm(-0.269887, -0.512545, 0.311801, -0.753153), 0.5);
		vars[ofs + 9*4 + 0] = q.x;
		vars[ofs + 9*4 + 1] = q.y;
		vars[ofs + 9*4 + 2] = q.z;
		vars[ofs + 9*4 + 3] = q.w;

		// lKneeA
		vars[ofs + rotOfs + 0] = -lud(dofNames, dofs, "l_knee_euler_z");
		// rKneeA
		vars[ofs + rotOfs + 1] = -lud(dofNames, dofs, "r_knee_euler_z");
		// lElbowA
		vars[ofs + rotOfs + 2] = lud(dofNames, dofs, "l_elbow_euler_z") * 0.8;
		// lForearmA
		vars[ofs + rotOfs + 3] = lud(dofNames, dofs, "l_elbow_euler_y");
		// lWristA
		vars[ofs + rotOfs + 4] = -lud(dofNames, dofs, "l_wrist_euler_x");
		// rElbowA
		vars[ofs + rotOfs + 5] = -lud(dofNames, dofs, "r_elbow_euler_z") * 0.8;
		// rForearmA
		vars[ofs + rotOfs + 6] = -lud(dofNames, dofs, "r_elbow_euler_y");
		// lWristA
		vars[ofs + rotOfs + 7] = -lud(dofNames, dofs, "r_wrist_euler_x");

		// lClavicle
		vars[ofs + rotOfs + 8] = lud(dofNames, dofs, "l_scapula_euler_z");
		vars[ofs + rotOfs + 9] = -(lud(dofNames, dofs, "l_scapula_euler_x") - 1.0);
		// rClavicle
		vars[ofs + rotOfs + 10] = lud(dofNames, dofs, "r_scapula_euler_z");
		vars[ofs + rotOfs + 11] = -(lud(dofNames, dofs, "r_scapula_euler_x") + 1.0);

		frame++;
	}

	cout << "loaded " << frame << " frames" << endl;
	in.close();

	init(frame, iNumVars);
	for (i=0; i < vars.size(); i++)
		data[i] = vars[i];

	delete []dofs;
	return true;
}



bool UMocapPoses::loadPoses(char *fname) {
	FILE *f;
	int frame, i, j;

	if (!openFile(&f, fname, "rb", "pose file"))
		return false;

	fread(&i, sizeof(int), 1, f);
	if (i != 0) {
		cout << "unknown version: " << i << endl;
		fclose(f);
		return false;
	}

	fread(&numFrames, sizeof(int), 1, f);
	fread(&i , sizeof(int), 1, f);
	if (numVars != i) {
		cout << "incompatible # of variables: expected " << numVars <<
			"; found " << i << endl;
		fclose(f);
		return false;
	}
	cout << "loading " << numFrames << " frames and " << numVars << " vars" << endl;
	init(numFrames, numVars);

	fread(data, sizeof(double), numFrames * numVars, f);
	fclose(f);

	for (i=0; i < NUM_QUATS; i++)
		fixQ[i] = QuatNorm();
	memset(fixR, 0, sizeof(double)*NUM_ROTS);

	return true;
}

bool UMocapPoses::savePoses(char *fname) {
	FILE *f;
	int frame, i, j;

	if (!openFile(&f, fname, "wb", "pose file"))
		return false;

	i = 0;
	fwrite(&i, sizeof(int), 1, f);
	fwrite(&numFrames, sizeof(int), 1, f);
	fwrite(&numVars , sizeof(int), 1, f);
	fwrite(data, sizeof(double), numFrames * numVars, f);
	fclose(f);

	return true;
}

void UMocapPoses::zeroFix() {
	int i;
	for (i=0; i < NUM_QUATS; i++)
		fixQ[i] = QuatNorm();
	memset(fixR, 0, sizeof(double)*NUM_ROTS);
	memset(fixT, 0, sizeof(Vec3d)*NUM_TRANS);
}

void UMocapPoses::calcFix(Skeleton *tSkel, Skeleton *mcSkel) {
	QuatNorm q;
	double d;
	Vec3d v;
	int i;

	for (i=0; i < NUM_QUATS; i++) {
		q = mcSkel->transforms.getT(QUAT_NAMES[i])->curCoord.q;
		q.w *= -1;
		fixQ[i] = tSkel->transforms.getT(QUAT_NAMES[i])->curCoord.q * q;
	}
	for (i=0; i < NUM_ROTS; i++) {
		d = ((SkelEulerRotation*)mcSkel->transforms.getT(ROT_NAMES[i]))->curAngle;
		fixR[i] = ((SkelEulerRotation*)tSkel->transforms.getT(ROT_NAMES[i]))->curAngle - d;
	}
	for (i=0; i < NUM_TRANS; i++) {
		v = mcSkel->transforms.getT(TRANS_NAMES[i])->curCoord.v;
		fixT[i] = tSkel->transforms.getT(TRANS_NAMES[i])->curCoord.v - v;
	}
}

void UMocapPoses::saveFix(const char *fname) {
	int i;
	FILE *f;
	if (!openFile(&f, fname, "wb", "mocap fix"))
		return;
	for (i=0; i < NUM_QUATS; i++) {
		fwrite(&fixQ[i].x, sizeof(double), 1, f);
		fwrite(&fixQ[i].y, sizeof(double), 1, f);
		fwrite(&fixQ[i].z, sizeof(double), 1, f);
		fwrite(&fixQ[i].w, sizeof(double), 1, f);
	}
	fwrite(fixR, sizeof(double), NUM_ROTS, f);
	fwrite(fixT, sizeof(Vec3d), NUM_TRANS, f);
	fclose(f);
}

void UMocapPoses::loadFix(const char *fname) {
	int i;
	FILE *f;
	if (!openFile(&f, fname, "rb", "mocap fix"))
		return;
	for (i=0; i < NUM_QUATS; i++) {
		fread(&fixQ[i].x, sizeof(double), 1, f);
		fread(&fixQ[i].y, sizeof(double), 1, f);
		fread(&fixQ[i].z, sizeof(double), 1, f);
		fread(&fixQ[i].w, sizeof(double), 1, f);
	}
	fread(fixR, sizeof(double), NUM_ROTS, f);
	fread(fixT, sizeof(Vec3d), NUM_TRANS, f);
	fclose(f);
}

QuatNorm eulerToQuat(double x, double y, double z) {
	QuatNorm q(x * DEG_TO_RAD, y * DEG_TO_RAD, -z * DEG_TO_RAD);
	return q;
}

void UMocapPoses::frameToSkel(int frame, Skeleton *skel) {
	QuatNorm q;
	double d;
	int i;
	int ind = frame * numVars;

	for (i=0; i < NUM_QUATS; i++) {
		if (strcmp(QUAT_NAMES[i]+1, "ShoulderCR") == 0) {
			Vec3d v = matToEulerXZY((fixQ[i] * 
				QuatNorm(data[ind+0], data[ind+1], 
				data[ind+2], data[ind+3])).toMatrixD());

/*			if (QUAT_NAMES[i][0] == 'r') {
				cout << QuatNorm(data[ind+0], data[ind+1], 
				data[ind+2], data[ind+3]) << endl;
				cout << QuatNorm(0, -v[2], 0) * QuatNorm(0, 0, -v[1]) * QuatNorm(-v[0], 0, 0) << endl << endl;
			}*/
			if (QUAT_NAMES[i][0] == 'l') {
				if (v[2] < -2.6) v[2] = -2.6;
				if (v[2] > 0.15) v[2] = 0.15;
			}
			else {
				if (v[2] < -2.2) v[2] = -2.2;
				if (v[2] > 0.3) v[2] = 0.3;
			}

			char tempStr[80];
			const char XZY[4] = "XZY";
			int j;
			strcpy(tempStr, QUAT_NAMES[i]);
			for (j=0; j < 3; j++) {
				tempStr[strlen(tempStr)-2] = XZY[j];
				SkelEulerRotation *trans = (SkelEulerRotation*)skel->transforms.getT(tempStr);
				if (!trans) {
					cout << "WARNING: unknown transform " << tempStr << endl;
					continue;
				}
				trans->curAngle = v[j];
			}
		}
		else {
			SkelQuatRotation *trans = (SkelQuatRotation*)skel->transforms.getT(QUAT_NAMES[i]);
			if (!trans) {
				cout << "WARNING: unknown transform " << QUAT_NAMES[i] << endl;
				continue;
			}
			trans->curQuat = fixQ[i] * QuatNorm(data[ind+0], data[ind+1], 
				data[ind+2], data[ind+3]);
		}
		ind += 4;
	}
	for (i=0; i < NUM_ROTS; i++) {
		SkelEulerRotation *trans = (SkelEulerRotation*)skel->transforms.getT(ROT_NAMES[i]);
		if (!trans) {
			cout << "WARNING: unknown transform " << ROT_NAMES[i] << endl;
			continue;
		}
		trans->curAngle = data[ind] + fixR[i];
		ind++;
	}
	for (i=0; i < NUM_TRANS; i++) {
		SkelTranslation *trans = (SkelTranslation*)skel->transforms.getT(TRANS_NAMES[i]);
		if (!trans) {
			cout << "WARNING: unknown transform " << TRANS_NAMES[i] << endl;
			continue;
		}
		trans->curVal = Vec3d(data[ind+0], data[ind+1], data[ind+2]) + fixT[i];
		if (strcmp(trans->name, "baseT") == 0) {
			trans->curVal[0] = 0;
			trans->curVal[1] = 0;		
		}
		ind += 3;
	}
}

void UMocapPoses::skelToFrame(int frame, Skeleton *skel) {
	QuatNorm q;
	double d;
	int i;
	int ind = frame * numVars;

	for (i=0; i < NUM_QUATS; i++) {
		SkelQuatRotation *trans = (SkelQuatRotation*)skel->transforms.getT(QUAT_NAMES[i]);
		if (!trans) {
			cout << "WARNING: unknown transform " << QUAT_NAMES[i] << endl;
			continue;
		}
		data[ind+0] = trans->curQuat.x;
		data[ind+1] = trans->curQuat.y;
		data[ind+2] = trans->curQuat.z;
		data[ind+3] = trans->curQuat.w;
		ind += 4;
	}
	for (i=0; i < NUM_ROTS; i++) {
		SkelEulerRotation *trans = (SkelEulerRotation*)skel->transforms.getT(ROT_NAMES[i]);
		if (!trans) {
			cout << "WARNING: unknown transform " << ROT_NAMES[i] << endl;
			continue;
		}
		data[ind] = trans->curAngle;
		ind++;
	}
	for (i=0; i < NUM_TRANS; i++) {
		SkelTranslation *trans = (SkelTranslation*)skel->transforms.getT(TRANS_NAMES[i]);
		if (!trans) {
			cout << "WARNING: unknown transform " << TRANS_NAMES[i] << endl;
			continue;
		}
		data[ind+0] = trans->curVal[0];
		data[ind+1] = trans->curVal[1];
		data[ind+2] = trans->curVal[2];
		ind += 3;
	}
}

/* old CSV version
void UMocapPoses::frameToSkel(int frame, Skeleton *skel, bool poseOnly) {
	SkelTranslation *transTrans;
	SkelQuatRotation *quatTrans;
	SkelEulerRotation *eulerTrans;
	int ofs = frame * numVars;

	quatTrans = (SkelQuatRotation*)skel->transforms.getT("baseQ");
	if (quatTrans) {
		quatTrans->curQuat = eulerToQuat(data[ofs+0], data[ofs+1], data[ofs+2]);
	}
	transTrans = (SkelTranslation*)skel->transforms.getT("baseT");
	if (transTrans) {
		Vec3d v = Vec3d(data[ofs+3], data[ofs+4], data[ofs+5]);
		transTrans->curVal = 0.001 * v;
	}
	ofs += 6;

	quatTrans = (SkelQuatRotation*)skel->transforms.getT("lHipQ");
	if (quatTrans) {
		quatTrans->curQuat = fixQ[0] * eulerToQuat(data[ofs+0], data[ofs+1], data[ofs+2]);
	}
	eulerTrans = (SkelEulerRotation*)skel->transforms.getT("lKneeA");
	if (eulerTrans) {
		eulerTrans->curAngle = -data[ofs+4] * DEG_TO_RAD + fixR[0];
	}
	quatTrans = (SkelQuatRotation*)skel->transforms.getT("lAnkleA");
	if (quatTrans) {
		quatTrans->curQuat = fixQ[1] * eulerToQuat(data[ofs+6], data[ofs+7], data[ofs+8]);
	}
	ofs += 15;

	quatTrans = (SkelQuatRotation*)skel->transforms.getT("rHipQ");
	if (quatTrans) {
		quatTrans->curQuat = fixQ[2] * eulerToQuat(data[ofs+0], data[ofs+1], data[ofs+2]);
	}
	eulerTrans = (SkelEulerRotation*)skel->transforms.getT("rKneeA");
	if (eulerTrans) {
		eulerTrans->curAngle = -data[ofs+4] * DEG_TO_RAD + fixR[1];
	}
	quatTrans = (SkelQuatRotation*)skel->transforms.getT("rAnkleA");
	if (quatTrans) {
		quatTrans->curQuat = fixQ[3] * eulerToQuat(data[ofs+6], data[ofs+7], data[ofs+8]);
	}
	ofs += 15;

	quatTrans = (SkelQuatRotation*)skel->transforms.getT("waistQ");
	if (quatTrans) {
		quatTrans->curQuat = fixQ[4] * eulerToQuat(data[ofs+0], data[ofs+1], data[ofs+2]);
	}
	quatTrans = (SkelQuatRotation*)skel->transforms.getT("abdomenQ");
	if (quatTrans) {
		quatTrans->curQuat = fixQ[5] * eulerToQuat(data[ofs+3], data[ofs+4], data[ofs+5]);
	}
	quatTrans = (SkelQuatRotation*)skel->transforms.getT("neckR");
	if (quatTrans) {
		quatTrans->curQuat = fixQ[6] * eulerToQuat(data[ofs+6], data[ofs+7], data[ofs+8]);
	}
	ofs += 12;

	quatTrans = (SkelQuatRotation*)skel->transforms.getT("lClavicleQ");
	if (quatTrans) {
		quatTrans->curQuat = fixQ[7] * eulerToQuat(data[ofs+0], data[ofs+1], data[ofs+2]);
	}
	quatTrans = (SkelQuatRotation*)skel->transforms.getT("lShoulderQ");
	if (quatTrans) {
		quatTrans->curQuat = fixQ[8] * eulerToQuat(data[ofs+3], data[ofs+4], data[ofs+5]);
	}
	eulerTrans = (SkelEulerRotation*)skel->transforms.getT("lElbowA");
	if (eulerTrans) {
		eulerTrans->curAngle = -data[ofs+6] * DEG_TO_RAD + fixR[2];
	}
	eulerTrans = (SkelEulerRotation*)skel->transforms.getT("lForearmA");
	if (eulerTrans) {
		eulerTrans->curAngle = data[ofs+10] * DEG_TO_RAD + fixR[3];
	}
	eulerTrans = (SkelEulerRotation*)skel->transforms.getT("lWristA");
	if (eulerTrans) {
		eulerTrans->curAngle = -data[ofs+14] * DEG_TO_RAD + fixR[4];
	}
	ofs += 18;
	
	quatTrans = (SkelQuatRotation*)skel->transforms.getT("rClavicleQ");
	if (quatTrans) {
		quatTrans->curQuat = fixQ[9] * eulerToQuat(data[ofs+0], data[ofs+1], data[ofs+2]);
	}
	quatTrans = (SkelQuatRotation*)skel->transforms.getT("rShoulderQ");
	if (quatTrans) {
		quatTrans->curQuat = fixQ[10] * eulerToQuat(data[ofs+3], data[ofs+4], data[ofs+5]);
	}
	eulerTrans = (SkelEulerRotation*)skel->transforms.getT("rElbowA");
	if (eulerTrans) {
		eulerTrans->curAngle = -data[ofs+6] * DEG_TO_RAD + fixR[5];
	}
	eulerTrans = (SkelEulerRotation*)skel->transforms.getT("rForearmA");
	if (eulerTrans) {
		eulerTrans->curAngle = data[ofs+10] * DEG_TO_RAD + fixR[6];
	}
	eulerTrans = (SkelEulerRotation*)skel->transforms.getT("rWristA");
	if (eulerTrans) {
		eulerTrans->curAngle = -data[ofs+14] * DEG_TO_RAD + fixR[7];
	}
	ofs += 18;
}*/