#ifndef U_MOCAP_H
#define U_MOCAP_H

#include "vec.h"

class Skeleton;

class UMocap {
public:
	int numFrames, numMarkers;
	Vec3d *data;

	UMocap();
	void init(int iNumFrames, int iNumMarkers);

	bool loadC3D(char *fname);
	bool loadCSV(char *fname);
};

class UMocapPoses {
public:
	int numFrames, numVars;
	double *data;
	QuatNorm *fixQ;
	double *fixR;
	Vec3d *fixT;

	UMocapPoses();
	void init(int iNumFrames, int iNumVars);

	bool load(const char *fname);
	bool loadAMC(const char *fname);
	bool loadCSV(const char *fname);
	bool loadDOF(const char *fname);
	bool loadPoses(char *fname);
	bool savePoses(char *fname);
	void zeroFix();
	void calcFix(Skeleton *tSkel, Skeleton *mcSkel);
	void saveFix(const char *fname);
	void loadFix(const char *fname);
	void frameToSkel(int frame, Skeleton *skel);
	void skelToFrame(int frame, Skeleton *skel);
};

#endif