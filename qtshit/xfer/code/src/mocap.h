#ifndef MOCAP_H
#define MOCAP_H

#define NUM_MAPPINGS 25

class Skeleton;

class MocapData {
public:
	int frames;
	int frameSize;
	double *data;
	int mappings[NUM_MAPPINGS];
	bool isText;

	MocapData();
	bool load(char *fname);
	bool loadTxt(char *fname);
	void toSkel(int frame, Skeleton *skel);
};

#endif
