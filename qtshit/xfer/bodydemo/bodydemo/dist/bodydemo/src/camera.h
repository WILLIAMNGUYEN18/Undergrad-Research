#ifndef CAMERA_H
#define CAMERA_H

#include "quatnorm.h"
#include "vec.h"

class Camera {
public:
	QuatNorm rot, lightRot;
	Vec3d trans;

	Camera();

	bool load(char *fname);
//	bool save(char *fname);

	void drawGL();
	void drawRIB(ostream &os);

	static Camera interp(Camera &cam0, Camera &cam1, double w);
};

#endif
