#ifndef VIEWER_H
#define VIEWER_H

#include <Fl/Fl_Gl_Window.H>
#include "camera.h"
#include "doppel2.h"
#include "tool.h"

extern Vec3d bkgColor;

class MainWin;

class Viewer : public Fl_Gl_Window {
	int handle(int);

protected:
	double transScale;
	double dragX, dragY, dragA;
	int dragButton;
	double ballW, ballH, ballX, ballY;
	float dragQuat[4];
	bool shiftDrag;
	bool controlDrag;
	Tool *controlTool;

	double pickX, pickY;
	void initGL();

public:
	Camera camera;
	GLdouble modelMatrix[16], projMatrix[16];
	GLint viewPort[4];

	QuatNorm dragQ;
	Vec3d rotCenter;
	int controlMode;

	bool relativeCam;
	Mat4d relCamMat;

	Mat4d curDrawMat, curModelViewMat, curProjectionMat;
	bool pickMode, reInit;

	MainWin *mainWin;

	Viewer(int X, int Y, int W, int H, const char *L);

	bool screenshot(char *fname = NULL);
	bool saveRIB(char *fname = NULL);
	void saveTransformRIB(ostream &rib);

	void draw();
	void redrawNow();

	void home();

	Tool *checkHit(double x, double y, int button);
	void *uiNotify(const char *msg);

	void loadView(char *fname);
};

#endif
