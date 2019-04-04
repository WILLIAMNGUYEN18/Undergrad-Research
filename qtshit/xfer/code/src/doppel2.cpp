/*
	doppel_main.cpp

	Main function for doppel project.

	Brett Allen
	May 2002
*/

#include "doppel2.h"
#include "cli.h"
#include "main_win.h"
#include "saveload.h"
#include "viewer.h"
#include "doppel_config.h"

#include "skeleton.h"
#include "dof.h"

#include "uMaster.h"
#include "uniUI.h"
#include "partMaster.h"
#include "partUI.h"

MainWin *mainWin;
SaveLoad saveLoad;
bool activeUI;


/*
MeshUI *meshUI;
MarkerUI *markerUI;
SkelUI *skelUI;
SurfaceUI *surfaceUI;
SolverUI *solverUI;
DispUI *dispUI;
MatchUI *matchUI;
SMatchUI *smatchUI;
MovieUI *movieUI;
*/
//ClothUI *clothUI;
UniUI *uniUI;
PartUI *partUI;

void vLoadView(const char *params) {
	if (!activeUI)
		return;
	char fname[80];
	params = extractString(params, fname, 80);

	ifstream in(fname);
	if (!in.good()) {
	  cout << "can't open " << fname << endl;
	  return;
	}
	saveLoad.loadInstance(in, &mainWin->viewer->camera);
	in.close();
	mainWin->viewer->invalidate();
	mainWin->viewer->redraw();
}

void vRedrawNow(const char *params) {
	redrawVNow();
}

void vScreenshot(const char * params) {
	char fname[80];
	params = extractString(params, fname, 80);

	uiScreenshot(fname);
}

void initData() {
	registerFunction(vLoadView, "vLoadView");
	registerFunction(vRedrawNow, "vRedrawNow");
	registerFunction(vScreenshot, "vScreenshot");

/*	initTrimeshMaster();
	initSurfaceMaster();
	initMovieMaster();
*/	initUMaster();
	initPartMaster();

/*	Marker::registerClass(&saveLoad);
	MarkerSet::registerClass(&saveLoad);
	DoppelConfig::registerClass(&saveLoad);
*/	Camera::registerClass(&saveLoad);

	Skeleton::registerClass(&saveLoad);
	SkelTranslation::registerClass(&saveLoad);
	SkelQuatRotation::registerClass(&saveLoad);
	SkelEulerRotation::registerClass(&saveLoad);
	SkelPolarAxisRotation::registerClass(&saveLoad);
	SkelSymmetricTranslation::registerClass(&saveLoad);
	SkelPartialTransform::registerClass(&saveLoad);
	SkelCombinedTransform::registerClass(&saveLoad);

/*	KeyframeMovie::registerClass(&saveLoad);
*/	DofSet::registerClass(&saveLoad);

	dispMode = VM_SMOOTH;
//	pickMode = false;

	showSkel = false;
	showTemplate = false;
	showDef = false;
}

#include "vl/VLd.h"

int main(int argc, char **argv) {
	char scriptName[256];
	scriptName[0] = 0;
	activeUI = true;

	Fl::visual(FL_DOUBLE|FL_RGB);

	if (argc > 1) {
		int i;
		for (i=1; i < argc; i++) {
			if (argv[i][0] == '-') {
/*				if (argv[i][1] == 'r') {
					RES = atoi(argv[1]+2);
					cout << "resolution is " << RES << endl;
				}
*/			}
			else {
				strncpy(scriptName, argv[i], 255);
				cout << "script is '" << scriptName << "'" << endl;
				activeUI = false;
			}
		}
	}

	mainWin = new MainWin();
	mainWin->window->show();

	mainWin->updateResButton();

	initData();

//	cvUI = new CVUI();
//	mainWin->addTool(cvUI);

//	scUI = new SCUI();
//	mainWin->addTool(scUI);

//	clothUI = new ClothUI();
//	mainWin->addTool(clothUI);

	uniUI = new UniUI();
	mainWin->addTool(uniUI);

	partUI = new PartUI();
	mainWin->addTool(partUI);

	if (config.load("data/default.cf.txt"))
		config.copyToGlobal();

	int i;
	for (i=0; i < mainWin->tools.size(); i++)
		mainWin->tools[i]->init();

/*	initSolverMaster();
	initSMatchMaster();
	initCamMaster();
*/
//	data.pose[data.model]->addWidgets(data.poseUI->poseWindow);
//	data.bodyPose->addWidgets(bodyUI->bodyScroll);
//	data.bodyPose->updateWidgets();

	redrawVNow();

	processFile("scripts/startup.txt");

//	showSkel = true;
	showDef = true;
	mainWin->bSkelButton->value(1);
	mainWin->bDefButton->value(1);

/*	// a test
	SkelEulerRotation r;
	r.curAngle = 0;
	for (i=0; i < 3; i++) {
		cout << "axis " << i << endl;
		r.axis = i;
		r.updateCoord();
		r.updateDerivs();
		cout << r.curCoord.mat << endl;
		cout << *r.curCoord.deriv[0] << endl;
	}*/

	if (scriptName[0] != 0) {
		processFile(scriptName);
		return 0;
	}
	else
		return Fl::run();
}
