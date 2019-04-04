#ifndef DOPPEL2_H
#define DOPPEL2_H

#ifdef WIN32
#include <windows.h>
#endif
#include "GL/glew.h"
#include "ba.h"
#include <math.h>
#include <iostream>
#include <vector>
using namespace std;
#include "vec.h"
#include "mat.h"
#include "quatnorm.h"
#include "trimesh_render.h"

class MainWin;
extern MainWin *mainWin;

void redrawV();
void redrawVNow();
void uiWait();
void uiScreenshot(char *fname = NULL);
void *uiNotify(const char *msg);

extern bool showSkel, showTemplate, showDef, lightingEnabled;
extern bool showMarkers, showFrames;
extern int dispMode;
extern Vec3d bkgColor;

#endif