#include <fstream>
using namespace std;

#include "doppel2.h"
#include "cli.h"
#include "liMaster.h"
#include "trimesh.h"
#include "main_win.h"
#include "viewer.h"

bool showColor = true;
bool specular = false;
vector<char*> meshSeq;

extern TriMesh *meshToRender;
extern MainWin *mainWin;

void applySpecular() {
	GLfloat mat_specular[4];
	GLfloat mat_shininess;
	if (specular) {
	    mat_specular[0] = mat_specular[1] = mat_specular[2] = 0.35f;
		mat_specular[3] = 1.0;
		mat_shininess = 2048.0;
		glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular );
		glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, &mat_shininess );
	} 
	else {
	    mat_specular[0] = mat_specular[1] = mat_specular[2] = 0.1f;
		mat_specular[3] = 1.0;
		mat_shininess = 50.0;
		glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular );
		glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, &mat_shininess );
	}
}

void setMesh(int frame) {
	if (meshToRender) {
		delete meshToRender;
		meshToRender = NULL;
	}

	if (frame < 0 || frame >= meshSeq.size())
		return;

	TriMesh *tm = new TriMesh();
	cout << "loading " << meshSeq[frame] << endl;
	tm->loadFile(meshSeq[frame]);
	int i;
	for (i=0; i < tm->numPts(); i++) {
		tm->getPt(i) -= Vec3d(25.124108, 69.291807, 1167.988814);
		tm->getPt(i) *= 1.0/1000.0;
	}
	tm->calcNormals();
	if (!showColor) {
		for (i=0; i < tm->numPts(); i++)
			tm->getPtColor(i) = Vec3d(0.6, 0.6, 0.6);
	}

	meshToRender = tm;
}

void liLoadPlySeq(const char *params) {
	char fname[256], curFN[256];
	int start = 0;
	int step = 1;
	int stop = 10000000;
	params = extractString(params, fname, 256);
	params = extractInt(params, &start);
	params = extractInt(params, &step);
	params = extractInt(params, &stop);

	int index;

	for (index = 0; index < meshSeq.size(); index++)
		delete meshSeq[index];
	meshSeq.clear();

	int count = 0;
	index = start;
	while (1) {
		sprintf(curFN, fname, index);

		FILE *f;
		if (!openFile(&f, curFN, "r")) {
			break;
		}
		fclose(f);
		meshSeq.push_back(strdup(curFN));

		count++;
		index += step;
		if (index >= stop)
			break;
	}
	cout << "loaded " << count << " meshes" << endl;
}

void liFrame(const char *params) {
	int frame = -1;
	params = extractInt(params, &frame);

	setMesh(frame);
	redrawV();
}

void liSaveAnim(const char *params) {
	char fname[256], curFN[256];
	vector<Camera> cameras;
	vector<int> frames;
	int i;

	fname[0] = 0;
	params = extractString(params, fname, 256);
	while (1) {
		i = -1;
		params = extractInt(params, &i);
		if (i == -1)
			break;

		curFN[0] = 0;
		params = extractString(params, curFN);
		if (curFN[0] == 0)
			break;

		Camera c = mainWin->viewer->camera;
		if (!c.load(curFN)) {
			cout << "WARNING: couldn't load " << curFN << endl;
		}
		cameras.push_back(c);
		frames.push_back(i);
	}

	if (fname[0] == 0)
		strcpy(fname, "frame%04d.tga");

	if (frames.size() == 0) {
		for (i=0; i < meshSeq.size(); i++) {
			setMesh(i);
			redrawVNow();

			sprintf(curFN, fname, i);
			uiScreenshot(curFN);
			uiWait();
		}
	}
	else {
		i = frames[0];
		int set;
		for (set=1; set < frames.size(); set++) {
			for (; i < frames[set]; i++) {
				double weight = 1.0 * (i - frames[set-1]) / (frames[set] - frames[set-1]);
				cout << i << " " <<  weight << endl;
				mainWin->viewer->camera = Camera::interp(cameras[set-1], cameras[set], weight);
//				mainWin->viewer->reInit = true;
//				applySpecular();

				setMesh(i);
				redrawVNow();

				sprintf(curFN, fname, i);
				uiScreenshot(curFN);
				uiWait();
			}
		}
	}
//	meshToRender = NULL;
}

void liColor(const char *params) {
	params = extractBool(params, &showColor);
	params = extractBool(params, &specular);
	applySpecular();
}

void initLiMaster() {
	registerFunction(liLoadPlySeq, "liLoadPlySeq");
	registerFunction(liFrame, "liFrame");
	registerFunction(liSaveAnim, "liSaveAnim");
	registerFunction(liColor, "liColor");
}