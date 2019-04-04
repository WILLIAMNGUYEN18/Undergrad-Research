#include "ba.h"
#include "main_win.h"
#include "bodydemo.h"
#include "viewer.h"
#include <GL/glu.h>
#include <Fl/fl_file_chooser.h>
#include "trackball.h"
//#include <Fl/fl_draw.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Multiline_Output.H>

bool pickMode;

bool lightingEnabled;
//bool blueBack = true;
bool blueBack = false;

Vec3d bkgColor;

extern MainWin *mainWin;

void redrawV() {
	if (mainWin)
		mainWin->viewer->redraw();
}

void redrawVNow() {
	if (mainWin)
		mainWin->viewer->redrawNow();
}

void uiWait() {
	if (mainWin)
		Fl::wait(0);
}

void uiScreenshot(char *fname) {
	if (mainWin)
		mainWin->viewer->screenshot(fname);
}


Fl_Window *messageWin = NULL;
Fl_Multiline_Output *messageO = NULL;
Fl_Value_Slider* hvs, *wvs;

void initMessage() {
	if (!messageWin) {
		messageWin = new Fl_Window(290, 95, "");
		messageWin->begin();
			hvs = new Fl_Value_Slider(85, 5, 200, 35, "Height (cm)");
			hvs->type(5);
			hvs->minimum(120);
			hvs->maximum(220);
			hvs->step(1);
			hvs->value(180);
			hvs->align(FL_ALIGN_LEFT);

			wvs = new Fl_Value_Slider(85, 45, 200, 35, "Weight (kg)");
			wvs->type(5);
			wvs->minimum(50);
			wvs->maximum(150);
			wvs->step(1);
			wvs->value(100);
			wvs->align(FL_ALIGN_LEFT);
		messageWin->end();
//		messageO = new Fl_Multiline_Output(200, 40, 100, 40, "test");
	}
	messageWin->show();
}

void setMessage(double h, double w) {
	hvs->value(h);
	wvs->value(w);
}

Viewer::Viewer(int X, int Y, int W, int H, const char *L) : Fl_Gl_Window(X, Y, W, H, L) {
	controlMode = 0;
	pickMode = false;
	reInit = true;
	lightingEnabled = true;
	home();
	rotCenter = Vec3d(0, 0, 0);
	transScale = 0.01;
//	bkgColor = Vec3d(1, 1, 1);
	bkgColor = Vec3d(0, 0, 0);

	if (!can_do(mode() | FL_ALPHA))
		cout << "warning: alpha mode not supported" << endl;
	else
		mode(mode() | FL_ALPHA);

	relativeCam = false;
}

void Viewer::draw() {
	if (!valid() || reInit) {
		initGL();
		reInit = false;
	//... window size is in w() and h().
	}

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode( GL_MODELVIEW );

    glPushMatrix();

	camera.drawGL();

//	glbAxes();

	if (relativeCam) {
		double mat[16];
		relCamMat.getGLMatrix(mat);
		glMultMatrixd(mat);
	}

	glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
	glGetIntegerv(GL_VIEWPORT, viewPort);
	curModelViewMat = 
		Mat4d(modelMatrix[0],modelMatrix[1],modelMatrix[2],modelMatrix[3],
		modelMatrix[4],modelMatrix[5],modelMatrix[6],modelMatrix[7],
		modelMatrix[8],modelMatrix[9],modelMatrix[10],modelMatrix[11],
		modelMatrix[12],modelMatrix[13],modelMatrix[14],modelMatrix[15]).transpose();
	curProjectionMat = 
		Mat4d(projMatrix[0],projMatrix[1],projMatrix[2],projMatrix[3],
		projMatrix[4],projMatrix[5],projMatrix[6],projMatrix[7],
		projMatrix[8],projMatrix[9],projMatrix[10],projMatrix[11],
		projMatrix[12],projMatrix[13],projMatrix[14],projMatrix[15]).transpose();
	curDrawMat = curProjectionMat * curModelViewMat;

	bodydemoDraw();
    
    glPopMatrix();
}

void Viewer::redrawNow() {
//	Fl_Window *last = Fl_Window::current();
	make_current();
	draw();
	swap_buffers();
//	if (last)
//		last->make_current();
}

void Viewer::initGL() {
	glClearColor(bkgColor[0],bkgColor[1],bkgColor[2],0.0f);

	Vec4d lightPos(0, 10, 0, 0);
	lightPos = camera.lightRot.toMatrixD() * lightPos;
	GLfloat light_position[] = { lightPos[0], lightPos[1], lightPos[2], 0.0f };

	GLfloat mat_ambient[4];
	GLfloat mat_diffuse[4];
	GLfloat mat_specular[4];
	GLfloat mat_shininess;

    mat_ambient[0] = mat_ambient[1] = mat_ambient[2] = 0.2f;
//    mat_ambient[0] = 0.7*127.0/255.0;
//	mat_ambient[1] = 0.7*64.0/255.0;
//	mat_ambient[2] = 0.7*47.0/255.0;
    mat_ambient[3] = 1.0;
    mat_diffuse[0] = mat_diffuse[1] = 0.2f;
	mat_diffuse[2] = 0.5;
    mat_diffuse[3] = 1.0;
    mat_specular[0] = mat_specular[1] = mat_specular[2] = 0.1f;
	mat_specular[3] = 1.0;
//    mat_specular[0] = 0.5*255.0/255.0;
//	mat_specular[1] = 0.5*133.0/255.0;
//	mat_specular[2] = 0.5*217.0/255.0;
//	mat_specular[3] = 1.0;
    mat_shininess = 50.0;

	glViewport(0, 0, w(), h());
	double xy_aspect = (double)w() / (double)h();

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	if (pickMode) {
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);
		gluPickMatrix(pickX, pickY, 2, 2, viewport);
	}
	//double fov = 41.82016;
	//double fov = 25.525;
	double fov = 45;
	//double fov = 10;
	gluPerspective( fov, xy_aspect, 0.05, 1000.0 );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    
    glEnable( GL_DEPTH_TEST );
	
//	glEnable( GL_POLYGON_SMOOTH );
//    glEnable (GL_BLEND);
//    glBlendFunc (GL_SRC_ALPHA_SATURATE, GL_ONE);
//    glDepthFunc( GL_LEQUAL );
//    glClearAccum( 0.0, 0.0, 0.0, 0.0 );

	if (pickMode) {
		glDisable( GL_LIGHTING );
		glDisable( GL_LIGHT0 );

		glDisable( GL_COLOR_MATERIAL );
		
		glCullFace( GL_BACK );
		glEnable( GL_CULL_FACE );
	}
	else {
		glShadeModel( GL_SMOOTH );
		glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_ambient );
		if (blueBack)
			glMaterialfv( GL_BACK, GL_DIFFUSE, mat_diffuse);
		glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular );
		glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, &mat_shininess );
		glLightfv( GL_LIGHT0, GL_POSITION, light_position );
		GLfloat zero[4] = {0, 0, 0, 1};
//		glLightfv( GL_LIGHT0, GL_SPECULAR, zero);
		GLfloat diff[4] = {0.9f, 0.9f, 0.9f, 1};
		glLightfv( GL_LIGHT0, GL_DIFFUSE, diff);
		glEnable( GL_LIGHTING );
		glEnable( GL_LIGHT0 );
		glEnable( GL_NORMALIZE );

		if (blueBack)
			glColorMaterial( GL_FRONT, GL_AMBIENT_AND_DIFFUSE );
		else
			glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
		glEnable( GL_COLOR_MATERIAL );
		
		glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1.0 /*GL_FALSE*/ );
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	//	glCullFace( GL_BACK );
		glDisable( GL_CULL_FACE ); 

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}

/*

        glLoadIdentity();
        glViewport(0,0,w(),h());
        glOrtho(-10,10,-10,10,-20010,10000);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);*/
}

int Viewer::handle(int event) {
	QuatNorm q;
	double curX, curY;

	curX = Fl::event_x();
	curY = Fl::event_y();

	switch(event) {
	case FL_PUSH:
		if (!bodydemoClick(curX, curY, Fl::event_state(FL_ALT))) {
			ballW = w() / 2.0;
			ballH = h() / 2.0;
			ballX = ballW;
			ballY = ballH;

			dragX = Fl::event_x();
			dragY = Fl::event_y();
			if (dragX == ballX && dragY == ballY)
				dragA = 0;
			else
				dragA = atan2(dragX - ballX, dragY - ballY);

			if (controlMode == 0 )
				dragQ = camera.rot;
			else 
				dragQ = camera.lightRot;

			shiftDrag = (Fl::event_state() & FL_SHIFT) > 0;
			controlDrag = (Fl::event_state() & FL_CTRL) > 0;
			dragButton = Fl::event_button();
		}

		return 1;

	case FL_DRAG:
		if (!bodydemoDrag(curX, curY)) {
			if (shiftDrag) {
				if (dragButton == FL_LEFT_MOUSE)
					q = QuatNorm(0, PI * (curX - dragX) / ballW, 0);
				else if (dragButton == FL_MIDDLE_MOUSE)
					q = QuatNorm(PI * (curY - dragY) / ballH, 0, 0);
				else if (dragButton == FL_RIGHT_MOUSE) {
					double curA;
					if (curX == ballX && curY == ballY)
						curA = 0;
					else
						curA = atan2(curX - ballX, curY - ballY);
					q = QuatNorm(0, 0, curA - dragA);
				}

				if (controlMode == 0 )
					camera.rot = dragQ * q;  // backwards
				else {
					camera.lightRot = dragQ * q; // backwards
					invalidate();
				}
			}
			else if (controlDrag) {
			}
			else {
				if (dragButton == FL_LEFT_MOUSE) {
					trackball(dragQuat, 
						-(dragX - ballX) / ballW, 
						(dragY - ballY) / ballH, 
						-(curX - ballX) / ballW, 
						(curY - ballY) / ballH);

					q.x = dragQuat[0];
					q.y = dragQuat[1];
					q.z = -dragQuat[2];
					q.w = -dragQuat[3];

					if (controlMode == 0 )
						camera.rot = camera.rot * q; // backwards
					else {
						camera.lightRot = camera.lightRot * q; // backwards
						invalidate();
					}
				}
				else if (dragButton == FL_MIDDLE_MOUSE) {
					camera.trans[0] += (curX - dragX) * transScale;
					camera.trans[1] -= (curY - dragY) * transScale;
				}
				else if (dragButton == (FL_RIGHT_MOUSE)) {
					camera.trans[2] -= (curY - dragY) * transScale;
				}

				dragX = curX;
				dragY = curY;
			}
			redraw();
		}
		return 1;

	case FL_RELEASE:
		bodydemoRelease();
//		if (Fl::event_is_click() && !controlDrag)
//			checkHit(Fl::event_x(), Fl::event_y(), Fl::event_button());
		return 1;
/*
	case FL_FOCUS:
		return 1;

	case FL_KEYBOARD:
		if (Fl::event_key() == FL_Enter || Fl::event_key() == FL_KP_Enter) {
			processCommand(curCmd);
			curCmd[0] = 0;
			return 1;
		}
		text = Fl::event_text();
		if (strlen(text) > 0) {
			if (strlen(curCmd) + strlen(text) < 1023)
				strcat(curCmd, text);
			cout << text;
			return 1;
		}
		return Fl_Gl_Window::handle(event);*/

	default:
		// pass other events to the base class...
		return Fl_Gl_Window::handle(event);
	}
}

bool Viewer::screenshot(char *fname) {
	// make sure OpenGL context is current
	make_current();

	// read the pixels
	int numPixels = w()*h();
	unsigned char *pixels = new unsigned char[numPixels*3*sizeof(unsigned char)];
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadBuffer(GL_FRONT);
	glReadPixels(0, 0, w(), h(), GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// swap red and blue, because TGA format is stupid
	int i;
	for (i=0; i < numPixels; i++) {
		pixels[i * 3 + 0] ^= pixels[i * 3 + 2];
		pixels[i * 3 + 2] ^= pixels[i * 3 + 0];
		pixels[i * 3 + 0] ^= pixels[i * 3 + 2];
	}

	// get file name
	if (fname == NULL) {
		fname = fl_file_chooser("L",
		"*.tga",
		NULL);
	}
	if (fname == NULL)
		return false;

	// open the file
	FILE *fptr;
	fptr = fopen(fname, "wb");
	if (fptr == NULL) {
		cout << "can't open " << fname << endl;
		return false;
	}

	// create tga header
	putc(0,fptr);
	putc(0,fptr);
	putc(2,fptr);                         // uncompressed RGB
	putc(0,fptr); putc(0,fptr);
	putc(0,fptr); putc(0,fptr);
	putc(0,fptr);
	putc(0,fptr); putc(0,fptr);           // X origin
	putc(0,fptr); putc(0,fptr);           // y origin
	putc((w() & 0x00FF),fptr);
	putc((w() & 0xFF00) / 256,fptr);
	putc((h() & 0x00FF),fptr);
	putc((h() & 0xFF00) / 256,fptr);
	putc(24,fptr);                        // 24 bit bitmap
	putc(0,fptr);

	// write the data
	fwrite(pixels, w()*h()*3*sizeof(char), 1, fptr);
	fclose(fptr);

	delete []pixels;

	return true;
}

bool Viewer::saveRIB(char *fname) {
	// get file name
	if (fname == NULL) {
		fname = fl_file_chooser("L",
		"*.rib",
		NULL);
	}
	if (fname == NULL)
		return false;

	// open the file
	ofstream out(fname);
	if (!out.good()) {
		cout << "can't open " << fname << endl;
		return false;
	}

	out << "Display \"test.tif\" \"file\" \"rgb\"" << endl;
	out << "Format 320 240 1" << endl;

	// negate y axis
	out << "Transform [ 1 0 0 0  0 -1 0 0  0 0 1 0  0 0 0 1 ]" << endl;

//	gluPerspective( 45.0, xy_aspect, 0.5, 1000.0 );
	out << "Projection \"perspective\" \"fov\" [ 45 ]" << endl;
//	glTranslatef(0, 0, -10);
//	gluLookAt(0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	out << "Rotate 180 1 0 0" << endl;
	out << "Translate 0 0 -15" << endl;

	out << "WorldBegin" << endl;
    out << "LightSource \"ambientlight\" 1 \"intensity\" [0.2]" << endl;
	Vec3d lightPos(0, 10, 0);
	lightPos = camera.lightRot.toMatrixD() * lightPos;
	out << "LightSource \"distantlight\" 2 \"intensity\" [1.2] \"from\" [ " << lightPos << " ] \"to\" [ 0 0 0 ]" << endl;
//	out << "LightSource \"pointlight\" 2 \"intensity\" [1.2] \"from\" [ " << lightPos << " ]" << endl;

	saveTransformRIB(out);

//	out << "Surface \"plastic\"" << endl;
	out << "Surface \"shinymetal\"" << endl;

	int i;
//	for (i=mainWin->tools.size()-1; i >= 0; i--)
//		mainWin->tools[i]->drawRIB(out);

	out << "WorldEnd" << endl;
    out.close();
	return true;
}

void Viewer::saveTransformRIB(ostream &out) {
	camera.drawRIB(out);

	if (relativeCam) {
		Mat4d rm = relCamMat.inverse();
		out << "ConcatTransform [";
		int r, c;
		for (r=0; r < 4; r++)
			for (c = 0; c < 4; c++)
				out << rm[c][r] << " ";
		out << "]" << endl;
	}
}

void Viewer::home() {
	camera.trans = Vec3d(0, 0, -3);
	camera.rot = QuatNorm();
	camera.lightRot = QuatNorm();
}

void Viewer::loadView(char *fname) {
	ifstream in(fname);
	if (!in.good()) {
		cout << "can't open " << fname << endl;
		return;
	}
	in >> camera.rot;
	in >> camera.lightRot;
	in >> camera.trans;
	in.close();
	invalidate();
	redraw();
}