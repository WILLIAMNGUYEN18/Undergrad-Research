#ifndef TOOL_H
#define TOOL_H

#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#include <FL/Fl_Group.H>
#include <FL/Fl_Widget.H>
#include <iostream>
using namespace std;

class Tool {
public:
	char *name;

	Tool() { name = "tool"; }
	virtual ~Tool() { flush(); }

	virtual void flush() { }

	virtual void init() {}
	virtual void drawGL() {}
	virtual void drawRIB(ostream &rib) {}
	virtual void clickGL(GLuint *nameBuf, double x, double y, int button) {}
	virtual bool startDragGL(GLuint *nameBuf, int x, int y) { return false; }
	virtual void dragGL(int x, int y) {}
	virtual void show() {};
	virtual void *notify(const char *msg) { return NULL; }

	virtual void updateData() {};

	static void cb_controls(Fl_Widget* o, void* v) {
		((Tool*)(o->parent()->user_data()))->cb_controls_i(o,v);
	}
	virtual void cb_controls_i(Fl_Widget*, void*) { }
};

#endif
