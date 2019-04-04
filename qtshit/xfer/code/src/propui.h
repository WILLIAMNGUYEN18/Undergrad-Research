/*
	propui.h

	Classes for user interfaces for property lists.

	Brett Allen
	May, 2002
*/

#ifndef PROPUI_H
#define PROPUI_H

#include <FL/Fl.H>
#include <FL/Fl_Widget.H>
#include <FL/Fl_Button.H>
#include "doppel2.h"
#include "saveload.h"

class Fl_Group;
class Fl_Multi_Browser;
class Fl_Slider;
class Fl_Value_Input;
class Fl_Roller;

class PropUI {
protected:
	Fl_Group *group;
	char *name;
	int x, y;
	PropUI *child;

public:
	void (*editCB)(PropUI*);

	PropUI();
	virtual ~PropUI();

	virtual void create(char *cName, int cX, int cY, Fl_Group *g);
	virtual void destroy() { }

	virtual void copyToData() { }
	virtual void copyFromData() { }

	static void cbGlobal(Fl_Widget* o, void* v) {
		if (o->user_data())
			((PropUI*)(o->user_data()))->callback(o,v);
	};
	virtual void callback(Fl_Widget* o, void* v) { }

	void expand(char *name, void *data, int type, SLInterface *parent = NULL);
	void closeChild();
};

class PNameTable : public PropUI {
private:
	NameTable *nt;
	int ntType;
	Fl_Button *allButton, *noneButton, *addButton, *saveButton;
	void *curExpanded;
	Fl_Multi_Browser *browser;

public:
	char loadMask[256];
	bool expandable;

	PNameTable(NameTable *cNT, int cType, char *cLoadMask = NULL);

	virtual void create(char *name, int x, int y, Fl_Group *g);
	virtual void destroy();

	virtual void copyFromData();

	bool isSelected(int i);
	void selectOne(int i);
	int getSelected();
	void add(char *name, void *val);

	virtual void callback(Fl_Widget* o, void* v);
};

class PObject : public PropUI {
protected:
	SLInterface *object;
	SLClassInfo *info;
	Fl_Multi_Browser *browser;
	void *curExpanded;

public:
	PObject(SLInterface *obj);

	virtual void create(char *name, int x, int y, Fl_Group *g);
	virtual void destroy();

	virtual void callback(Fl_Widget* o, void* v);
};

class PDouble : public PropUI {
protected:
	double *d;
	Fl_Value_Input *input;
	SLInterface *parent;

public:
	PDouble(double *dd, SLInterface *iParent = NULL);

	virtual void create(char *name, int x, int y, Fl_Group *g);
	virtual void destroy();

	virtual void callback(Fl_Widget* o, void* v);
};

class PVec3d : public PropUI {
protected:
	Vec3d *vec;
	Fl_Value_Input *inputs[3];
	SLInterface *parent;

public:
	PVec3d(Vec3d *v, SLInterface *iParent = NULL);

	virtual void create(char *name, int x, int y, Fl_Group *g);
	virtual void destroy();

	virtual void callback(Fl_Widget* o, void* v);
};

class PQuatNorm : public PropUI {
protected:
	QuatNorm *quat;
	Fl_Roller *rollers[3];
	double lastVals[3];
	SLInterface *parent;

public:
	PQuatNorm(QuatNorm *q, SLInterface *iParent = NULL);

	virtual void create(char *name, int x, int y, Fl_Group *g);
	virtual void destroy();

	virtual void callback(Fl_Widget* o, void* v);
};

class PNormalizedNameTableD : public PropUI {
protected:
	NameTableX<double> *nt;
	double *ranges;
	vector<Fl_Slider *> inputs;
	SLInterface *parent;
	bool normalized;

public:
	PNormalizedNameTableD(NameTableX<double> *n, SLInterface *iParent = NULL, double *r = NULL);

	virtual void create(char *name, int x, int y, Fl_Group *g);
	virtual void destroy();
	void loadVals();
	void setNormalized(bool isNorm);
	void update();

	virtual void callback(Fl_Widget* o, void* v);
};

#endif
