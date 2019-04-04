/*
	propui.cpp

	Classes for user interfaces for property lists.

	Brett Allen
	May, 2002
*/

#include <fstream>
#include "propui.h"
#include <FL/Fl_Group.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_Multi_Browser.H>
#include <FL/Fl_Roller.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Value_Input.H>
#include <FL/fl_file_chooser.h>

PropUI::PropUI() {
	name = NULL;
	child = NULL;
	editCB = NULL;
}

PropUI::~PropUI() {
	if (name)
		free(name);
}

void PropUI::create(char *cName, int cX, int cY, Fl_Group *g) {
	if (name)
		free(name);
	name = strdup(cName);
	x = cX;
	y = cY;
	group = g;

//	if (g->w() < x + 150)
//		g->size(x + 150, g->h());
}

void PropUI::expand(char *name, void *data, int type, SLInterface *eParent) {
	closeChild();

	switch (type) {
	case SL_VOID:
		child = new PObject(*((SLInterface**)data));
		break;

	case SL_NAMETABLE_T:
		child = new PNameTable((NameTable*)data, -1);
		((PNameTable*)child)->expandable = true;
		break;

	case SL_DOUBLE:
		child = new PDouble((double*)data, eParent);
		break;

	case SL_QUATNORM:
		child = new PQuatNorm((QuatNorm*)data, eParent);
		break;

	case SL_VEC3D:
		child = new PVec3d((Vec3d*)data, eParent);
		break;
	}

	if (child) {
		child->create(name, x + 160, y, group);
		child->editCB = editCB;
	}

	group->redraw();
}

void PropUI::closeChild() {
	if (child) {
		child->closeChild();
		child->destroy();
		child = NULL;
	}
}

// PNameTable =============================================

PNameTable::PNameTable(NameTable *cNT, int cType, char *cLoadMask) {
	nt = cNT;
	ntType = cType;
	browser = NULL;

	loadMask[0] = 0;
	if (cLoadMask)
		strcpy(loadMask, cLoadMask);

	expandable = false;
	curExpanded = NULL;
}

void PNameTable::create(char *cName, int cX, int cY, Fl_Group *g) {
	g->begin();

	PropUI::create(cName, cX, cY, g);
	browser = new Fl_Multi_Browser(x, y + 10, 150, 215, name);
	browser->user_data(this);
	browser->callback((Fl_Callback*)PropUI::cbGlobal);
	browser->align(FL_ALIGN_TOP | FL_ALIGN_LEFT);
//	g->add(browser);

	allButton = new Fl_Button(x + 25, y + 230, 100, 20, "Select all");
	allButton->user_data(this);
	allButton->callback((Fl_Callback*)PropUI::cbGlobal);
//	g->add(allButton);

	noneButton = new Fl_Button(x + 25, y + 255, 100, 20, "Select none");
	noneButton->user_data(this);
	noneButton->callback((Fl_Callback*)PropUI::cbGlobal);
//	g->add(noneButton);

	if (loadMask[0]) {
		addButton = new Fl_Button(x + 25, y + 280, 100, 20, "Load...");
		addButton->user_data(this);
		addButton->callback((Fl_Callback*)PropUI::cbGlobal);
//		g->add(addButton);
		saveButton = new Fl_Button(x + 25, y + 280, 100, 20, "Save...");
		saveButton->user_data(this);
		saveButton->callback((Fl_Callback*)PropUI::cbGlobal);
	}
	else {
		addButton = NULL;
		saveButton = NULL;
	}

	g->end();

	copyFromData();
}

void PNameTable::destroy() {
	if (browser) {
		group->remove(browser);
		delete browser;
		browser = NULL;
	}
	if (allButton) {
		group->remove(allButton);
		delete allButton;
		allButton = NULL;
	}
	if (noneButton) {
		group->remove(noneButton);
		delete noneButton;
		noneButton = NULL;
	}
	if (addButton) {
		group->remove(addButton);
		delete addButton;
		addButton = NULL;
	}
}

void PNameTable::copyFromData() {
	int i;

	browser->clear();
	for (i=0; i < nt->size(); i++) {
		browser->add(nt->getName(i));
	}
}

bool PNameTable::isSelected(int i) {
	return (browser->selected(i+1) != 0);
}

void PNameTable::selectOne(int i) {
	browser->deselect();
	if (i > -1)
		browser->select(i+1);
}

int PNameTable::getSelected() {
	return browser->value() - 1;
}

void PNameTable::add(char *name, void *val) {
	browser->add(name);
	nt->add(name, val);
}

void PNameTable::callback(Fl_Widget* o, void* v) {
	if (o == addButton) {
		char *fname = fl_file_chooser("Load...", loadMask, NULL);
		if (fname) {
			add(fname, NULL);
		}
	}
	else if (o == saveButton) {
		if (browser->value() == 0)
			return;
		char *fname = fl_file_chooser("Save...", loadMask, NULL);
		if (fname) {
			SLInterface *curItem = (SLInterface*)nt->getItem(browser->value() - 1);
			ofstream os;
			if (!openOFStream(&os, fname, "saved object"))
				return;
			saveLoad.save(os, curItem);
			os.close();
		}
	}
	else if (o == allButton) {
		int i;
		for (i=0; i < browser->size(); i++)
			browser->select(i+1, 1);
	}
	else if (o == noneButton) {
		browser->deselect();
	}
	else {
		if (expandable && browser->value() > 0) {
			void *newExpand = nt->getItem(browser->value() - 1);
			if (newExpand != curExpanded) {
				curExpanded = newExpand;
				expand(nt->getName(browser->value()-1), &curExpanded, SL_VOID);
			}
		}
		redrawV();
	}

}

// PObject ================================================

PObject::PObject(SLInterface *obj) {
	object = obj;
	info = saveLoad.getClassInfo(obj);
	browser = NULL;
	curExpanded = NULL;
}

void PObject::create(char *cName, int cX, int cY, Fl_Group *g) {
	g->begin();

	PropUI::create(cName, cX, cY, g);
	browser = new Fl_Multi_Browser(x, y+10, 150, 290, name);
	browser->user_data(this);
	browser->callback((Fl_Callback*)PropUI::cbGlobal);
	browser->align(FL_ALIGN_TOP | FL_ALIGN_LEFT);
	g->end(); //add(browser);

	if (info) {
		int i;
		for (i=0; i < info->members.size(); i++) {
			browser->add(info->members.getName(i));
		}
	}
}

void PObject::destroy() {
	if (browser) {
		group->remove(browser);
		delete browser;
		browser = NULL;
	}
}

void PObject::callback(Fl_Widget* o, void* v) {
	SLProperty *newExpand = info->members.getT(browser->value() - 1);
	if (newExpand != curExpanded) {
		curExpanded = newExpand;
		expand(info->members.getName(browser->value() - 1), newExpand->getPtr(object), newExpand->type, object);
	}	
}

// PDouble ------------------------------------------------

PDouble::PDouble(double *dd, SLInterface *iParent) {
	d = dd;
	parent = iParent;
	input = NULL;
}

void PDouble::create(char *cName, int cX, int cY, Fl_Group *g) {
	PropUI::create(cName, cX, cY, g);

	g->begin();
	input = new Fl_Value_Input(x+20, y+10, 120, 25);
	input->user_data(this);
	input->callback((Fl_Callback*)PropUI::cbGlobal);
	input->align(FL_ALIGN_LEFT);
	input->minimum(-1000);
	input->maximum(1000);
	input->step(0.01);

	if (d)
		input->value(*d);

	g->end();
}

void PDouble::destroy() {
	if (input) {
		group->remove(input);
		delete input;
		input = NULL;
	}
}

void PDouble::callback(Fl_Widget* o, void* v) {
	if (!d)
		return;

	if (o == input) {
		(*d) = input->value();

		if (editCB)
			editCB(this);
//		redrawV();
	}
}

// PVec3d -------------------------------------------------

PVec3d::PVec3d(Vec3d *v, SLInterface *iParent) {
	vec = v;
	parent = iParent;

	int i;
	for (i=0; i < 3; i++)
		inputs[i] = NULL;
}

void PVec3d::create(char *cName, int cX, int cY, Fl_Group *g) {
	PropUI::create(cName, cX, cY, g);

	static const char LABELS[3][5] = {"X:", "Y:", "Z:"};
	int i;
	g->begin();
	for (i=0; i < 3; i++) {
		inputs[i] = new Fl_Value_Input(x+20, y+10+30*i, 120, 25, LABELS[i]);
		inputs[i]->user_data(this);
		inputs[i]->callback((Fl_Callback*)PropUI::cbGlobal);
		inputs[i]->align(FL_ALIGN_LEFT);
		inputs[i]->minimum(-10);
		inputs[i]->maximum(10);
		inputs[i]->step(0.01);

		if (vec)
			inputs[i]->value((*vec)[i]);
	}
	g->end();
}

void PVec3d::destroy() {
	int i;
	for (i=0; i < 3; i++) {
		if (inputs[i]) {
			group->remove(inputs[i]);
			delete inputs[i];
			inputs[i] = NULL;
		}
	}
}

void PVec3d::callback(Fl_Widget* o, void* v) {
	int i;

	if (!vec)
		return;

	for (i=0; i < 3; i++) {
		if (o == inputs[i]) {
			(*vec)[i] = ((Fl_Valuator*)o)->value();

			if (editCB)
				editCB(this);
//			redrawV();
		}
	}
}


PQuatNorm::PQuatNorm(QuatNorm *q, SLInterface *iParent) {
	quat = q;
	parent = iParent;

	int i;
	for (i=0; i < 3; i++)
		rollers[i] = NULL;
}

void PQuatNorm::create(char *cName, int cX, int cY, Fl_Group *g) {
	PropUI::create(cName, cX, cY, g);

	static const char LABELS[3][5] = {"X:", "Y:", "Z:"};
	int i;
	g->begin();
	for (i=0; i < 3; i++) {
		rollers[i] = new Fl_Roller(x+20, y+10+30*i, 120, 25, LABELS[i]);
		rollers[i]->user_data(this);
		rollers[i]->callback((Fl_Callback*)PropUI::cbGlobal);
		rollers[i]->type(FL_HORIZONTAL);
		rollers[i]->align(FL_ALIGN_LEFT);
		rollers[i]->minimum(-10);
		rollers[i]->maximum(10);
		rollers[i]->step(0.02);
//		g->add(rollers[i]);

		lastVals[i] = 0;
	}
	g->end();
}

void PQuatNorm::destroy() {
	int i;
	for (i=0; i < 3; i++) {
		if (rollers[i]) {
			group->remove(rollers[i]);
			delete rollers[i];
			rollers[i] = NULL;
		}
	}
}

void PQuatNorm::callback(Fl_Widget* o, void* v) {
	int i;

	if (!quat)
		return;

	for (i=0; i < 3; i++) {
		if (o == rollers[i]) {
			double angle = ((Fl_Valuator*)o)->value() - lastVals[i];
			lastVals[i] = ((Fl_Valuator*)o)->value();

			if (angle != 0) {
				QuatNorm deltaQ;

				switch (i) {
				case 0:
					deltaQ = QuatNorm(angle, 0, 0);
					break;
				case 1:
					deltaQ = QuatNorm(0, angle, 0);
					break;
				default:
					deltaQ = QuatNorm(0, 0, angle);
					break;
				}

				*quat = *quat * deltaQ;

				if (editCB)
					editCB(this);
//				redrawV();
			}
		}
	}
}


// PNormalizedNameTableD ----------------------------------

PNormalizedNameTableD::PNormalizedNameTableD(NameTableX<double> *n, SLInterface *iParent, double *r) {
	nt = n;
	parent = iParent;

	ranges = r;
	if (ranges) {
		normalized = false;
	}
	else
		normalized = true;
}

void PNormalizedNameTableD::create(char *cName, int cX, int cY, Fl_Group *g) {
	PropUI::create(cName, cX, cY, g);

	int i;
	inputs.resize(nt->size());
	g->begin();
	for (i=0; i < nt->size(); i++) {
		inputs[i] = new Fl_Value_Slider(g->x() + x+200, g->y() + y+18*i, 150, 18, nt->getName(i));
		inputs[i]->user_data(this);
		inputs[i]->callback((Fl_Callback*)PropUI::cbGlobal);
		inputs[i]->align(FL_ALIGN_LEFT);
		inputs[i]->type(FL_HORIZONTAL);
		if (ranges) {
			inputs[i]->minimum(-2.0*sqrt(ranges[i]));
			inputs[i]->maximum(2.0*sqrt(ranges[i]));
			inputs[i]->step(sqrt(ranges[i]) / 100);
		}
		else {
			inputs[i]->minimum(0);
			inputs[i]->maximum(1);
			inputs[i]->step(0.01);
		}
		inputs[i]->value(nt->getX(i));
	}
	g->end();
}

void PNormalizedNameTableD::destroy() {
	int i;
	for (i=0; i < inputs.size(); i++) {
		if (inputs[i]) {
			group->remove(inputs[i]);
			delete inputs[i];
		}
	}
	inputs.clear();
}

void PNormalizedNameTableD::loadVals() {
	int i;
	for (i=0; i < nt->size(); i++)
		inputs[i]->value(nt->getX(i));
}

void PNormalizedNameTableD::setNormalized(bool isNorm) {
	normalized = isNorm;
}

void PNormalizedNameTableD::callback(Fl_Widget* o, void* v) {
	int i;

	if (normalized) {
		for (i=0; i < inputs.size(); i++) {
			if (o == inputs[i])
				break;
		}

		if (i >= inputs.size())
			return;

		int primary = i;
		double remainder = 1.0 - inputs[i]->value();
		double sum = 0;

		for (i=0; i < inputs.size(); i++) {
			if (i != primary)
				sum += inputs[i]->value();
		}

		if (sum < 1e-6 || inputs.size() == 2) {
			for (i=0; i < inputs.size(); i++) {
				if (i != primary)
					inputs[i]->value(remainder / (inputs.size() - 1));
			}
		}
		else {
			double mult = remainder / sum;
			for (i=0; i < inputs.size(); i++) {
				if (i != primary)
					inputs[i]->value(inputs[i]->value() * mult);
			}
		}
	}

	for (i=0; i < nt->size(); i++)
		nt->getX(i) = inputs[i]->value();
	
	if (editCB)
		editCB(this);
//	redrawV();
}

void PNormalizedNameTableD::update() {
	int i;
	for (i=0; i < nt->size(); i++) {
		inputs[i]->value(nt->getX(i));
	}
}
