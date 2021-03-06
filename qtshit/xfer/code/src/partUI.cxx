// generated by Fast Light User Interface Designer (fluid) version 1.0103

#include "partUI.h"
#include "cli.h"
#include "skeleton.h"
#include "skinCalc.h"
#include "trimesh.h"
#include "uMaster.h"

inline void PartUI::cb_partCB_i(Fl_Check_Button*, void*) {
  showPart = partCB->value();
redrawV();
}
void PartUI::cb_partCB(Fl_Check_Button* o, void* v) {
  ((PartUI*)(o->parent()->user_data()))->cb_partCB_i(o,v);
}

inline void PartUI::cb_compVS0_i(Fl_Value_Slider*, void*) {
  updateFromSliders();
}
void PartUI::cb_compVS0(Fl_Value_Slider* o, void* v) {
  ((PartUI*)(o->parent()->user_data()))->cb_compVS0_i(o,v);
}

inline void PartUI::cb_compVS1_i(Fl_Value_Slider*, void*) {
  updateFromSliders();
}
void PartUI::cb_compVS1(Fl_Value_Slider* o, void* v) {
  ((PartUI*)(o->parent()->user_data()))->cb_compVS1_i(o,v);
}

inline void PartUI::cb_compVS2_i(Fl_Value_Slider*, void*) {
  updateFromSliders();
}
void PartUI::cb_compVS2(Fl_Value_Slider* o, void* v) {
  ((PartUI*)(o->parent()->user_data()))->cb_compVS2_i(o,v);
}

inline void PartUI::cb_compVS3_i(Fl_Value_Slider*, void*) {
  updateFromSliders();
}
void PartUI::cb_compVS3(Fl_Value_Slider* o, void* v) {
  ((PartUI*)(o->parent()->user_data()))->cb_compVS3_i(o,v);
}

inline void PartUI::cb_compVS4_i(Fl_Value_Slider*, void*) {
  updateFromSliders();
}
void PartUI::cb_compVS4(Fl_Value_Slider* o, void* v) {
  ((PartUI*)(o->parent()->user_data()))->cb_compVS4_i(o,v);
}

inline void PartUI::cb_Load_i(Fl_Button*, void*) {
  char *fname = fl_file_chooser("Load mesh",
 "*.{ply,obj}",
 NULL);
if (fname != NULL)
 loadTargetMesh(fname);
}
void PartUI::cb_Load(Fl_Button* o, void* v) {
  ((PartUI*)(o->parent()->user_data()))->cb_Load_i(o,v);
}

inline void PartUI::cb_targetCB_i(Fl_Check_Button*, void*) {
  showTarget = targetCB->value();
redrawV();
}
void PartUI::cb_targetCB(Fl_Check_Button* o, void* v) {
  ((PartUI*)(o->parent()->user_data()))->cb_targetCB_i(o,v);
}

PartUI::PartUI() {
  Fl_Window* w;
  { Fl_Window* o = partWindow = new Fl_Window(380, 377, "Parts");
    w = o;
    o->user_data((void*)(this));
    { Fl_Check_Button* o = partCB = new Fl_Check_Button(10, 5, 25, 25, "Show part");
      o->down_box(FL_DOWN_BOX);
      o->callback((Fl_Callback*)cb_partCB);
    }
    { Fl_Value_Slider* o = compVS0 = new Fl_Value_Slider(15, 65, 355, 25, "Component 0:");
      o->type(5);
      o->minimum(-2);
      o->maximum(2);
      o->callback((Fl_Callback*)cb_compVS0);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    { Fl_Value_Slider* o = compVS1 = new Fl_Value_Slider(15, 105, 355, 25, "Component 1:");
      o->type(5);
      o->minimum(-2);
      o->maximum(2);
      o->callback((Fl_Callback*)cb_compVS1);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    { Fl_Value_Slider* o = compVS2 = new Fl_Value_Slider(15, 145, 355, 25, "Component 2:");
      o->type(5);
      o->minimum(-2);
      o->maximum(2);
      o->callback((Fl_Callback*)cb_compVS2);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    { Fl_Value_Slider* o = compVS3 = new Fl_Value_Slider(15, 185, 355, 25, "Component 3:");
      o->type(5);
      o->minimum(-2);
      o->maximum(2);
      o->callback((Fl_Callback*)cb_compVS3);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    { Fl_Value_Slider* o = compVS4 = new Fl_Value_Slider(15, 225, 355, 25, "Component 4:");
      o->type(5);
      o->minimum(-2);
      o->maximum(2);
      o->callback((Fl_Callback*)cb_compVS4);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    { Fl_Button* o = new Fl_Button(120, 280, 135, 25, "Load target mesh...");
      o->callback((Fl_Callback*)cb_Load);
    }
    { Fl_Check_Button* o = targetCB = new Fl_Check_Button(15, 280, 25, 25, "Show target");
      o->down_box(FL_DOWN_BOX);
      o->callback((Fl_Callback*)cb_targetCB);
    }
    o->end();
  }
  name="Parts";
  { Fl_Window* o = skelWindow = new Fl_Window(1039, 373, "Skels");
    w = o;
    o->user_data((void*)(this));
    o->end();
  }
}

void PartUI::show() {
  partWindow->show();
}

void PartUI::updateFromSliders() {
  double vals[5];
vals[0] = compVS0->value();
vals[1] = compVS1->value();
vals[2] = compVS2->value();
vals[3] = compVS3->value();
vals[4] = compVS4->value();
updatePartWeights(vals, 5);
}
