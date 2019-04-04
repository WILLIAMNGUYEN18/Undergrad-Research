// generated by Fast Light User Interface Designer (fluid) version 1.0103

#ifndef cvUI_h
#define cvUI_h
#include <FL/Fl.H>
#include "doppel2.h"
#include "propui.h"
extern TriMesh *meshToRender;
#include <FL/Fl_Window.H>
#include "tool.h"
#include "Fl/fl_file_chooser.h"
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Scroll.H>

class CVUI : public Tool {
public:
  CVUI();
  Fl_Window *cvWindow;
  Fl_Round_Button *showCV;
private:
  inline void cb_showCV_i(Fl_Round_Button*, void*);
  static void cb_showCV(Fl_Round_Button*, void*);
public:
  Fl_Browser *exampleB;
private:
  inline void cb_exampleB_i(Fl_Browser*, void*);
  static void cb_exampleB(Fl_Browser*, void*);
public:
  Fl_Round_Button *showSelEx;
private:
  inline void cb_showSelEx_i(Fl_Round_Button*, void*);
  static void cb_showSelEx(Fl_Round_Button*, void*);
  inline void cb_Update_i(Fl_Button*, void*);
  static void cb_Update(Fl_Button*, void*);
public:
  Fl_Value_Slider *comp0VS;
private:
  inline void cb_comp0VS_i(Fl_Value_Slider*, void*);
  static void cb_comp0VS(Fl_Value_Slider*, void*);
public:
  Fl_Value_Slider *comp1VS;
private:
  inline void cb_comp1VS_i(Fl_Value_Slider*, void*);
  static void cb_comp1VS(Fl_Value_Slider*, void*);
  inline void cb_Edit_i(Fl_Button*, void*);
  static void cb_Edit(Fl_Button*, void*);
  inline void cb_Copy_i(Fl_Button*, void*);
  static void cb_Copy(Fl_Button*, void*);
public:
  Fl_Round_Button *showMarkers;
private:
  inline void cb_showMarkers_i(Fl_Round_Button*, void*);
  static void cb_showMarkers(Fl_Round_Button*, void*);
public:
  Fl_Value_Slider *comp2VS;
private:
  inline void cb_comp2VS_i(Fl_Value_Slider*, void*);
  static void cb_comp2VS(Fl_Value_Slider*, void*);
public:
  Fl_Value_Slider *comp3VS;
private:
  inline void cb_comp3VS_i(Fl_Value_Slider*, void*);
  static void cb_comp3VS(Fl_Value_Slider*, void*);
public:
  Fl_Value_Slider *comp4VS;
private:
  inline void cb_comp4VS_i(Fl_Value_Slider*, void*);
  static void cb_comp4VS(Fl_Value_Slider*, void*);
public:
  Fl_Value_Slider *comp5VS;
private:
  inline void cb_comp5VS_i(Fl_Value_Slider*, void*);
  static void cb_comp5VS(Fl_Value_Slider*, void*);
  inline void cb_Set_i(Fl_Button*, void*);
  static void cb_Set(Fl_Button*, void*);
  inline void cb_Set1_i(Fl_Button*, void*);
  static void cb_Set1(Fl_Button*, void*);
public:
  Fl_Window *skelWindow;
  Fl_Scroll *skelScroll;
private:
  inline void cb_Load_i(Fl_Button*, void*);
  static void cb_Load(Fl_Button*, void*);
  inline void cb_Save_i(Fl_Button*, void*);
  static void cb_Save(Fl_Button*, void*);
public:
  void show();
  virtual void drawGL();
  void * notify(const char *msg);
  void updateExampleB();
  void setComponents();
  void showSkel();
  static void editCB(PropUI*);
  virtual void clickGL(GLuint *nameBuf, double x, double y, int button);
  virtual void dragGL(int x, int y);
  virtual bool startDragGL(GLuint *nameBuf, int x, int y);
};
#endif
