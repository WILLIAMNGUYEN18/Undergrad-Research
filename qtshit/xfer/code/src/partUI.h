// generated by Fast Light User Interface Designer (fluid) version 1.0103

#ifndef partUI_h
#define partUI_h
#include <FL/Fl.H>
#include "doppel2.h"
#include "partMaster.h"
#include <FL/Fl_Window.H>
#include "tool.h"
#include "Fl/fl_file_chooser.h"
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Button.H>

class PartUI : public Tool {
public:
  PartUI();
  Fl_Window *partWindow;
  Fl_Check_Button *partCB;
private:
  inline void cb_partCB_i(Fl_Check_Button*, void*);
  static void cb_partCB(Fl_Check_Button*, void*);
public:
  Fl_Value_Slider *compVS0;
private:
  inline void cb_compVS0_i(Fl_Value_Slider*, void*);
  static void cb_compVS0(Fl_Value_Slider*, void*);
public:
  Fl_Value_Slider *compVS1;
private:
  inline void cb_compVS1_i(Fl_Value_Slider*, void*);
  static void cb_compVS1(Fl_Value_Slider*, void*);
public:
  Fl_Value_Slider *compVS2;
private:
  inline void cb_compVS2_i(Fl_Value_Slider*, void*);
  static void cb_compVS2(Fl_Value_Slider*, void*);
public:
  Fl_Value_Slider *compVS3;
private:
  inline void cb_compVS3_i(Fl_Value_Slider*, void*);
  static void cb_compVS3(Fl_Value_Slider*, void*);
public:
  Fl_Value_Slider *compVS4;
private:
  inline void cb_compVS4_i(Fl_Value_Slider*, void*);
  static void cb_compVS4(Fl_Value_Slider*, void*);
  inline void cb_Load_i(Fl_Button*, void*);
  static void cb_Load(Fl_Button*, void*);
public:
  Fl_Check_Button *targetCB;
private:
  inline void cb_targetCB_i(Fl_Check_Button*, void*);
  static void cb_targetCB(Fl_Check_Button*, void*);
public:
  Fl_Window *skelWindow;
  void show();
  void updateFromSliders();
  virtual void drawGL();
};
#endif