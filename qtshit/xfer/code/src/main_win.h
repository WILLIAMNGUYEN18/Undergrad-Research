// generated by Fast Light User Interface Designer (fluid) version 1.0103

#ifndef main_win_h
#define main_win_h
#include <FL/Fl.H>
#include "doppel2.h"
#include "viewer.h"
extern int dispMode;
#include <FL/Fl_Window.H>
#include "doppel2.h"
#include "Fl/fl_file_chooser.h"
#include <iostream>
#include <fstream>
#include <FL/Fl_Menu_Button.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Input.H>

class MainWin {
public:
  MainWin();
  Fl_Window *window;
  Fl_Menu_Button *resButton;
  static Fl_Menu_Item menu_resButton[];
  static Fl_Menu_Item *r320MI;
private:
  inline void cb_r320MI_i(Fl_Menu_*, void*);
  static void cb_r320MI(Fl_Menu_*, void*);
public:
  static Fl_Menu_Item *r640MI;
private:
  inline void cb_r640MI_i(Fl_Menu_*, void*);
  static void cb_r640MI(Fl_Menu_*, void*);
public:
  static Fl_Menu_Item *r720MI;
private:
  inline void cb_r720MI_i(Fl_Menu_*, void*);
  static void cb_r720MI(Fl_Menu_*, void*);
public:
  static Fl_Menu_Item *r290MI;
private:
  inline void cb_r290MI_i(Fl_Menu_*, void*);
  static void cb_r290MI(Fl_Menu_*, void*);
public:
  static Fl_Menu_Item *rOtherMI;
private:
  inline void cb_rOtherMI_i(Fl_Menu_*, void*);
  static void cb_rOtherMI(Fl_Menu_*, void*);
public:
  Viewer *viewer;
  Fl_Box *sepBox;
  Fl_Check_Button *whiteCB;
private:
  inline void cb_whiteCB_i(Fl_Check_Button*, void*);
  static void cb_whiteCB(Fl_Check_Button*, void*);
  inline void cb_Save_i(Fl_Button*, void*);
  static void cb_Save(Fl_Button*, void*);
  inline void cb_Load_i(Fl_Button*, void*);
  static void cb_Load(Fl_Button*, void*);
public:
  Fl_Button *smButton;
private:
  inline void cb_smButton_i(Fl_Button*, void*);
  static void cb_smButton(Fl_Button*, void*);
public:
  Fl_Button *flButton;
private:
  inline void cb_flButton_i(Fl_Button*, void*);
  static void cb_flButton(Fl_Button*, void*);
public:
  Fl_Button *hwButton;
private:
  inline void cb_hwButton_i(Fl_Button*, void*);
  static void cb_hwButton(Fl_Button*, void*);
public:
  Fl_Button *wfButton;
private:
  inline void cb_wfButton_i(Fl_Button*, void*);
  static void cb_wfButton(Fl_Button*, void*);
public:
  Fl_Button *rsButton;
private:
  inline void cb_rsButton_i(Fl_Button*, void*);
  static void cb_rsButton(Fl_Button*, void*);
  inline void cb__i(Fl_Button*, void*);
  static void cb_(Fl_Button*, void*);
public:
  Fl_Button *rlButton;
private:
  inline void cb_rlButton_i(Fl_Button*, void*);
  static void cb_rlButton(Fl_Button*, void*);
  inline void cb_1_i(Fl_Button*, void*);
  static void cb_1(Fl_Button*, void*);
public:
  Fl_Check_Button *videoCB;
private:
  inline void cb_videoCB_i(Fl_Check_Button*, void*);
  static void cb_videoCB(Fl_Check_Button*, void*);
  inline void cb_RIB_i(Fl_Button*, void*);
  static void cb_RIB(Fl_Button*, void*);
public:
  Fl_Button *bSkelButton;
private:
  inline void cb_bSkelButton_i(Fl_Button*, void*);
  static void cb_bSkelButton(Fl_Button*, void*);
public:
  Fl_Button *bTemplateButton;
private:
  inline void cb_bTemplateButton_i(Fl_Button*, void*);
  static void cb_bTemplateButton(Fl_Button*, void*);
public:
  Fl_Button *bDefButton;
private:
  inline void cb_bDefButton_i(Fl_Button*, void*);
  static void cb_bDefButton(Fl_Button*, void*);
public:
  Fl_Check_Button *lightingCB;
private:
  inline void cb_lightingCB_i(Fl_Check_Button*, void*);
  static void cb_lightingCB(Fl_Check_Button*, void*);
public:
  Fl_Input *cmdInput;
private:
  inline void cb_cmdInput_i(Fl_Input*, void*);
  static void cb_cmdInput(Fl_Input*, void*);
  inline void cb_Run_i(Fl_Button*, void*);
  static void cb_Run(Fl_Button*, void*);
public:
  Fl_Menu_Button *presetMB;
  void viewerResize(int x, int y, bool video);
  void updateResButton();
  void displayMode(int mode);
  void saveView();
  void loadView();
  vector<Tool*> tools;
  static void cb_button(Fl_Button* o, void* v);
  inline void cb_button_i(Fl_Button*o, void*);
  void addTool(Tool *tool);
  bool videoAspect;
  static void cb_procLabel(Fl_Menu_*m, void*s);
  char * fixLabel(const char *s);
};
#endif