#ifndef CLOTH_MASTER
#define CLOTH_MASTER

#include "tool.h"

class ClothUI : public Tool {
  virtual void drawGL();
};

void initClothMaster();

#endif