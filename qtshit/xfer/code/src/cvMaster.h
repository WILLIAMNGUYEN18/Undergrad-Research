#ifndef CV_MASTER
#define CV_MASTER

void initCVMaster();
void cvClick(int x, int y, int button);
bool cvStartDrag(GLuint *nameBuf, int x, int y);
void cvDrag(int x, int y);

#endif