#ifndef U_MASTER
#define U_MASTER

#include "tool.h"
#include "uSolver.h"
#include "uMocap.h"

const int NUM_CUR_COMPS = 20;
extern UDataSet dataSet;
extern double curComps[NUM_CUR_COMPS];
extern UMocapPoses mocapPoses;
extern USkin skin;
const int MAX_MAPPING_SIZE = 8;
extern int visSkel;
extern bool showTex;
extern bool softwareNM;

extern Vec3d *normalMap;

// skin painting stuff
extern int spTransform, spMode;
extern double spIntensity, spInnerRadius, spOuterRadius;
extern bool spAutoUpdate, spAutoExtend, spGeodesic;


void uShow(int toShow = -1000000);
void uColor(const char *params);
void uUpdateFromSkel();
void uShowComps(bool updateSkin = true);
void initUMaster();
bool loadMapping(char *fname, int &numSmallPts, int &numBigPts, int expectedPts, int *&mapPts, float *&mapWeights);
void updateMesh();

void updatePDDW();
void setPDDWJoint(const char *name);

#endif
