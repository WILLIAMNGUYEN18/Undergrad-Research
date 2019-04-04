#ifndef _TRIMESH_RENDER_H_
#define _TRIMESH_RENDER_H_

const int VM_SURF_ON       = 1;
const int VM_SURF_FLAT     = 2;
const int VM_SURF_SMOOTH   = 4;
const int VM_SURF_COLOR_V  = 8;
const int VM_SURF_COLOR_F  = 16;
const int VM_SURF_COLOR_S  = 32;
const int VM_WF_ON         = 64;
const int VM_WF_FLAT       = 128;
const int VM_WF_SMOOTH     = 256;
const int VM_WF_COLOR_V    = 512;
const int VM_WF_COLOR_F    = 1024;
const int VM_WF_COLOR_S    = 2048;
const int VM_TEX           = 4096;

const int VM_WIREFRAME     = VM_WF_ON;
const int VM_HIDDEN_LINE   = VM_SURF_ON | VM_WF_ON;
const int VM_FLAT          = VM_SURF_ON | VM_SURF_FLAT | VM_SURF_COLOR_V;
const int VM_SMOOTH        = VM_SURF_ON | VM_SURF_SMOOTH | VM_SURF_COLOR_V;

class TriMesh;

void renderTriMesh(TriMesh *tm, int viewMode, Vec3d bkg = Vec3d(0, 0, 0));

#endif