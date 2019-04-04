/*
	File:		VLgl.h
	
	Purpose:	Provides some handy wrappers for using vl with
				OpenGL.
 */

#ifndef __VLgl__
#define __VLgl__

inline Void glVertex(const Vec2f &a) 
{ glVertex2fv(a.Ref()); }

inline Void glVertex(const Vec3f &a) 
{ glVertex3fv(a.Ref()); }

inline Void glVertex(const Vec4f &a) 
{ glVertex4fv(a.Ref()); }

inline Void glColor(const Vec3f &a) 
{ glColor3fv(a.Ref()); }

inline Void glColor(const Vec4f &a) 
{ glColor4fv(a.Ref()); }

inline Void glNormal(const Vec3f &a) 
{ glNormal3fv(a.Ref()); }

inline Void glLoadMatrix(const Mat4f &m)
{ glLoadMatrixf(m.Ref()); }

inline Void glVertex(const Vec2d &a) 
{ glVertex2dv(a.Ref()); }

inline Void glVertex(const Vec3d &a) 
{ glVertex3dv(a.Ref()); }

inline Void glVertex(const Vec4d &a) 
{ glVertex4dv(a.Ref()); }

inline Void glColor(const Vec3d &a) 
{ glColor3dv(a.Ref()); }

inline Void glColor(const Vec4d &a) 
{ glColor4dv(a.Ref()); }

inline Void glNormal(const Vec3d &a) 
{ glNormal3dv(a.Ref()); }

inline Void glLoadMatrix(const Mat4d &m)
{ glLoadMatrixd(m.Ref()); }
// Note: glLoadMatrix[fd] expects matrices in column-major
// order, not row-order, so for things to work correctly, you
// should build with VL_ROW_ORIENT defined. Interestingly, 
// OpenGL internally operates with row vectors (just like
// the original GL), and transformation matrices stored in
// row-major order. However, externally they pretend that
// they use column vectors; everything still works, but
// because trans(Av) = trans(v) trans(A), it appears that
// matrices are stored in column major order.

#endif
