/*
	File:			VLfd.h

	Function:		Master header for a version of the VL library based on floats and doubles:
					vectors are made up of floats, and matrices of doubles. The various classes
					are named Vecf, Mat3d, SparseVecf, etc. To use this header you should link
					with -lvl.
					
	Author(s):		Andrew Willmott

	Copyright:		Copyright (c) 1995-1996, Andrew Willmott
 */

#ifndef __VLfd__
#define __VLfd__

#ifdef __VL__
#include "VLUndef.h"
#endif

#define VL_V_REAL Double
#define VL_V_SUFF(X) X ## d

#include "Mat2.h"
#include "Mat3.h"
#include "Mat4.h"
#include "Mat.h"
#include "SparseMat.h"
#include "Solve.h"

#include "Transform.h"

#include "VLUndef.h"
#define VL_V_REAL Float
#define VL_V_SUFF(X) X ## f
#define VL_M_REAL Double
#define VL_M_SUFF(X) X ## d

#include "Vec2.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Vec.h"
#include "SparseVec.h"
#include "Mixed.h"
#include "Solve.h"

#endif
