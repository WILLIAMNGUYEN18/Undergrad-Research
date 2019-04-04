/*
	File:			VLd.h

	Function:		Master header for a version of the VL library based on doubles.
					The various classes are named Vecd, Mat3d, SparseVecd, etc.
					
	Author(s):		Andrew Willmott

	Copyright:		Copyright (c) 1995-1996, Andrew Willmott
 */


#ifndef __VLd__
#define __VLd__

#ifdef __VL__
#include "VLUndef.h"
#endif

#define VL_V_REAL Double
#define VL_V_SUFF(X) X ## d

#include "Vec2.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Vec.h"
#include "SparseVec.h"

#include "Mat2.h"
#include "Mat3.h"
#include "Mat4.h"
#include "Mat.h"
#include "SparseMat.h"
#include "Solve.h"
#include "Transform.h"

#endif
