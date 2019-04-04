/*
	File:			VLf.h

	Function:		Master header for a version of the VL library based on
					floats. The various classes are named Vecf, Mat3f, 
					SparseVecf, etc.
					
	Author(s):		Andrew Willmott

	Copyright:		Copyright (c) 1995-1996, Andrew Willmott
 */


#ifndef __VLf__
#define __VLf__

#ifdef __VL__
#include "VLUndef.h"
#endif

#define VL_V_REAL Float
#define VL_V_SUFF(X) X ## f

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
