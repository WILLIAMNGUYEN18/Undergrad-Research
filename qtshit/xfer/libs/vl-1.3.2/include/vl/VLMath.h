/*
	File:			VLMath.h

	Function:		Various math definitions for VL
					
	Author(s):		Andrew Willmott

	Copyright:		(c) 1995-2000, Andrew Willmott
 */

#ifndef __VL_MATH__
#define __VL_MATH__

#include <stdlib.h>

// --- Inlines ----------------------------------------------------------------

// additions to arithmetic functions

#ifdef VL_HAS_IEEEFP
#include <ieeefp.h>
#define vl_is_finite(X) finite(X)
#elif defined (__GNUC__) && defined(__USE_MISC)
#define vl_is_finite(X) finite(X)
#else
#define vl_is_finite(X) (1)
#endif

#ifdef VL_HAS_DRAND
inline Double vl_rand()
{ return(drand48()); }
#else
inline Double vl_rand()
{ return(rand() / (RAND_MAX + 1.0)); }
#endif

#ifndef __CMATH__
// GNU's complex.h defines its own abs(double)
// (as does Visual Studio .NET  -- BA)
/*#ifdef VL_HAS_ABSF
inline Float abs(Float x)
{ return (fabsf(x)); }
#endif
inline Double abs(Double x)
{ return (fabs(x)); }*/
#endif
#ifdef VL_HAS_ABSF
inline Float len(Float x)
{ return (fabsf(x)); }
#endif
inline Double len(Double x)
{ return (fabs(x)); }

inline Float sqrlen(Float r)
{ return(sqr(r)); }
inline Double sqrlen(Double r)
{ return(sqr(r)); }

inline Float mix(Float a, Float b, Float s)
{ return((1.0 - s) * a + s * b); }
inline Double mix(Double a, Double b, Double s)
{ return((1.0 - s) * a + s * b); }

inline Double sign(Double d)
{
	if (d < 0)
		return(-1.0);
	else
		return(1.0);
}

// useful routines

inline Void SetReal(Float &a, Double b)
{ a = b; }
inline Void SetReal(Double &a, Double b)
{ a = b; }

inline Bool IsPowerOfTwo(Int a)
{ return((a & -a) == a); };

template <class S, class T> inline Void ConvertVec(const S &u, T &v)
{
	for (Int i = 0; i < u.Elts(); i++)
		v[i] = u[i];
}

template <class T> inline Void ConvertVec(const T &u, T &v)
{ v = u; }

template <class S, class T> inline Void ConvertMat(const S &m, T &n)
{
	for (Int i = 0; i < m.Rows(); i++)
		ConvertVec(m[i], n[i]);
}

template <class T> inline Void ConvertMat(const T &m, T &n)
{ n = m; }


#endif
