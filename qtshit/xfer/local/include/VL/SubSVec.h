/*
	File:			SubSVec.h

	Function:		Defines a scatter-gather sparse vector, i.e., a subvector
					of another vector or the row, column or diagonal of a
					sparse matrix.
					
	Author(s):		Andrew Willmott

	Copyright:		Copyright (c) 1995-1996, Andrew Willmott
 */

#ifndef __SubSVec__
#define __SubSVec__

#include "Vec.h"


class TSparseVec;

class TSubSVec
{
public:
					TSubSVec(Int start, Int length, TSparseVec *target,
							 Int span = 1);
					TSubSVec(const TSubSVec &v);
	
	inline Int		Elts() const { return(elts); };

	TSubSVec 		&operator = (const TSubSVec &v);
	TSubSVec		&operator = (const TSparseVec &v);
	Void			Store(TSparseVec &sv) const;

protected: 
	
	Int				elts;		// # of elements
	Int				start;		// start element in the target vector
	Int				colSpan;	// columns between elements
	Int				rowSpan;	// rows between elements
	TSparseVec		*target; 	// the sparse vector this is a subvector of
};

// --- Sub-vector operator ----------------------------------------------------
	
TSubSVec	sub(const TSparseVec &v, Int start, Int length); 

#endif
