//----------[ File Info ]----------------------------------
// File: solver.h
// Auth: Charles Gordon
// Date: 7/18/2000
// Copy: copyright (c) 2000, Charles Gordon
// Desc: 
//---------------------------------------------------------

#ifndef __SOLVER_HEADER__
#define __SOLVER_HEADER__

#include "vec.h"

//==========[ class IFunction ]================================================

class IFunction {
public:
		// returns the value of the function given the
		// current set of independent variables
	virtual double evaluateFunction( Vecd& variables ) = 0;

	virtual void startSolver() {}
	virtual void endSolver() {}
	virtual void solverStep() {}
};

//==========[ class IDifferentiableFunction ]==================================
// Desc: interface class used for defining functions that have first partial
//		 derivatives.  Used with NonlinearConjugateGradient and
//		 SteepestDescent.
//=============================================================================

class IDifferentiableFunction : public IFunction {
public:
		// returns the partial derivatives of the function with
		// respect to its independent variables
	virtual void evaluateGradient( Vecd& variables, Vecd& gradient ) = 0;
};

//==========[ class LBFGSSolver ]==============================================

class LBFGSSolver {
	IDifferentiableFunction*		mFunction;

	int								mNumVars;

	Vecd							mLowerBounds;
	Vecd							mUpperBounds;
	int*							mActiveBounds;
	Vecd							mGradient;

	double*							mWorkspace;
	int*							mIntWorkspace;
	char*							mTask;
	char*							mCharSave;
	int*							mBoolSave;
	int*							mIntSave;
	double*							mDoubleSave;


public:
	bool stopNow;

	LBFGSSolver( IDifferentiableFunction* function );
	~LBFGSSolver();

	void resize( int size );

	double solve( double factr, double pgtol, Vecd& x, int maxIterations = 1000000 );

	void setBound(int index, double lower, double upper);
};

#endif
