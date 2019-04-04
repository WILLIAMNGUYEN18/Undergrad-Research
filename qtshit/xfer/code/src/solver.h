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

//==========[ class ISolverUpdate ]============================================

class ISolverUpdate {
public:
	virtual bool solverUpdate() = 0;
};

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

//==========[ class ILinearSystem ]============================================
// Desc: interface class for defining linear systems of the form Ax = b, where
//		 A is a matrix, b is a vector of knowns and x is a vector of
//		 independent variables.  Used with LinearConjugateGradient solver.
//=============================================================================

class ILinearSystem {
public:

		// called before the solver begins to cache any necessary values
	virtual void startSolver() {}
		// called when the solver is done to remove any temporary memory
	virtual void endSolver() {}
		// called every time through the solver loop
	virtual void solverStep() {}

		// returns the b vector
	virtual void getB( Vecd& b ) = 0;
		// computes (A * x) and returns the result
	virtual void evaluateMatrixVectorProduct( Vecd& x, Vecd& result ) = 0;
};

//==========[ class LinearConjugateGradientSolver ]============================
// Desc: solves systems of the form Ax = b for values of x, the independent
//		 variables.
//=============================================================================

class LinearConjugateGradientSolver {
	
	ILinearSystem*	mLinearSystem;
	Vecd			mResidual;
	Vecd			mQ;
	Vecd			mDirection;
	Vecd			mTemp;
	int				mMaxVars;

	void resize( int size );

public:

	int			totalIterations;
	double		errorBound;

	LinearConjugateGradientSolver( ILinearSystem* ls );

	void solve( int maxIterations, double epsilon, Vecd& variables );
};

//==========[ class NonlinearConjugateGradientSolver ]=========================
// Desc: solves statements of the form f(x) = 0 for the optimal values of x,
//		 where f(x) is a nonlinear function.
//=============================================================================

class NonlinearConjugateGradientSolver {

	IDifferentiableFunction*	mFunction;

	Vecd		mResidual;
	Vecd		mOldResidual;
	Vecd		mDirection;
	Vecd		mTempVars;
	Vecd		mTempGradient;
	int			mMaxVars;

	void resize( int size );

public:

	NonlinearConjugateGradientSolver( IDifferentiableFunction* function );

	void solve( int maxIterations, double epsilon, Vecd& variables );
};

//==========[ class SteepestDescentSolver ]====================================
// Desc: solves statements of the form f(x) = 0 for the optimal values of x.
//		 Uses a simple gradient descent method.
//=============================================================================

class SteepestDescentSolver {

	IDifferentiableFunction*	mFunction;

	int			mMaxVars;
	Vecd		mResidual;

	void resize( int size );

public:

	int			totalIterations;
	double		errorBound;

	SteepestDescentSolver( IDifferentiableFunction* function );

	void solve( int maxIterations, double epsilon, Vecd& variables );
};

//==========[ class LBFGSSolver ]==============================================



class LBFGSSolver {
	IDifferentiableFunction*		mFunction;

	int								mNumVars;

	Vecd							mLowerBounds;
	Vecd							mUpperBounds;
	long*						mActiveBounds;
	Vecd							mGradient;

	double*							mWorkspace;
	long*							mIntWorkspace;
	char*							mTask;
	char*							mCharSave;
	long*							mBoolSave;
	long*							mIntSave;
	double*							mDoubleSave;


public:
	bool stopNow;

	LBFGSSolver( IDifferentiableFunction* function );
	~LBFGSSolver();

	void resize( int size );

	double solve( double factr, double pgtol, Vecd& x, int maxIterations = 1000000 );

	void setBound(int index, double lower, double upper);
};

//==========[ class SimulatedAnnealingSolver ]=================================

class SimulatedAnnealingSolver {
	IFunction*		mFunction;

	Vecd			mXOpt;
	Vecd			mXUpdate;
	Veci			mDimensionFreq;
	Vecd			mFTerm;
	Vecd			mStepVector;
	Vecd			mVaryingCriteria;

	double			mFOpt;
	double			mF;
	double			mNumCyclesBeforeTermTest;
	double			mNumCyclesBeforeTempUpdate;
	double			mNumDirectionCycles;
	double			mTemperature;
	double			mCoolingRate;

	void resize( int size );

public:

	SimulatedAnnealingSolver( IFunction* function );

	double solve( double epsilon, Vecd& variables );
};

#endif
