#include "solver.h"
#include "ba.h"

//==========[ Function Prototypes ]============================================

	// Fortran method for LBFGS solver
#include "f2c.h"
extern "C" { int setulb_(integer *n, integer *m, doublereal *x, 
	doublereal *l, doublereal *u, integer *nbd, doublereal *f, doublereal 
	*g, doublereal *factr, doublereal *pgtol, doublereal *wa, integer *
	iwa, char *task, integer *iprint, char *csave, logical *lsave, 
	integer *isave, doublereal *dsave, ftnlen task_len, ftnlen csave_len); }
/*#ifndef WIN32
#define SETULB setulb_
#endif
extern "C" { extern void SETULB(int*,int*,double*,double*,double*,int*,
								double*,double*,double*,double*,double*,
								int*,char*,int*,char*,int*,int*,double* ); }
*/
//==========[ LinearConjugateGradientSolver methods ]==========================

void LinearConjugateGradientSolver::resize( int size ) {
	if( mMaxVars >= size )
		return;

	mResidual.resize( size );
	mTemp.resize( size );
	mQ.resize( size );
	mDirection.resize( size );
	mMaxVars = size;
}

LinearConjugateGradientSolver::LinearConjugateGradientSolver( ILinearSystem* ls ) {
	mLinearSystem = ls;
	mMaxVars = 0;
}

void LinearConjugateGradientSolver::solve( int maxIterations, double epsilon, Vecd& variables ) {
	mLinearSystem->startSolver();

	int		numVars = variables.size();
	double	delta, deltaOld, deltaOrig, errorEpsilon;
	int		i = 0;
	double	alpha, beta;

	// make sure the cached vectors are the right length
	resize( numVars );

	// set the residual to b - Ax
	mLinearSystem->getB( mResidual );
	mLinearSystem->evaluateMatrixVectorProduct( variables, mTemp );
	mResidual -= mTemp;
	// the initial search direction is just the residual
	mDirection = mResidual;
	// calculate our error metric (the length of the residual squared)
	delta = mResidual * mResidual;
	deltaOrig = delta;
	errorEpsilon = epsilon * epsilon * deltaOrig;
	//errorEpsilon = epsilon;

	// repeat until we exceed maxIterations or meet our error bound
	while( i < maxIterations && delta > errorEpsilon ) {
		// set q = Ad
		mLinearSystem->evaluateMatrixVectorProduct( mDirection, mQ );
		alpha = delta / (mDirection * mQ);
		// update our estimate for x
		variables += alpha * mDirection;

		// restart if we have passed 50 iterations
		//if( i % 50 == 0 ) {
		//	mLinearSystem->getB( mResidual );
		//	mLinearSystem->evaluateMatrixVectorProduct( variables, mTemp );
		//	mResidual -= mTemp;
		//} else {
			mResidual -= alpha * mQ;
		//}

		// calculate the new search direction
		deltaOld = delta;
		delta = mResidual * mResidual;
		beta = delta / deltaOld;
		mDirection = mResidual + (beta * mDirection);

		i++;
	}

	totalIterations = i;
	errorBound = delta;

	mLinearSystem->endSolver();
}

//==========[ NonlinearConjugateGradientSolver methods ]=======================

void NonlinearConjugateGradientSolver::resize( int size ) {
	if( mMaxVars >= size )
		return;

	mMaxVars = size;
	mResidual.resize( size );
	mOldResidual.resize( size );
	mDirection.resize( size );
	mTempVars.resize( size );
	mTempGradient.resize( size );
}

NonlinearConjugateGradientSolver::NonlinearConjugateGradientSolver( IDifferentiableFunction* function ) {
	mFunction = function;
	mMaxVars = 0;
}

void NonlinearConjugateGradientSolver::solve( int maxIterations, double epsilon, Vecd& variables ) {
	static const int maxSecantIterations = 10;

	// parameters for the conjugate gradient method
	int		iterations = 0;
	int		secantIterations;
	double	deltaNew;
	double	deltaOld;
	double	deltaMid;
	double	deltaOriginal;
	double	beta;
	double	restartCondition;
	double	alpha;
	double	sigma;
	int		numVars = variables.size();

	resize( numVars );

	// various temporary variables;
	double	residuali;
	int		i;
	double	gradDotDir;
	double	secgradDotDir;

	mFunction->evaluateFunction( variables );
	mFunction->evaluateGradient( variables, mResidual );
	
	deltaOriginal = 0.0;
	for( i=0;i<numVars;i++ ) {
		residuali = -mResidual[i];

		mResidual[i] = residuali;
		mDirection[i] = residuali;
		deltaOriginal += residuali * residuali;
	}

	deltaNew = deltaOriginal;

	double errorEpsilon = epsilon * epsilon * deltaOriginal;
	while( iterations < maxIterations && deltaNew > errorEpsilon ) {
		for( i=0;i<numVars;i++ )
			mOldResidual[i] = mResidual[i];

		// do a line search using the secant method to find the value
		// that gets us closest to the minimum along this direction
		secantIterations = 0;
		sigma = 0.1;
		while( secantIterations < maxSecantIterations ) {
			for( i=0;i<numVars;i++ ) {
				mTempVars[i] = variables[i] + sigma * mDirection[i];
			}

			mFunction->evaluateFunction( variables );
			mFunction->evaluateGradient( mTempVars, mTempGradient );

			gradDotDir = 0.0;
			secgradDotDir = 0.0;
			for( i=0;i<numVars;i++ ) {
				if( secantIterations == 0 )
					gradDotDir += -mResidual[i] * mDirection[i];
				else
					gradDotDir += mResidual[i] * mDirection[i];
				secgradDotDir += mTempGradient[i] * mDirection[i];
			}

			if( secgradDotDir != gradDotDir )
				alpha = (-sigma * gradDotDir) / (secgradDotDir - gradDotDir);

			for( i=0;i<numVars;i++ ) {
				variables[i] += alpha * mDirection[i];
			}

			sigma = -alpha;

			mFunction->evaluateFunction( variables );
			mFunction->evaluateGradient( variables, mResidual );

			secantIterations++;
		}

		// calculate the new direction using Polak-Ribierre
		deltaOld = deltaNew;
		deltaNew = 0.0;
		deltaMid = 0.0;
		for( i=0;i<numVars;i++ ) {
			residuali = -mResidual[i];

			mResidual[i] = residuali;
			deltaNew += residuali * residuali;
			deltaMid += residuali * mOldResidual[i];
		}

		beta = (deltaNew-deltaMid) / deltaOld;
		restartCondition = 0.0;
		for( i=0;i<numVars;i++ ) {
			mDirection[i] = mResidual[i] + beta * mDirection[i];
			restartCondition += mDirection[i] * mResidual[i];
		}

		iterations++;

		// do a restart if we have tried more directions than we have
		// dimensions
		if( iterations % numVars == 0 || restartCondition <= 0.0 ) {
			mDirection = mResidual;
		}
	}
}

//==========[ SteepestDescentSolver methods ]==================================

void SteepestDescentSolver::resize( int size ) {
	if( mMaxVars >= size )
		return;

	mResidual.resize( size );
	mMaxVars = size;
}

SteepestDescentSolver::SteepestDescentSolver( IDifferentiableFunction* function ) {
	mFunction = function;
	mMaxVars = 0;
}

void SteepestDescentSolver::solve( int maxIterations, double epsilon, Vecd& variables ) {
	// parameters for calculation of steepest descent
	int					iterations = 0;
	double				deltaOrig = 0.0;
	double				delta = 0.0;
	double				alpha;
	int					numVars = variables.size();

	// temporary variables
	int					i;
	double				residuali;

	// make sure the cached vectors are correctly sized
	resize( numVars );

	mFunction->evaluateFunction( variables );
	mFunction->evaluateGradient( variables, mResidual );

	for( i=0;i<numVars;i++ ) {
		residuali = -mResidual[i];

		mResidual[i] = residuali;
		deltaOrig += residuali * residuali;
	}

	delta = deltaOrig;

	double errorEpsilon = epsilon * epsilon * deltaOrig;
	while( iterations < maxIterations && delta > errorEpsilon ) {
		// take a step in the opposite direction of the gradient
		alpha = 0.004; //0.016;
		variables += alpha * mResidual;

		// get the new residual
		mFunction->evaluateFunction( variables );
		mFunction->evaluateGradient( variables, mResidual );

		delta = 0.0;
		for( i=0;i<numVars;i++ ) {
			residuali = -mResidual[i];

			mResidual[i] = residuali;
			delta += residuali * residuali;
		}

		iterations++;
		mFunction->solverStep();
	}

	errorBound = delta;
	totalIterations = iterations;
}

//==========[ LBFGSSolver Methods ]============================================

LBFGSSolver::LBFGSSolver( IDifferentiableFunction* function ) {
	mFunction = function;

	mWorkspace = NULL;
	mIntWorkspace = NULL;
	mActiveBounds = NULL;

	mTask = new char[60];
	mCharSave = new char[60];
	mBoolSave = new long[4];
	mIntSave = new long[44];
	mDoubleSave = new double[29];

	mNumVars = 0;

	stopNow = false;
}

LBFGSSolver::~LBFGSSolver() {
	delete [] mTask;
	delete [] mCharSave;
	delete [] mBoolSave;
	delete [] mIntSave;
	delete [] mDoubleSave;

	if( mWorkspace != NULL )
		delete [] mWorkspace;

	if( mIntWorkspace != NULL )
		delete [] mIntWorkspace;

	if( mActiveBounds != NULL )
		delete [] mActiveBounds;
}

void LBFGSSolver::resize( int size ) {
	if( mNumVars < size ) {
		mGradient.resize( size );
		mLowerBounds.resize( size );
		mUpperBounds.resize( size );

		if( mWorkspace != NULL ) delete [] mWorkspace;
		if( mIntWorkspace != NULL ) delete [] mIntWorkspace;
		if( mActiveBounds != NULL ) delete [] mActiveBounds;

		mWorkspace = new double[14*size+315];
		baAssert(mWorkspace != NULL, "can't allocate LBFGS workspace", true);
		mIntWorkspace = new long[3*size];
		mActiveBounds = new long[size];

		int i;
		for( i=0;i<size;i++ )
			mActiveBounds[i] = 0;

		mNumVars = size;
	}
}

double LBFGSSolver::solve( double factr, double pgtol, Vecd& x, int maxIterations ) {
	int			i;
	integer			n = x.size();
	integer			m = 5;
	double		f = 0;
	integer			iprint = -1;
	bool		converged = false;
	int         numIterations = 0;

	resize( n );

	strcpy( mTask, "START" );
	for( i=5;i<60;i++ )
		mTask[i] = ' ';

//	for( i=0;i<n;i++ )
//		mActiveBounds[i] = 0;

	while( !converged ) {
//		cout << ".";
		// note: I have no idea what the last two parameters are for
		setulb_( &n, &m, x.getPointer(), mLowerBounds.getPointer(), mUpperBounds.getPointer(),
				mActiveBounds, &f, mGradient.getPointer(), &factr, &pgtol,
				mWorkspace, mIntWorkspace, mTask, &iprint, mCharSave,
				mBoolSave, mIntSave, mDoubleSave, 255, 255 );

		if( mTask[0] == 'F' && mTask[1] == 'G' ) {
			f = mFunction->evaluateFunction( x );
			mFunction->evaluateGradient( x, mGradient );
		} else if( mTask[0] == 'N' && mTask[1] == 'E' &&
					mTask[2] == 'W' && mTask[3] == '_' &&
					mTask[4] == 'X' ) {
			converged = false;
			mFunction->solverStep();
			if (stopNow)
				break;
			numIterations++;
			if (numIterations >= maxIterations) {
				break;
			}
		} else {
			converged = true;
		}
	}
/*
	if (converged)
		cout << "solver converged." << endl;
	else
		cout << "solve stopped." << endl;
*/
	return f;
}

//==========[ SimulatedAnnealingSolver Methods ]===============================

SimulatedAnnealingSolver::SimulatedAnnealingSolver( IFunction* function ) {
	mFunction = function;

	mFOpt = 0.0;
	mF = 0.0;
	mNumCyclesBeforeTermTest = 4;
	mNumCyclesBeforeTempUpdate = 100;
	mNumDirectionCycles = 20;
	mCoolingRate = 0.85;
	mTemperature = 10e13;
}

void SimulatedAnnealingSolver::resize( int size ) {
	mXOpt.resize( size );
	mXUpdate.resize( size );
	mDimensionFreq.resize( size );
	mStepVector.resize( size );
	mVaryingCriteria.resize( size );

	mFTerm.resize( mNumCyclesBeforeTermTest );
}

double SimulatedAnnealingSolver::solve( double epsilon, Vecd& x ) {
	// initialization
	int		j, k, m, u;
	bool	terminate = false;
	int		nDimensions = x.size();
	resize( nDimensions );

	for( u=0;u<nDimensions;u++ )
		mStepVector[u] = 10.0;

	mXOpt = x;
	mFOpt = mF = mFunction->evaluateFunction( x );
	
	for( u=0;u<nDimensions;u++ )
		mDimensionFreq[u] = 0;

	for( u=0;u<mNumCyclesBeforeTermTest;u++ )
		mFTerm[u] = mFOpt;

	j = k = m = 0;

	do {
		m = 0;
		while( m < mNumCyclesBeforeTempUpdate ) {

			j = 0;
			while( j < mNumDirectionCycles ) {

				mXUpdate = x;

				// choose a new direction
				for( int h=0;h<nDimensions;h++ ) {
					double r;
					
					double oldXVal = mXUpdate[h];
					do {
						r = boundedRand( -1.0f, 1.0f );
						mXUpdate[h] = oldXVal + r * mStepVector[h];
					} while( mXUpdate[h] > 2000 );

					double fUpdate = mFunction->evaluateFunction( mXUpdate );

					if( fUpdate <= mF ) { // if this direction is an improvement, accept it
						x = mXUpdate;
						mF = fUpdate;
						mDimensionFreq[h]++;
						if( mF < mFOpt ) {
							mFOpt = mF;
							mXOpt = x;
						}
					} else { // otherwise accept it with some probability
						double probOfAcceptance = exp( (mF - fUpdate) / mTemperature );
						r = boundedRand( 0.0f, 1.0f );
						if( r < probOfAcceptance ) {
							x = mXUpdate;
							mF = fUpdate;
							mDimensionFreq[h]++;
						}
					}
				}

				j++;
			}

			// calculate the new step vector
			for( u=0;u<nDimensions;u++ ) {
				if( mDimensionFreq[u] > 0.6 * mNumDirectionCycles ) {
					mStepVector[u] = mStepVector[u] * (1 + (/*mVaryingCriteria[u]*/2.0 * 
									 ((mDimensionFreq[u]/mNumDirectionCycles-0.6)/0.4)));
				} else if( mDimensionFreq[u] < 0.4 * mNumDirectionCycles ) {
					mStepVector[u] = mStepVector[u] / (1 + (/*mVaryingCriteria[u]*/2.0 *
									 ((0.4-(mDimensionFreq[u]/mNumDirectionCycles))/0.4)));
				}
			}

			for( u=0;u<nDimensions;u++ )
				mDimensionFreq[u] = 0;

			m++;
		}

		// update the temperature
		mTemperature *= mCoolingRate;
		mFTerm[k] = mF;

		x = mXOpt;
		mF = mFOpt;

		// check termination criterion
		if( mFTerm[k] - mFOpt <= epsilon ) {
			terminate = true;
			for( u=0;u<mNumCyclesBeforeTermTest;u++ ) {
				if( abs( mFTerm[k] - mFTerm[u] ) > epsilon ) {
					terminate = false;
					break;
				}
			}
		}

		k++;

		if( k == mNumCyclesBeforeTermTest )
			k = 0;

	} while( !terminate );

	x = mXOpt;
	return mFOpt;
}

void LBFGSSolver::setBound(int index, double lower, double upper) {
	mActiveBounds[index] = 2;
	mLowerBounds[index] = lower;
	mUpperBounds[index] = upper;
}