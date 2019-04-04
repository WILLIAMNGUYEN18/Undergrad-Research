#include "keeper.h"
#include "trimesh.h"
#include "cvgf.h"

extern CVGoalFunction *curCVGF;

Keeper::Keeper(int iNumVars) {
	numVars = iNumVars;
	curGradient = new double[numVars];

	int numObsPts = curCVGF->numObsPts();

	// calculate number of functions
	numFuncs = 3 * numObsPts * MAX_MM_PTS + 
		3 * 3 * curCVGF->cv->tm->numTris() * (curCVGF->numComponents+1);
	tempM = new double[numFuncs];
	if (!tempM) {
		cout << "WARNING: can't allocate " << (numFuncs*sizeof(double)) << " bytes for tempM!" << endl;
	}

	cout << "observed points: " << numObsPts << "; number of functions: " <<
		(3 * numObsPts * MAX_MM_PTS) << " + " << (3 * 3 * curCVGF->cv->tm->numTris() * (curCVGF->numComponents+1)) << 
		" = " << numFuncs << endl;
	zeroHessian = true;
}

Keeper::~Keeper() {
	delete []curGradient;
	delete []tempM;
}

double Keeper::calcHG(double vars[]) {
	zeroHessian = false;
	return curCVGF->updateHG(vars, this);
}

void Keeper::matVecMult(const double x[], double r[]) const {
	int i;

	memset(tempM, 0, sizeof(double) * numFuncs);

	if (!zeroHessian) {
		// tempM = J * x
		curCVGF->jacMult(x, tempM, this);
		// r = J' * tempM
		curCVGF->jacTMult(tempM, r, this);
	}
}
  
int Keeper::numVar() const {
	return numVars;
}

void Keeper::refresh() {
	// set Hessian, gradient to 0
	memset(curGradient, 0, sizeof(double) * numVar());
	zeroHessian = true;
}

const double* Keeper::g() const {
	return curGradient;
}

double Keeper::calculateRo(const double* p, const double newTheta, const double oldTheta) const {
  double num = oldTheta - newTheta;
  double denom = vecDot(numVar(), g(), p);
  double *Bp = new double[numVar()];
  matVecMult(p, Bp);
  denom += .5 * vecDot(numVar(), p, Bp);
  denom *= -1.;
  assert(denom > 0);
  
  double result = num/denom;
  return result; 
}


