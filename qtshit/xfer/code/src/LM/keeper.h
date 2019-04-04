#ifndef KEEPER_H
#define KEEPER_H

#include <assert.h>
#include <stdio.h>
#include "linearSolver.h"




// This is a stub class that needs to be filled in by YOU.  The constructor takes the number of variables 
// (which we will refer to as n).  It hold the Hessian matrix (nxn), and the gradient (nx1).  The function
// refresh allows us to reuse this class without reallocating storage.
// Only the function calculateRo is already implemented.  
// Of course, the class is designed to allow you to store the Hessian sparsely!


class Keeper : public implicitMatrix {

 public:
	 int numVars, numFuncs;
	 double *curGradient, *tempM;
	 bool zeroHessian;

  Keeper(int iNumVars);
  ~Keeper();

  double calcHG(double vars[]);

  void matVecMult(const double x[], double r[]) const ;  // r = B*x
  
  int numVar() const; // return number of variables n

  void refresh(); // set Hessian, gradient to 0

  const double* g() const; // return gradient 

  double calculateRo(const double* p, const double newTheta, double const oldTheta) const; // already implemented

  
  
};

#endif
