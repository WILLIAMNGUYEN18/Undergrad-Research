#include <vector>
#include <float.h>
#include "keeper.h"

using namespace std;

int numVec;
double *Z, *Z2;

inline double lmMin(double a, double b) {
	if (a < b)
		return a;
	return b;
}


// The following two functions project a solution and search direction to the boundary of the
// trust region.  USe the first for L2 norm, second for L-infinity norm.  Gets used by solve routine.

// This one is for L-2 norm
// d is search direction, p is current location
// KLT_TrackingContext::
void projectToTR(double *x, const double* d, const double* p, double trustRadius, int n) {

  double a,b,c, det, r1, tau;
  a = vecSqrLen(n,d);
  b = 2. * vecDot(n,d,p);
  c = vecSqrLen(n,p) - trustRadius*trustRadius;

  det = b*b-4.*a*c;
  assert(det >=0);
  det = sqrt(det);

  r1 = (-b + det) / (2.*a);
  if (r1>0) 
    tau = r1;
  else {
    tau = (-b - det) / (2.*a);
    assert(tau > 0);
  }

  assert(tau >0 && tau <1);
  vecAssign(n,x,d);
  vecTimesScalar(n,x,tau);
  vecAddEqual(n,x,p);
}


// This one is for L-infinity norm
// d is search direction, p is current location
void projectToTR2(double *x, const double* d, const double* p, double trustRadius, int n) {
  int i;
  double min = DBL_MAX, k1, k2;
  
  for (i=0; i<n; i++) {

    k1 = (trustRadius-p[i]) / d[i];
    k2 = (-trustRadius-p[i]) / d[i];
    if (k1>0 && k1<min)
	min = k1;
    if (k2>0 && k2<min)
      	min = k2;

  }

  
  vecAssign(n,x,d);
  vecTimesScalar(n,x,min);
  vecAddEqual(n,x,p);
}

// This routine is a modified conjugate-gradient to solve for one step of the optimization
// args:
// B: Holds Hessian, gradient
// x: we put the resultant calculated step here
// trustRadius: self-explanatory
// boundaryHit: We set this to true if the boundary of the trust region is hit, false otherwise
#define DELETE_SS delete[] p; delete[] r; delete[] d; delete[] Bd; delete[] pn;
void steihaugSolver(Keeper* B, double* x, const double trustRadius, bool* boundaryHit) {
	int n = B->numVar(), j=0;
	double *p = new double[n], *r = new double[n], *d = new double[n], *Bd = new double[n], *pn = new double[n];
	const double* g = B->g();
	double dBd, rr, alpha, beta, rrn, tmp;
	double error2; 
	double gradnorm2 = vecSqrLen(n,g);
	error2 = lmMin(.5, pow(gradnorm2,.25));
	error2 *= error2*gradnorm2;
	error2 = lmMin(1., error2);
	//  printf("error threshold = %.7f, grad norm : %.5f\n",error2, sqrt(gradnorm2));

	*boundaryHit = false;
	memset(p,0,n*sizeof(double));
	vecAssign(n, r, g);
	vecTimesScalar(n,r,-1.); 
	vecAssign(n, d, r);

	if ((tmp=vecSqrLen(n, r)) <= error2) {
		vecAssign(n,x,p);
		//    printf("error low %f, inner converge\n",tmp);
		DELETE_SS;
		return;
	}


	rr = vecDot(n,r,r);  
	while (1) {
		j++;
		B->matVecMult(d,Bd);
		dBd = vecDot(n,d,Bd);

		// case 1: negative curvature encountered, project to trust radius
		if (dBd < 0) {
			//      printf("negative curvature, projecting, inner converge\n");
			*boundaryHit = true;
			projectToTR(x,d,p,trustRadius, n); // use TR() for L-2 norm
			//      projectToTR2(x,d,p,trustRadius, n); // use TR() for L-2 norm
			DELETE_SS;
			//      printf("solved in %d  iterations \n",j);
			return;
		}

		alpha = rr / dBd;
		vecAssign(n,pn,d);
		vecTimesScalar(n,pn,alpha);
		vecAddEqual(n,pn,p);


		// case 2: trust region boundary hit
		if (vecSqrLen(n,pn) >= trustRadius*trustRadius) {  // USE THIS for L-2 norm
			//    if (fabs(vecMax(n,pn)) >= trustRadius) {
			//      printf("Hit trust radius boundary %f, inner converge\n", trustRadius);
			*boundaryHit = true;
			//      projectToTR2(x,d,p,trustRadius, n); // use TR() for L-2 norm
			projectToTR(x,d,p,trustRadius, n); // use TR() for L-2 norm
			DELETE_SS;
			//      printf("solved in %d  iterations \n",j);
			return;
		}


		vecTimesScalar(n, Bd, alpha); // Bd now alpha*Bd
		vecDiffEqual(n,r,Bd);

		// case 3: Converged to full quadratic step without hitting trust boundary
		if (vecSqrLen(n,r) < error2) {
			//      printf("gradient near nil, inner converge\n");
			vecAssign(n,x,pn);
			//      printf("final step, 2length %f\n",(tmp=vecSqrLen(n,pn)));
			DELETE_SS;
			//      printf("solved in %d  iterations \n",j);
			return;
		}

		rrn = vecDot(n,r,r);
		beta = rrn/rr;
		rr = rrn;
		vecTimesScalar(n,d,beta);
		vecAddEqual(n, d, r);
		vecAssign(n,p,pn);
	} 
}






// This is the actual solver that you call
// You will probably need other arguments to fill Hessian, gradient
// args:
// numVar: number of variables
// prec: Hoped for precision in solving for parameters (I use .01)
// double: initial AND maximum trust radius
// currSol: array of numVar doubles, initial solution for optimization (I use 10)
void solve(const int numVar, const double prec, double trustRadius, double* currSol, bool reAllocate, int maxIteration) {
	int iteration=0;
	bool toContinue = true;

	double *x = new double[numVar]; 
	double *newSol = new double[numVar];

	double* sol;
	static Keeper *keep1, *keep2;
	if (reAllocate || keep1 == NULL) {
		if (keep1)
			delete keep1;
		if (keep2)
			delete keep2;
		keep1 = new Keeper(numVar);
		keep2 = new Keeper(numVar); 
	}

	double currTheta, newTheta;
	currTheta = keep1->calcHG(currSol); // HERE, FILL HESSIAN, GRADIENT, RETURN FUNCTION VALUE, USING currSol

	do {
		//    fprintf(stdout,"starting iteration %d, %d active frames, %d variables\n",iteration, numVar);

		double ro;
		double maxStep = 1.; //  Initial value just forces into loop initially

		bool boundaryHit;
		while (maxStep > .05 && trustRadius >= prec) {
			steihaugSolver(keep1, x, trustRadius, &boundaryHit); 

			// eval solution
			maxStep = vecMax(numVar,x);
//			printf("max step %f\n",maxStep);
			vecAssign(numVar, newSol, currSol);
			vecAddEqual(numVar, newSol, x);

			newTheta = keep2->calcHG(newSol);  // HERE, FILL HESSIAN, GRADIENT, RETURN FUNCTION VALUE, USING newSol

			ro = keep1->calculateRo(x, newTheta, currTheta);
//			printf("ro: %f\n",ro);

			printf("old %f, new %f\n", currTheta, newTheta);
			if (ro < .25) {
				trustRadius *= .25;
				printf("reducing trustRadius to %f\n",trustRadius);
			}
			else if (ro > .75 &&  boundaryHit) 
				trustRadius = lmMin(2.*trustRadius, 1000.);

			if (ro <= 0) { // don't take step
				keep2->refresh();
			}
			else { // step is fine
				iteration++;
				printf("Taking step, iteration %d\n\n", iteration);
				Keeper* kswap = keep1;
				keep1 = keep2; 
				keep2 = kswap;
				keep2->refresh();
				currTheta = newTheta;
				vecAssign(numVar, currSol, newSol);
			}
		}

		if (trustRadius < prec || maxStep < prec || iteration > maxIteration)
			toContinue=false;
		else
			toContinue=true;
//		printf("\n\n");
	} while (toContinue);

	//  vecAssign(numVar, currSol, x);  
	delete[] x; delete[] newSol; 
}
