#include "linearSolver.h"
#include "cgsolver.h"
#include "ba.h"


void conjGradSolver(double *hessian, double *g, int n, double* x) {
	static int oldN = 0;
	static double *p = NULL, *r = NULL, *d = NULL, *Bd = NULL, *pn = NULL;
	int j=0;
	double dBd, rr, alpha, beta, rrn, tmp;
	double error2; 
	double gradnorm2 = vecSqrLen(n,g);
	error2 = baMin(.5, pow(gradnorm2,.25));
	error2 *= error2*gradnorm2;
	error2 = baMin(1., error2);
	int numIterations = 0;

	if (n > oldN) {
		if (p) delete []p;
		p = new double[n];
		if (r) delete []r;
		r = new double[n];
		if (d) delete []d;
		d = new double[n];
		if (Bd) delete []Bd;
		Bd = new double[n];
		if (pn) delete []pn;
		pn = new double[n];
		oldN = n;
	}

	memset(p,0,n*sizeof(double));
	vecAssign(n, r, g);
	vecTimesScalar(n,r,-1.); 
	vecAssign(n, d, r);

	if ((tmp=vecSqrLen(n, r)) <= error2) {
		vecAssign(n,x,p);
//		printf("error low %f, inner converge\n",tmp);
		return;
	}

	rr = vecDot(n,r,r);  
	while (1) {
		j++;
		// Bd = H * d
		int hr, hc, index = 0;
		for (hr = 0; hr < n; hr++) {
			double sum = 0;
			for (hc = 0; hc < n; hc++) {
				sum += hessian[index++] * d[hc];
			}
			Bd[hr] = sum;
		}
		//B->matVecMult(d,Bd);
		dBd = vecDot(n,d,Bd);

		// case 1: negative curvature encountered, project to trust radius
		if (dBd < 0) {
			printf("negative curvature; aborting solve!\n");
			return;
		}

		alpha = rr / dBd;
		vecAssign(n,pn,d);
		vecTimesScalar(n,pn,alpha);
		vecAddEqual(n,pn,p);


		vecTimesScalar(n, Bd, alpha); // Bd now alpha*Bd
		vecDiffEqual(n,r,Bd);

		// Converged to full quadratic step
		double curErr = vecSqrLen(n, r);
		if (curErr < error2) {
			//      printf("gradient near nil, inner converge\n");
			vecAssign(n,x,pn);
			//      printf("final step, 2length %f\n",(tmp=vecSqrLen(n,pn)));
			//      printf("solved in %d  iterations \n",j);
			return;
		}
		else if (!_finite(curErr)) {
			printf("warning: error not finite!!\n");
			break;
		}

		rrn = vecDot(n,r,r);
		beta = rrn/rr;
		rr = rrn;
		vecTimesScalar(n,d,beta);
		vecAddEqual(n, d, r);
		vecAssign(n,p,pn);
		
		numIterations++;
	} 
}

