#include <iostream>
#include "jactest.h"
#include "vec.h"

using namespace std;

void runTest(IDifferentiableFunction *func, Vecd &data, int minV, int maxV) {
	double val, newVal, newEst;
	int i, n = data.size();
	Vecd grad(n), orig(n);

	for (i=0; i < n; i++) {
		orig[i] = data[i];
	}
	
	val = func->evaluateFunction(data);
	func->evaluateGradient(data, grad);

	double max = 0;
	int maxInd = 0;
	for (i=0; i < n; i++) {
		if (fabs(grad[i]) > fabs(max)) {
			max = grad[i];
			maxInd = i;
		}
	}
	cout << "derivative magnitude: " << grad.length() << "; maximum element: " << max << " (" << maxInd << ")" << endl;

	int MOD_DOF;
	double DELTA;
	double err;

	if (minV < 0)
		minV = 0;
	if (maxV < 0 || maxV > n)
		maxV = n;

	for (MOD_DOF = minV; MOD_DOF < maxV; MOD_DOF++) {
//		if (data[MOD_DOF] == 0) {
//			continue;
//		}

		//cout << "testing degree of freedom #" << MOD_DOF << endl;
		cout << MOD_DOF << ": ";
		// now tweak the dofs
		DELTA = 0.0001;

		data[MOD_DOF] += DELTA;

		newVal = func->evaluateFunction(data);
		newEst = val + grad[MOD_DOF] * DELTA;

		//err = ((newVal - newEst) / DELTA);
		//cout << err;
		//cout << grad[MOD_DOF] << endl;
		//cout << (newVal - val) << endl;

		err = newVal - newEst;
//		err = (newVal - val) / DELTA - grad[MOD_DOF];
		cout << "err = " << err << "; est. delta = " << (grad[MOD_DOF] * DELTA) << "; true delta = " << (newVal - val) << endl;

		data[MOD_DOF] = orig[MOD_DOF]; //-= DELTA;

		cout << MOD_DOF << ": ";
		DELTA = 0.00001;
		data[MOD_DOF] += DELTA;
		newVal = func->evaluateFunction(data);
		newEst = val + grad[MOD_DOF] * DELTA;
		err = newVal - newEst;
		cout << "err = " << err << "; est. delta = " << (grad[MOD_DOF] * DELTA) << "; true delta = " << (newVal - val) << endl;
		data[MOD_DOF] = orig[MOD_DOF]; //-= DELTA;
		
		/*

		DELTA = 0.0001;
		dofs[MOD_DOF] += DELTA;

		newVal = func->evaluateFunction(dofs);
		newEst = val + grad[MOD_DOF] * DELTA;

		err2 = ((newVal - newEst) / DELTA);
		cout << " -- " << err2 << " -- " << (err / err2);
		if (fabs((err / err2) - 10.0) > 0.1)
			cout << " ******";
		cout << endl;
		
		dofs[MOD_DOF] -= DELTA;*/
	}

}
