#ifndef RBF_H
#define RBF_H

#include "quatnorm.h"
#include <float.h>
#include "vl/VLd.h"


class RBF { 
public:
	int mN, mD;
	double *mRadius;
	bool mQuatDist;

	VLMatd samplePts;
	VLMatd radialBasisWeights;

	double *curWeights, *curDerivs;

	RBF(int n, int d, bool quatDist);
	~RBF();
	void init(double *samples, bool firstIsZero = false);

	double radialEval(int comp, double *pt);
	double radialDist(double *pt1, double *pt2, double radius);
	double radialDist(double dist, double radius);
	
	void eval(double *pt, double *weights = NULL);
	void evalDerivs(double *pt);
};


Vec3d rbfQuatInterp(const QuatNorm &q, Vec4d *interpQuats, int numInterpQuats,
				Vec3d *samples, double *weights = NULL);

Vec3d rbfQuatInterp(const QuatNorm &q, Vec4d *interpQuats, int numInterpQuats,
				Vec3d *samples, double *weights, Vec3d *derivs);

#endif
