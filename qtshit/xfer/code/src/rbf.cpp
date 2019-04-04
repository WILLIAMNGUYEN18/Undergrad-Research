#include "rbf.h"

static const double EPSILON = 1e-10;


RBF::RBF(int n, int d, bool quatDist) {
	int i, j;

	mN = n;
	mD = d;
	mRadius = new double[n];
	for (i=0; i < n; i++)
		mRadius[i] = 1;
	mQuatDist = quatDist;

	samplePts.SetSize(mN, mD+1);
	radialBasisWeights.SetSize(mN, mN);

	curWeights = new double[mN];
	curDerivs = new double[mN * mD];
}

RBF::~RBF() {
	if (curWeights)
		delete []curWeights;
	if (curDerivs)
		delete []curDerivs;
}

void RBF::init(double *samples, bool firstIsZero) {
	int i, j;

	// initialize sample points
	samplePts = vl_zero;
	if (firstIsZero) {
		samplePts[0][mD] = 1;
		if (mQuatDist)
			samplePts[0][mD-1] = 1;

		for (i=1; i < mN; i++) {
			double sum = 0;

			for (j=0; j < mD; j++) {
				samplePts[i][j] = samples[(i-1)*mD + j];
				sum += sqr(samplePts[i][j]);
			}

			if (mQuatDist && sum > 0) {
				sum = sqrt(sum);
				for (j=0; j < mD; j++) {
					samplePts[i][j] /= sum;
				}
			}

			samplePts[i][j] = 1;
		}
	}
	else {
		for (i=0; i < mN; i++) {
			double sum = 0;

			for (j=0; j < mD; j++) {
				samplePts[i][j] = samples[i*mD + j];
				sum += sqr(samplePts[i][j]);
			}

			if (mQuatDist && sum > 0) {
				sum = sqrt(sum);
				for (j=0; j < mD; j++) {
					samplePts[i][j] /= sum;
				}
			}

			samplePts[i][j] = 1;
		}
	}

	for (i=0; i < mN; i++) {
		mRadius[i] = 1e6;
		for (j=0; j < mN; j++) {
			if (j == i)
				continue;
			double dist = 0;
			if (mQuatDist) {
				dist = quatDist(samplePts[i].Ref(), samplePts[j].Ref());
			}
			else {
				int k;
				for (k = 0; k < mD; k++)
					dist += sqr(samplePts[i][k] - samplePts[j][k]);
				dist = sqrt(dist);
			}
			if (dist * 2.0 < mRadius[i])
				mRadius[i] = dist * 2.0;
		}
	}

	// calculate radial basis weights
	for (i=0; i < mN; i++) {
		for (j=0; j < mN; j++) {
			radialBasisWeights[i][j] = radialDist(
				samplePts.Ref() + i*(mD+1), 
				samplePts.Ref() + j*(mD+1), mRadius[i]);
		}
	}
//	cout << radialBasisWeights << endl;
	radialBasisWeights = inv(radialBasisWeights);
}

double RBF::radialEval(int comp, double *pt) {
	double ret = 0;
	int i;

	if (mQuatDist) {
		double sum = 0;
		for (i=0; i < mD; i++)
			sum += sqr(pt[i]);
		if (sum > 0) {
			sum = sqrt(sum);
			for (i=0; i < mD; i++)
				pt[i] /= sum;
		}
	}

	for (i=0; i < mN; i++) {
		ret += radialBasisWeights[comp][i] * 
			radialDist(samplePts.Ref() + i*(mD+1), pt, mRadius[i]);
	}
	return ret;
}

double RBF::radialDist(double *pt1, double *pt2, double radius) {
	int i;
	double dist = 0;

	if (mQuatDist) {
		dist = quatDist(pt1, pt2);
	}
	else {
		for (i = 0; i < mD; i++)
			dist += sqr(pt1[i] - pt2[i]);
		dist = sqrt(dist);
	}
	return radialDist(dist, radius);
}

double RBF::radialDist(double dist, double radius) {
	dist = dist/(radius / 2.0);

	if (dist <= 1.0) {
		dist = 1.0 - dist;
		return -dist*dist*dist/2.0 + dist*dist/2.0 + dist/2.0 + 1.0/6.0;
	}
	else if (dist <= 2.0) {
		dist = 2.0 - dist;
		return dist * dist * dist / 6.0;
	}
	else 
		return 0;
}


void RBF::eval(double *point, double *result) {
	int i, j;

	if (result == NULL)
		result = curWeights;

	double sum = 0;
	for (i=0; i < mN; i++) {
		result[i] = radialEval(i, point);
		sum += result[i];
	}

	if (sum != 0) {
		for (i=0; i < mN; i++) {
			result[i] /= sum;
		}
	}
}

void RBF::evalDerivs(double *point) {
	int i, j;
	const static double DELTA = 0.00001;

	// store current weights
	eval(point, curWeights);

	// estimate derivatives by finite differences
	for (i=0; i < mD; i++) {
		double oldV = point[i];
		point[i] += DELTA;
		eval(point, curDerivs + i*mN);
		point[i] = oldV;

		for (j=0; j < mN; j++) {
			curDerivs[i*mN + j] = (curDerivs[i*mN + j] - curWeights[j]) / DELTA;
		}
	}
}

/*

linear basis stuff

	VLMatd linearBasisWeights;

in init:

	// calculate linear basis using pseudoinverse
	VLMatd U(n, d+1), V(d+1, d+1);
	VLVecd diagonal(d+1);
	SVDFactorization(samplePts, U, V, diagonal);
	linearBasisWeights.SetSize(d+1, n);
	VLMatd diag(d+1, d+1);
	diag = vl_0;
	for (i = 0; i < d+1; i++) {
		if (fabs(diagonal[i]) > 0.00001)
			diag[i][i] = 1.0 / diagonal[i];
	}
	linearBasisWeights = V * diag * trans(U);

	// calculate residuals
	VLMatd residuals(n, n);
	for (i=0; i < n; i++) {
		for (j=0; j < n; j++) {
			int kronDelta = (i == j)?1:0;
			residuals[i][j] = kronDelta - linearPart(j, samplePts.Ref()+i*(d+1));
		}
	}


	radialBasisWeights.SetSize(n, n);


double RBF::linearPart(int comp, double *pt) {
	double ret;
	int i;
	ret = linearBasisWeights[comp][mD];
	for (i=0; i < mD; i++) {
		ret += linearBasisWeights[i][comp] * pt[i];
	}
	return ret;
}

*/
