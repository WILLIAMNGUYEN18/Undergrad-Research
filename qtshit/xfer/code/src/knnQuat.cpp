#include "knnQuat.h"

static const double EPSILON = 1e-10;

Vec3d knnQuatInterp(const QuatNorm &q, Vec4d *interpQuats, int numInterpQuats,
					Vec3d *samples, double *weights, Vec3d *derivs) {

	int i;

	if (numInterpQuats < 1) {
		for (i=0; i < 4; i++)
			derivs[i] = Vec3d();
		return Vec3d();
	}

	QuatNorm q2;
	const double DELTA = 1e-6;
	Vec3d ret = knnQuatInterp(q, interpQuats, numInterpQuats, samples, weights);

	// calculate derivatives using finite differences
	q2 = q;
	q2.x += DELTA;
	derivs[0] = (knnQuatInterp(q2, interpQuats, numInterpQuats, samples) - ret) / DELTA;
	q2 = q;
	q2.y += DELTA;
	derivs[1] = (knnQuatInterp(q2, interpQuats, numInterpQuats, samples) - ret) / DELTA;
	q2 = q;
	q2.z += DELTA;
	derivs[2] = (knnQuatInterp(q2, interpQuats, numInterpQuats, samples) - ret) / DELTA;
	q2 = q;
	q2.w += DELTA;
	derivs[3] = (knnQuatInterp(q2, interpQuats, numInterpQuats, samples) - ret) / DELTA;
	return ret;
}

// original knn method
Vec3d knnQuatInterp(const QuatNorm &q, Vec4d *interpQuats, int numInterpQuats,
					Vec3d *samples, double *weights) {

	// take care of the trivial case
	if (numInterpQuats < 1) {
		return Vec3d();
	}

	int i, j, k;
	double minDist[KNN_Q_K+1];
	int minIndex[KNN_Q_K+1];
	Vec3d ret;

	// clear the weights
	if (weights)
		memset(weights, 0, sizeof(double) * numInterpQuats);

	Vec4d curQ(q.x, q.y, q.z, q.w);
	double len = curQ.length();
	curQ = curQ / len;

	// if there are less than K samples, then use IDW
	if (numInterpQuats < KNN_Q_K+1) {
		double sum = quatDist(curQ, Vec4d(0, 0, 0, 1));
		if (sum < EPSILON) {
			// we're at the zero quat!
			return Vec3d();
		}
		sum = 1.0 / sum;
		for (i=0; i < numInterpQuats; i++) {
			double dist = quatDist(curQ, interpQuats[i]);
			if (dist < EPSILON) {
				// we're at one if the key quats
				if (weights) weights[i] = 1;
				return samples[i];
			}
			sum += 1.0/dist;
		}
		for (i=0; i < numInterpQuats; i++) {
			double w = (1.0 / quatDist(curQ, interpQuats[i])) / sum;
			if (weights) weights[i] = w;
			ret += w * samples[i];
		}

		return ret;
	}

	// first, find the top K+1 candidates and their distance
	for (i=0; i < KNN_Q_K+1; i++) {
		minDist[i] = DBL_MAX;
		minIndex[i] = -2;
	}
	minDist[0] = quatDist(curQ, Vec4d(0, 0, 0, 1));
	minIndex[0] = -1;
	for (i=0; i < numInterpQuats; i++) {
		double curDist = quatDist(curQ, interpQuats[i]);
		int curIndex = i;

		for (j=0; j < KNN_Q_K+1; j++) {
			if (curDist < minDist[j]) {
				// push old values down
				for (k = KNN_Q_K; k > j; k--) {
					minDist[k] = minDist[k-1];
					minIndex[k] = minIndex[k-1];
				}
				// insert new value
				minDist[j] = curDist;
				minIndex[j] = curIndex;
				break;
			}
		}
	}

	// now do the interpolation
	if (minDist[0] < EPSILON) {
		// if we're very close to one point, just use it
		if (minIndex[0] < 0)
			ret = Vec3d();
		else {
			ret = samples[minIndex[0]];
			if (weights) weights[minIndex[0]] = 1;
		}
	}
	else {
		double cur, sum = 0;
		int last = KNN_Q_K;
		while (minIndex[last] == -2 && last > 0)
			last--;
		double thresh;
		thresh = minDist[last];

		if (fabs(minDist[0] - thresh) < EPSILON) {
			// in the unlikely case that all weights are equal, then 
			// just evenly interpolate
			for (i=0; i < last; i++) {
				if (minIndex[i] >= 0) {
					ret += (1.0 / last) * samples[minIndex[i]];
					if (weights) weights[minIndex[i]] = 1.0 / last;
				}
			}
		}
		else {
			for (i=0; i < last; i++) {
				cur = (thresh / minDist[i]) - 1.0;
				sum += cur;
			}
			for (i=0; i < last; i++) {
				cur = (thresh / minDist[i]) - 1.0;
				if (minIndex[i] >= 0) {
					if (weights) weights[minIndex[i]] = cur/sum;
					ret += (cur / sum) * samples[minIndex[i]];
				}
			}
		}
	}
	return ret;
}

