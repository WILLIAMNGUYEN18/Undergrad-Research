#ifndef PART_MATCH_H
#define PART_MATCH_H

#include "vec.h"
#include "solver.h"

class PCAPart;
class TriMesh;
class MarkerSet;

class PartMatchGF : public IDifferentiableFunction {
public:
	double lastErr;
	double markerVariance;

	Vecd vars, grad;

	PCAPart *part;
	TriMesh *targetMesh;
	MarkerSet *targetMkr;

	bool optTransform, optPCA;

	PartMatchGF(PCAPart *iPart, TriMesh *iTargetMesh, MarkerSet *iTargetMkr);
	~PartMatchGF();

	void applyDef(Vecd &variables, bool updateAll = true);
	void zeroDeformation(bool recalc = true);

	void updateCurrent(int pt);
	void addGradientVec(int index, Vec3d v);

	virtual double evaluateFunction(Vecd& variables);
	virtual void evaluateGradient(Vecd& variables, Vecd& gradient);
//	virtual void solverStep();
};

#endif