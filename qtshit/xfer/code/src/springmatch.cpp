#include "doppel2.h"
#include "springmatch.h"
#include "trimesh.h"
#include "solver.h"
#include "ba.h"
#include <vector>

using namespace std;

class SpringMatchGF : public IDifferentiableFunction {
public:
	TriMesh *mesh;
	int numVerts;
	vector<int> springStart, springEnd;
	vector<double> springLen;
	double lastErr;
	int innerStart;

	Vecd vars, grad, targetVerts;

	SpringMatchGF() {
		mesh = NULL;
		targetVerts = NULL;
		numVerts = 0;
		innerStart = 0;
	}

	void initSize(int iNumVerts) {
		numVerts = iNumVerts;
		springStart.clear();
		springEnd.clear();
		springLen.clear();

		vars.resize(numVerts * 3);
		grad.resize(numVerts * 3);
		targetVerts.resize(numVerts * 3);
	}

	void varsToMesh() {
		int i;
		for (i=0; i < numVerts; i++) {
			mesh->getPt(i) = Vec3d(vars[i*3+0], vars[i*3+1], vars[i*3+2]);
		}
		mesh->calcNormals();
		redrawV();
	}

	virtual double evaluateFunction(Vecd& variables) {
		lastErr = 0;
		grad.zeroElements();

		if (&variables != &vars)
			vars = variables;
		varsToMesh();

		int i;
		for (i = 0; i < springLen.size(); i++) {
			int v0 = springStart[i];
			int v1 = springEnd[i];
			Vec3d delta = 
				Vec3d(vars[3*v0+0], vars[3*v0+1], vars[3*v0+2]) -
				Vec3d(vars[3*v1+0], vars[3*v1+1], vars[3*v1+2]);

			double dist = delta.length();
			double err = dist - springLen[i];
			lastErr += sqr(err);

			grad[3*v0+0] += 2.0 * err * (1.0 / dist) * delta[0];
			grad[3*v0+1] += 2.0 * err * (1.0 / dist) * delta[1];
			grad[3*v0+2] += 2.0 * err * (1.0 / dist) * delta[2];
			grad[3*v1+0] -= 2.0 * err * (1.0 / dist) * delta[0];
			grad[3*v1+1] -= 2.0 * err * (1.0 / dist) * delta[1];
			grad[3*v1+2] -= 2.0 * err * (1.0 / dist) * delta[2];
		}

		for (i = 0; i < vars.size(); i++) {
//			lastErr += 0.001 * sqr(vars[i] - targetVerts[i]);
//			grad[i] += 0.001 * 2.0 * (vars[i] - targetVerts[i]);
		}

		uiWait();

		cout << "error: " << lastErr << endl;
		return lastErr;
	}

	virtual void evaluateGradient(Vecd& variables, Vecd& gradient) {
		gradient = grad;
	}

	virtual void solverStep() {
	}
};


SpringMatchGF smGF;
static LBFGSSolver *lbfgs = NULL;

void smInitMesh(TriMesh *mesh) {
	smGF.mesh = mesh;
	smGF.initSize(mesh->numPts());

	// surface springs
	int tri, edge;
	for (tri = 0; tri < mesh->numTris(); tri++) {
		for (edge = 0; edge < 3; edge++) {
			int v0 = mesh->getTri(tri, edge);
			int v1 = mesh->getTri(tri, (edge + 1) % 3);
			smGF.springStart.push_back(v0);
			smGF.springEnd.push_back(v1);
			smGF.springLen.push_back((mesh->getPt(v0) - mesh->getPt(v1)).length());
		}
	}

	// inner springs
	smGF.innerStart = (int)smGF.springLen.size();
/*	mesh->calcHBB(16);
	int pt;
	for (pt = 0; pt < mesh->numPts(); pt++) {
		mesh->calcRayIntersection(mesh->getPt(pt) + mesh->getPtNormal(pt) * 0.005, mesh->getPtNormal(pt));
		if (mesh->ptPos > 0) {
			double dist = (mesh->getPt(mesh->ptPos) - mesh->getPt(pt)).length();
			if (dist > 1e-6 && dist < 0.75) {
				smGF.springStart.push_back(pt);
				smGF.springEnd.push_back(mesh->ptPos);
				smGF.springLen.push_back(dist);
			}
		}
	}*/
}

void smInitTarget(TriMesh *mesh) {
	if (mesh->numPts() != smGF.numVerts) {
		cout << "WARNING: point size mismatch" << endl;
		return;
	}

	int pt;
	for (pt = 0; pt < mesh->numPts(); pt++) {
		Vec3d v = mesh->getPt(pt);
		smGF.targetVerts[pt * 3 + 0] = v[0];
		smGF.targetVerts[pt * 3 + 1] = v[1];
		smGF.targetVerts[pt * 3 + 2] = v[2];
	}
	smGF.vars = smGF.targetVerts;
}

void smRun(int maxIter) {
	if (lbfgs) {
		cout << "WARNING: solver already running" << endl;
		return;
	}
	if (!smGF.mesh) {
		cout << "WARNING: spring match not initialized" << endl;
		return;
	}

	lbfgs = new LBFGSSolver(&smGF);
	lbfgs->solve(1e+3, 1e-5, smGF.vars, maxIter);
	cout << "final error: " << smGF.lastErr << endl;

	delete lbfgs;
	lbfgs = NULL;
}

void smDrawInnerSprings() {
	if (!smGF.mesh)
		return;

	int i;
	glColor3f(0, 1, 0);
	glDisable(GL_LIGHTING);
	for (i = smGF.innerStart; i < smGF.springLen.size(); i++) {
		glBegin(GL_LINES);
		smGF.mesh->getPt(smGF.springStart[i]).glVertex();
		smGF.mesh->getPt(smGF.springEnd[i]).glVertex();
		glEnd();
	}
	glEnable(GL_LIGHTING);
}