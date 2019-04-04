#include <float.h>
#include "doppel2.h"
#include "cvgf.h"
#include "skeleton.h"
#include "trimesh.h"
#include "viewer.h"
#include "jactest.h"
#include "skinCalc.h"
#include "cgsolver.h"
#include "vl/VLd.h"


const bool useConf = true;
bool verboseSolver = false;


// CVExample methods ================================================

void CVExample::init(int nPoints) {
	numPoints = nPoints;
	numTPoints = nPoints;
	trans = new Mat4d[numTPoints];
	iTrans = new Mat4d[numTPoints];
	points = new Vec3d[numPoints];
	conf = new double[numPoints];
	numDofs = 42;
	dofs = new double[numDofs];
}

void CVExample::buildExample(TriMesh *tm, Skin *skin) {
	int i, targetPt = 0;

	numPoints = tm->numPts();
	points = new Vec3d[numPoints];
	conf = new double[numPoints];
	for (i=0; i < numPoints; i++) {
		points[targetPt] = tm->getPt(i);
		conf[targetPt] = 1;

		// cut off the platform
//		if (points[targetPt][2] > -0.980)
			targetPt++;
	}
	cout << "eliminated " << (numPoints-targetPt) << " points" << endl;
	numPoints = targetPt;

	numTPoints = skin->numPts;
	trans = new Mat4d[numTPoints];
	iTrans = new Mat4d[numTPoints];
	for (i=0; i < numTPoints; i++) {
#ifdef USE_POSITIONS
		Mat4d m = skin->points[i].localMat;
#else
		Mat4d m;
#endif
		trans[i] = m;
		iTrans[i] = m.inverse();
	}

	CharVec::dofsFromSkel(skin->skel, dofs);
}

void CVExample::save(char *fname) {
	FILE *f;

	if (!openFile(&f, fname, "wb", "example file"))
		return;

	fwrite(&numPoints, sizeof(int), 1, f);
	fwrite(points, sizeof(Vec3d), numPoints, f);
	fwrite(conf, sizeof(double), numPoints, f);
	fwrite(&numTPoints, sizeof(int), 1, f);
	fwrite(trans, sizeof(Mat4d), numTPoints, f);
	fwrite(iTrans, sizeof(Mat4d), numTPoints, f);
	fwrite(&numDofs, sizeof(int), 1, f);
	fwrite(dofs, sizeof(double), numDofs, f);
//	fwrite(&charID, sizeof(int), 1, f);

	fclose(f);
}

void CVExample::load(char *fname) {
	FILE *f;

	if (!openFile(&f, fname, "rb", "example file"))
		return;

	fread(&numPoints, sizeof(int), 1, f);
	points = new Vec3d[numPoints];
	fread(points, sizeof(Vec3d), numPoints, f);
	conf = new double[numPoints];
	fread(conf, sizeof(double), numPoints, f);
	fread(&numTPoints, sizeof(int), 1, f);
	trans = new Mat4d[numTPoints];
	fread(trans, sizeof(Mat4d), numTPoints, f);
	iTrans = new Mat4d[numTPoints];
	fread(iTrans, sizeof(Mat4d), numTPoints, f);
	fread(&numDofs, sizeof(int), 1, f);
	dofs = new double[numDofs];
	fread(dofs, sizeof(double), numDofs, f);
//	fread(&charID, sizeof(int), 1, f);

	fclose(f);
}


// CVGoalFunction methods ===========================================

CVGoalFunction::CVGoalFunction() {
	cv = NULL;
	solver = NULL;

	distDeviation = 0.02;
	matchDeviation = 0.05;

	lockShape = false;
	lockPose = false;

	numExamples = 0;

	useRegularization = true;
}

void CVGoalFunction::init(int iNumExamples, int iNumChars) {
	int i;

	numExamples = iNumExamples;
	numChars = iNumChars;

	cout << "examples: " << numExamples << endl;
	cout << "components: " << cv->numComponents << endl;
	cout << "characters: " << numChars << endl;

	// initialize mu and phi
	curMu = new VLVecd[numChars];
	curPhi = new VLMatd[numChars];
	for (i=0; i < numChars; i++) {
		curMu[i].SetSize(cv->numComponents);
		curMu[i].MakeZero();
		curMu[i][0] = 1;

		curPhi[i].SetSize(cv->numComponents, cv->numComponents);
		curPhi[i].MakeZero();
		curPhi[i][0][0] = 1;
	}
}

int CVGoalFunction::numObsPts() {
	int ex, ret = 0;
	for (ex = 0; ex < numExamples; ex++) {
		ret += examples[ex].numPoints;
	}
	return ret;
}

void CVGoalFunction::runEM(int numIterations) {
	int i, j;

	if (solver != NULL) {
		cout << "solver already running! "  << endl;
		return;
	}

	solver = (LBFGSSolver*)1;

	lockPose = true;
	lockShape = false;

	calcFreeEnergy(true);

	/*
	// save W
	ofstream out("W.txt");
	int k;
	for (j=0; j < cv->numPts; j++) {
		for (k=0; k < 3; k++) {
			for (i=1; i < cv->numComponents; i++) {
				if (i > 1)
					 out << ", ";
				out << cv->cvPts[j].curComponents[i][k];
			}
			out << ";" << endl;
		}
	}
	out.close();
	// save b
	out.open("b.txt");
	for (j=0; j < cv->numPts; j++) {
		for (k=0; k < 3; k++) {
			out << cv->cvPts[j].data[k] << ";" << endl;
		}
	}
	out.close();
	// save y0
	out.open("y0.txt");
	for (j=0; j < cv->numPts; j++) {
		for (k=0; k < 3; k++) {
			out << examples[0].points[j][k] << ";" << endl;
		}
	}
	out.close();
	*/

	for (i=0; i < numIterations; i++) {
		useCovariance = true;
//	distDeviation = 0.002;
//	matchDeviation = 0.05;

		cout << "-------------------- E step " << i << " --------------------" << endl;
		baTimer();
		eStep();
		cout << "E-step time: " << baTimer() << endl;
		calcFreeEnergy(true);

//		cout << "normalizing..." << endl;
//		normalizeN();

		cout << "-------------------- M step " << i << " --------------------" << endl;
		baTimer();
//		if (i > 0)
//			optimizeWb(true);
//		verboseSolver = (i > 1);
/*		if (i < 2)
			useCovariance = false;
		else
			useCovariance = true;
*/		optimizeWb();
//		if (i < numIterations-1) {
//			cout << "regularizing..." << endl;
//			for (j=0; j < 5; j++)
//				regularize();
//		}
		cout << "M-step time: " << baTimer() << endl;

		// update UI
		calcFreeEnergy(true);
//		cv->updateCurComponents(examples[numExamples-1].dofs);
//		cv->updateLocalPos(curMu[numChars-1].Ref());
		cv->updateTM(examples[numExamples-1].trans);
		redrawVNow();
		uiWait();

//		if (solver->stopNow == true)
//			break;
	}

	solver = NULL;
}

void CVGoalFunction::eStep() {
	// i = scan (example) index
	// j = index of a point in scan i
	// jComp = component (among x, y, z)
	int i, j, jComp;
	// y_ij = vertex with index j from example i
	Vec3d y_ij;
	// ii, jj, kk are used to iterate over vectors and matrices
	int ii, jj, kk;
	double sigma2 = distDeviation * distDeviation;
	// current and last character id
	int charID = -1, lastID = -1;

	// static data
	static bool doneInit = false;
	static VLVecd muAcc;


	// initialize static data
	if (!doneInit) {
		muAcc.SetSize(cv->numComponents);
		doneInit = true;
	}


	charID = examples[0].charID;
	for (i = 0; i < numExamples; i++) {
		// calculate the current reconstruction
		cv->updateCurComponents(examples[i].dofs);
		cv->updateLocalPos(curMu[charID].Ref());
//		cv->updateTM(examples[i].trans);

		if (lastID != charID) {
			// clear the current average and covariance estimates
			curMu[charID].MakeZero();
			curMu[charID][0] = 1;
			curPhi[charID].MakeDiag(sigma2);
			curPhi[charID][0][0] = 1;
			
			muAcc.MakeZero();
		}

		if (cv->numComponents < 2)
			continue;

		// iterate over all observed points
		for (j = 0; j < examples[i].numPoints; j++) {
			y_ij = examples[i].points[j];

			for (jComp = 0; jComp < 3; jComp++) {
				// accumulate WffW into phi
				for (ii=1; ii < cv->numComponents; ii++) {
					for (jj=1; jj < cv->numComponents; jj++) {
						curPhi[charID][ii][jj] += cv->cvPts[j].curComponents[ii][jComp] * cv->cvPts[j].curComponents[jj][jComp];
					}
				}

				// accumulate (y-bf)f'W into muAcc
				for (ii=1; ii < cv->numComponents; ii++) {
					muAcc[ii] += (y_ij[jComp] - cv->cvPts[j].curComponents[0][jComp]) * cv->cvPts[j].curComponents[ii][jComp];
				}
			}
		}

		lastID = charID;
		if (i < numExamples-1) {
			charID = examples[i+1].charID;
		}
		if ((i == numExamples-1) || (lastID != charID)) {

//			if (lastID == 0) {
//				cout << "phi" << endl << curPhi[lastID];
//			}

			// actual phi uses the inverse
			curPhi[lastID] = inv(curPhi[lastID]);

//			if (lastID == 0) {
//				cout << "inv(phi)" << endl << curPhi[lastID];
//			}

			// mu = muAcc * phi
//			if (lastID == 0)
//				cout << "component values for character " << lastID << ": " << endl;
			for (ii=1; ii < cv->numComponents; ii++) {
				curMu[lastID][ii] = 0;
				for (jj = 1; jj < cv->numComponents; jj++) {
					curMu[lastID][ii] += (muAcc[jj] * curPhi[lastID][jj][ii]);
				}
//				cout << curN[lastID][ii] << " ";
				curPhi[lastID][0][ii] = curMu[lastID][ii];
				curPhi[lastID][ii][0] = curMu[lastID][ii];
			}
			curPhi[lastID][0][0] = 1;

			// add mu'mu into curPhi
			for (ii=1; ii < cv->numComponents; ii++) {
				for (jj=1; jj < cv->numComponents; jj++) {
					curPhi[lastID][ii][jj] = sigma2 * curPhi[lastID][ii][jj] + curMu[lastID][ii] * curMu[lastID][jj];
				}
			}

//			if (lastID == 0) {
//				cout << "E(n^2)" << endl << curPhi[lastID];
//			}

//			if (lastID == 0) {
//				cout << curPhi[lastID];
//				cout << endl;
//			}
			cout << curMu[lastID] << endl;
		}

//		delete tree;
	}
}

void CVGoalFunction::optimizeWb(bool test)  {
	int i;
	double distF = 1.0 / (2.0 * sqr(distDeviation));
	double matchF = 1.0 / (2.0 * sqr(matchDeviation));


	double totalDErr = 0;
	double totalMErr = 0;
	double totalPts = 0;
	int totalNeighbors = 0;
	cv->backupData();
	solvingPtIndex = 0;
	for (solvingPt = 0; solvingPt < cv->numPts; solvingPt++) {
		if (solveGrad.size() < cv->cvPts[solvingPt].data.size())
			solveGrad.resize(cv->cvPts[solvingPt].data.size());

		solver = new LBFGSSolver(this);
		if (test && solvingPt == 0) {
			runTest(this, cv->cvPts[solvingPt].data);
			exit(0);
		}
		solver->solve(1e-7, 1e-5, cv->cvPts[solvingPt].data, 20);
		totalDErr += curDErr;
		totalMErr += curMErr;

		
		if (useConf)
			for (i=0; i < numExamples; i++)
				totalPts += examples[i].conf[solvingPt];
		else
			for (i=0; i < numExamples; i++)
				totalPts += 1;
		totalNeighbors += (int)cv->neighbors[solvingPt].size();
		delete solver;
		if (verboseSolver)
			cout << endl << endl;

		solvingPtIndex += cv->cvPts[solvingPt].data.size();

		if (solvingPt % 1000 == 0)
			cout << ".";
	}
	cout << endl;

	cout << "variational free energy: " << totalDErr << " + " << totalMErr << " = " << 
		(totalDErr + totalMErr) << endl; // + 
//		 1.5*totalPts*log(sqr(distDeviation)) + 
//		 1.5*totalNeighbors*log(sqr(matchDeviation))) << endl;

	// autoupdate
//	distDeviation = sqrt((2.0 * totalDErr) / (3.0 * totalPts * MAX_MM_PTS));
//	matchDeviation = sqrt((2.0 * totalMErr) / (3.0 * totalNeighbors));
//	cout << "sigma_d = " << distDeviation << "; sigma_m = " << matchDeviation << endl;

	/*
	// for each charVec point, move to the optimal position
	cout << "optimizing..." << endl;
	int numFuncs, maxFuncs = maxMatches + 10 * (numComponents+1);
	double *jac = new double[(numComponents+1)*maxFuncs];
	double *hessian = new double[(numComponents+1)*(numComponents+1)];
	double *gradient = new double[maxFuncs];
	double *vars = new double[numComponents+1];
	double *residual = new double[maxFuncs];
	curGrad = curVars;  // we'll co-opt curGrad to store the new values
	for (cPt = 0; cPt < cv->numPts; cPt++) {
		Vec3d *color = &cv->tm->getPtColor(cPt);
		*color = Vec3d();

		int numInf = matches[cPt].size();
		if (numInf < 1) {
//			continue;
			*color = Vec3d(1, 0, 1);
		}
		else
			*color = Vec3d(0.8, 0.8, 0.8);
//		*color = Vec3d(0,0,0);

		for (xyz=0; xyz < 3; xyz++) {
			// calculate jacobian
			memset(jac, 0, sizeof(double)*maxFuncs*(numComponents+1));
			numFuncs = 0;
			for (i = 0; i < numInf; i++) {
				MMRecord *curMatch = matches[cPt][i];

				for (comp = 0; comp < numComponents+1; comp++) {
					if (comp == 0)
						jac[numFuncs*(numComponents+1) + comp] = -curMatch->weight;
					else
						jac[numFuncs*(numComponents+1) + comp] = -curMatch->weight * curN[curMatch->ex*numComponents + comp-1];
				}
				residual[numFuncs] = matches[cPt][i]->err[xyz];
//				(*color)[xyz] += min(1.0, sqr(matches[cPt][i]->err[xyz]));

				numFuncs++;
			}
			for (i = 0; i < cv->neighbors[cPt].size(); i++) {
				int neigh = cv->neighbors[cPt][i];

				for (comp = 0; comp < numComponents+1; comp++) {
					jac[numFuncs*(numComponents+1) + comp] = -matchF;

					double dist;
					int nInd = getIndex(comp-1, neigh, 0)+xyz;
					int ind = getIndex(comp-1, cPt, 0)+xyz;
					if (comp > 0)
						dist = (curVars[nInd] - curVars[ind]);
					else {
						dist = ((curVars[nInd] - templateVars[nInd]) - (curVars[ind] - templateVars[ind]));
//						(*color)[xyz] = min(1.0, sqr(matchF * dist));
					}
					residual[numFuncs] = matchF * dist;

					numFuncs++;
				}
			}

			// hessian = jac' * jac; gradient = jac' * residual
			for (comp = 0; comp < numComponents+1; comp++) {
				double sum = 0;

				for (i = 0; i < numFuncs; i++) {
					sum += jac[i*(numComponents+1) + comp] * residual[i];
				}
				gradient[comp] = sum;

				int comp2;
				for (comp2 = 0; comp2 < numComponents+1; comp2++) {
					sum = 0;
					for (i = 0; i < numFuncs; i++) {
						sum += jac[i*(numComponents+1) + comp] * jac[i*(numComponents+1) + comp2];
					}
					hessian[comp*(numComponents+1) + comp2] = sum;
				}
			}

			conjGradSolver(hessian, gradient, (numComponents+1), vars);
			// add delta into solution
			for (i = 0; i < numComponents+1; i++) {
				int index = getIndex(i-1, cPt, 0) + xyz;
				curGrad[index] += vars[i];
//				(*color)[xyz] += min(vars[i] / 0.01, 1.0);
			}

//			cv->tm->getPtColor(cPt) = Vec3d(1, 0, 1);
		}
	}
	curVars = curGrad;
	delete []gradient;
	delete []hessian;
	delete []jac;
	delete []vars;
	delete []residual;
	*/
}

double CVGoalFunction::evaluateFunction(Vecd& variables) {
	int ex, index, comp, inf, degree, i;
	int ii, jj;
	double trace;
	double distF = 1.0 / (2.0 * sqr(distDeviation));
	double matchF = 1.0 / (2.0 * sqr(matchDeviation));
	int solvingXYZ;

	CharVecPt *cvPt = &(cv->cvPts[solvingPt]);

	if (&variables != &(cv->cvPts[solvingPt].data))
		cv->cvPts[solvingPt].data = variables;

	curDErr = 0;
	curMErr = 0;

	solveGrad.zeroElements();
	for (ex = 0; ex < numExamples; ex++) {
		double *mult = cvPt->calcMultipliers(examples[ex].dofs);
		cvPt->updateCurComponents();
		cvPt->updateLocalPos(curMu[examples[ex].charID].Ref());
		
		Vec3d &target = examples[ex].points[solvingPt];

		for (solvingXYZ = 0; solvingXYZ < 3; solvingXYZ++) {
			if (cv->numComponents < 2 || !useCovariance) {
				double err = target[solvingXYZ] - cvPt->localPos[solvingXYZ];
				double weight = distF;
				if (useConf)
					weight *= examples[ex].conf[solvingPt];
				curDErr += weight * err * err;


				// calc derivs of err * err
				index = 0;
				for (inf = 0; inf < cvPt->dataParts; inf++) {
					for (comp = 0; comp < cv->numComponents; comp++) {
						if (comp == 0)
							solveGrad[index+solvingXYZ] += -2.0 * err * mult[inf] * weight;
						else
							solveGrad[index+solvingXYZ] += -2.0 * err * mult[inf] * curMu[examples[ex].charID][comp] * weight;
						index += 3;
					}
				}
			}
			else {
				double trace = 0;
				double weight = distF;
				if (useConf)
					weight *= examples[ex].conf[solvingPt];

				// calculate tr(W'ff'W phi)
				for (ii=0; ii < cv->numComponents; ii++) {
					for (jj=0; jj < cv->numComponents; jj++) {
						trace += 
							cvPt->curComponents[ii][solvingXYZ] * 
							cvPt->curComponents[jj][solvingXYZ] *
							curPhi[examples[ex].charID][jj][ii];

						// update gradients for stuff that made curComponents[ii] and curComponents[jj]
						double iMult = cvPt->curComponents[jj][solvingXYZ] *
							curPhi[examples[ex].charID][jj][ii];
						double jMult = cvPt->curComponents[ii][solvingXYZ] *
							curPhi[examples[ex].charID][jj][ii];
						index = 0;
						for (inf = 0; inf < cvPt->dataParts; inf++) {
							solveGrad[index+(3*ii)+solvingXYZ] += mult[inf] * iMult * weight;
							solveGrad[index+(3*jj)+solvingXYZ] += mult[inf] * jMult * weight; // iMult?

							index += cv->numComponents*3;
						}
						baAssert(index == cvPt->data.size());
					}
				}

				// calc derivs of - 2.0 * target[solvingXYZ] * cvPt->localPos[solvingXYZ]
				index = 0;
				for (inf = 0; inf < cvPt->dataParts; inf++) {
					for (comp = 0; comp < cv->numComponents; comp++) {
						if (comp == 0)
							solveGrad[index+solvingXYZ] += -2.0 * target[solvingXYZ] * mult[inf] * weight;
						else
							solveGrad[index+solvingXYZ] += -2.0 * target[solvingXYZ] * mult[inf] * curMu[examples[ex].charID][comp] * weight;
						index += 3;
					}
				}

/*cout << 					(
					target[solvingXYZ] * target[solvingXYZ]
					- 2.0 * target[solvingXYZ] * cvPt->localPos[solvingXYZ]
					+ trace) << " " << sqr(target[solvingXYZ] - cvPt->localPos[solvingXYZ]) << endl;*/

				curDErr += 
					weight * (
					target[solvingXYZ] * target[solvingXYZ]
					- 2.0 * target[solvingXYZ] * cvPt->localPos[solvingXYZ]
					+ trace);

				if (!_finite(curDErr)) {
					cout << "data error went infinite!!" << endl;
				}
			}
		}
	}

	
	// regularization
	if (useRegularization) {
		int varsPerInf = cv->numComponents*(cvPt->dataParts / cvPt->numInfluences)*3;
		int neigh, j;
		for (inf = 1; inf < cvPt->numInfluences; inf++) {
			int index = inf*varsPerInf;

			int numNeigh = (int)cv->neighbors[solvingPt].size();
			for (neigh = 0; neigh < numNeigh; neigh++) {
				double neighVal = 0;
				int neighIndex = -1;
				int neighInf = cv->cvPts[cv->neighbors[solvingPt][neigh]].getInfIndex(cv->cvPts[solvingPt].influences[inf]);
				if (neighInf >= 0) {
					neighIndex = neighInf * varsPerInf;

					for (j = 0; j < varsPerInf; j++) {
						double dist = cv->cvPts[solvingPt].data[inf*varsPerInf + j] - cv->cvPts[cv->neighbors[solvingPt][neigh]].backup[neighIndex + j];
						curMErr += matchF * sqr(dist);
						solveGrad[inf*varsPerInf + j] += 2.0 * matchF * dist;
					}
				}
				else {
					for (j = 0; j < varsPerInf; j++) {
						double dist = cv->cvPts[solvingPt].data[inf*varsPerInf + j];
						curMErr += matchF * sqr(dist);
						solveGrad[inf*varsPerInf + j] += 2.0 * matchF * dist;
					}
				}
			}
		}
	}

	curErr = curDErr + curMErr;
	if (verboseSolver)
		cout << curErr << " ";

	return curErr;
}

void CVGoalFunction::evaluateGradient(Vecd& variables, Vecd& gradient) {
	memcpy(gradient.n, solveGrad.n, sizeof(double)*gradient.size());
}

void CVGoalFunction::regularize() {
	int solvingPt, inf, j, neigh;

	int varsPerInf = cv->numComponents*(cv->cvPts[0].dataParts / cv->cvPts[0].numInfluences)*3;

	cv->backupData();

	for (solvingPt=0; solvingPt < cv->numPts; solvingPt++) {

		int numNeigh = (int)cv->neighbors[solvingPt].size();
		for (j = varsPerInf; j < cv->cvPts[solvingPt].data.size(); j++) {
			cv->cvPts[solvingPt].data[j] /= (numNeigh+1.0);
		}

		for (inf = 1; inf < cv->cvPts[solvingPt].numInfluences; inf++) {
			int index = inf*varsPerInf;

			int numNeigh = (int)cv->neighbors[solvingPt].size();
			for (neigh = 0; neigh < numNeigh; neigh++) {
				double neighVal = 0;
				int neighIndex = -1;
				int neighInf = cv->cvPts[cv->neighbors[solvingPt][neigh]].getInfIndex(cv->cvPts[solvingPt].influences[inf]);
				if (neighInf >= 0) {
					neighIndex = neighInf * varsPerInf;

					for (j = 0; j < varsPerInf; j++) {
						cv->cvPts[solvingPt].data[inf*varsPerInf + j] += (1.0 / (numNeigh+1)) * 
							cv->cvPts[cv->neighbors[solvingPt][neigh]].backup[neighIndex + j];
					}
				}
			}
		}
	}
}

void CVGoalFunction::normalizeN()  {
	/*
	int comp, ex, i;
	double sum, sum2;

	if (numExamples < 1)
		return;

	int compSize = cv->numPts * NUM_P_VALUES * 3;

	for (comp = 0; comp < numComponents; comp++) {
		// calculate standard deviation of n_comp
		sum = 0;
		sum2 = 0;
		for (ex = 0; ex < numExamples; ex++) {
			sum += curN[ex * numComponents + comp];
			sum2 += sqr(curN[ex * numComponents + comp]);
		}

		if (sum2 < 1e-6)
			continue;

		double avg = sum / numExamples;
		double var = sqrt((1.0 / (numExamples-1)) * (sum2 - (1.0 / numExamples) * sum * sum));
		cout << "avg = " << avg << "; var = " << var << endl;
		if (var < 1e-6 || !_finite(var))
			continue;

		// center
		for (i=0; i < compSize; i++) {
			curVars[i] += curVars[(comp+1)*compSize + i] * avg;
		}
		// fix N
		for (ex = 0; ex < numExamples; ex++) {
			curN[ex * numComponents + comp] -= avg;
			curN[ex * numComponents + comp] /= var;
		}
		// scale
		for (i=0; i < compSize; i++) {
			curVars[(comp+1)*compSize + i] *= var;
		}
	}*/
}

double CVGoalFunction::calcFreeEnergy(bool verbose) {
	int ex, pt, inf;
	int curChar = -1;
	double distF = 1.0 / (2.0 * sqr(distDeviation));
	double matchF = 1.0 / (2.0 * sqr(matchDeviation));
	double ret;
	double distErr = 0;
	double regErr = 0;

	for (ex = 0; ex < numExamples; ex++) {
		// calculate the current reconstruction
		cv->updateCurComponents(examples[ex].dofs);
		cv->updateLocalPos(curMu[examples[ex].charID].Ref());

		for (pt = 0; pt < examples[ex].numPoints; pt++) {
			double conf = 1.0;
			if (useConf)
				conf = examples[ex].conf[pt];
			distErr += conf * distF * (examples[ex].points[pt] - cv->cvPts[pt].localPos).length2();
		}
	}

	if (useRegularization) {
		// regularization
		int varsPerInf = cv->numComponents*(cv->cvPts[0].dataParts / cv->cvPts[0].numInfluences)*3;
		int neigh, j;
		for (pt = 0; pt < cv->numPts; pt++) {
			for (inf = 1; inf < cv->cvPts[pt].numInfluences; inf++) {
				int index = inf*varsPerInf;

				int numNeigh = (int)cv->neighbors[pt].size();
				for (neigh = 0; neigh < numNeigh; neigh++) {
					double neighVal = 0;
					int neighIndex = -1;
					int neighInf = cv->cvPts[cv->neighbors[pt][neigh]].getInfIndex(cv->cvPts[pt].influences[inf]);
					if (neighInf >= 0) {
						neighIndex = neighInf * varsPerInf;

						for (j = 0; j < varsPerInf; j++) {
							double dist = cv->cvPts[pt].data[inf*varsPerInf + j] - cv->cvPts[cv->neighbors[pt][neigh]].backup[neighIndex + j];
							regErr += matchF * sqr(dist);
						}
					}
					else {
						for (j = 0; j < varsPerInf; j++) {
							double dist = cv->cvPts[pt].data[inf*varsPerInf + j];
							regErr += matchF * sqr(dist);
						}
					}
				}
			}
		}
	}

	ret = distErr + regErr;

	if (verbose) {
		cout << "free energy: " << ret << " = " << distErr << " + " << regErr << endl;
	}

	return ret;
}


// CVSkinningGF methods =============================================

void CVSkinningGF::init(CVGoalFunction *gf, Skin *sk, Skeleton *mp) {
	int i, index;

	cvGF = gf;
	skin = sk;
	matchPoses = mp;

	// initialize vars
	numVars = 0;
	for (i=0; i < skin->numPts; i++) {
		numVars += skin->points[i].numVars;
	}
	cout << "number of variables: " << numVars << endl;

	vars.resize(numVars);
	grad.resize(numVars);

	index = 0;
	for (i=0; i < skin->numPts; i++) {
		skin->points[i].copyToVars(vars.n + index);
		index += skin->points[i].numVars;
	}
}

void CVSkinningGF::varsToSkin(Vecd &variables) {
	int index, i;

	index = 0;
	for (i=0; i < skin->numPts; i++) {
		skin->points[i].copyFromVars(vars.n + index);
		index += skin->points[i].numVars;
	}
}

double CVSkinningGF::evaluateFunction(Vecd& variables) {
	double err = 0;
	double curErr;
	int index;
	int ex, pt, inf;

	grad.zeroElements();
	varsToSkin(variables);

	for (ex=0; ex < cvGF->numExamples; ex++) {
		skin->skel->copyVals(&matchPoses[ex]);
		skin->skel->updateCoords();
		skin->updatePoints();

		index = 0;

		for (pt=0; pt < skin->numPts; pt++) {
			err += skin->points[pt].calcGrad(grad.n + index, cvGF->examples[ex].points[pt], skin->skel);
			index += skin->points[pt].numVars;
		}
	}

	lastErr = err;
	return err;
}

void CVSkinningGF::evaluateGradient(Vecd& variables, Vecd& gradient) {
	memcpy(gradient.n, grad.n, sizeof(double)*gradient.size());
}

void CVSkinningGF::solverStep() {
//	cout << "val: " << skin->points[0].tWeights[9] << endl;
	cout << "err: " << lastErr << endl;

	int i, j;
	for (i=0; i < scMesh->numPts(); i++) {
		Vec3d color;
		for (j=0; j < skin->points[i].numTInf; j++) 
			if (skin->points[i].tWeights[j] > 0)
				color += skin->points[i].tWeights[j] *
					skin->skel->transforms.getT(skin->points[i].tTransforms[j])->color;
		scMesh->getPtColor(i) = color;
	}
	redrawVNow();
}