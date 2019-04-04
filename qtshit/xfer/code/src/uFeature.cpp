#include "doppel2.h"
#include "uFeature.h"
#include "ba.h"
#include "trimesh.h"
#include "uMaster.h"
#include "vl/VLd.h"

class FeatureMat {
public:
	int numPCA, numFeatures;
	float *data;
	float *invData;
	char80 *name, *iName;
	float *minV, *maxV, *mult, *iMult, *power;

	void load(char *fname, char *invFName);
	float v(int i, int f);
	void init(int iNumFeatures);
};


int maxFeatures = 10;
FeatureMat fCirc;
float *featureBase;
int numPCA;
float *curPCA, *curBase;



// FeatureMat class =================================================
// This class encapsulates the data for editing a particular set of
// features (sliders).

// load a feature matrix (and its inverse) from text files
void FeatureMat::load(char *fname, char *invFName) {
	int i, j;
	ifstream in;
	if (!openIFStream(&in, fname, "feature matrix"))
		return;

	in >> numPCA >> numFeatures;
	numFeatures--;
	data = new float[numPCA*(numFeatures+1)];
	for (i=0; i < numPCA; i++) {
		for (j=0; j < numFeatures+1; j++) {
			in >> data[i * (numFeatures+1) + j];
		}
	}
	in.close();

	if (!openIFStream(&in, invFName, "inverse feature matrix"))
		return;
	invData = new float[numPCA*numFeatures];
	for (i=0; i < numFeatures; i++) {
		for (j=0; j < numPCA; j++) {
			in >> invData[i * numPCA + j];
		}
	}
	in.close();

	name = new char80[numFeatures];
	iName = new char80[numFeatures];
	minV = new float[numFeatures];
	maxV = new float[numFeatures];
	mult = new float[numFeatures];
	iMult = new float[numFeatures];
	power = new float[numFeatures];
	for (i=0; i < numFeatures; i++) power[i] = 1;
}

void FeatureMat::init(int n) {
	int i;
	numFeatures = n;
	name = new char80[numFeatures];
	iName = new char80[numFeatures];
	minV = new float[numFeatures];
	maxV = new float[numFeatures];
	mult = new float[numFeatures];
	iMult = new float[numFeatures];
	power = new float[numFeatures];
	for (i=0; i < numFeatures; i++) power[i] = 1;
}

// look up the i'th entry for feature f in the feature matrix
float FeatureMat::v(int i, int f) {
	if (f < 0)
		return data[i * (numFeatures+1) + numFeatures];

	if (f >= numFeatures)
		return 0;

	return data[i * (numFeatures+1) + f];
}


// ==================================================================


double convertf(int f, double d) {
	return pow((1.0 / fCirc.iMult[f]) * d, 1.0/fCirc.power[f]);
}

double unconvertf(int f, double d) {
	return fCirc.iMult[f] * pow(d, (double)fCirc.power[f]);
}

// given a set of feature (slider) values, update the mesh relative to curBase
void setFeatures(float *fVals, double *pca) {
	int i, j;

	// update curPCA
	for (i=0; i < numPCA; i++) {
		pca[i+1] = curBase[i];

		for (j=0; j < fCirc.numFeatures; j++) {
			double val = convertf(j, fVals[j]) - featureBase[j];
			pca[i+1] += fCirc.v(i, j) * val;
		}
	}
}

void setFeatureBase(float *fVals, double *pca) {
	int i;
	for (i=0; i < numPCA; i++)
		curBase[i] = pca[i+1];
	for (i=0; i < fCirc.numFeatures; i++)
		featureBase[i] = convertf(i, fVals[i]);
}

/*
// update the current feature (slider) values based on curBase
void updateFromBase() {
	int i, j;
	memset(featureBase, 0, sizeof(float)*maxFeatures);

	int numF = fCirc.numFeatures;
	for (i=0; i < numF; i++) {
		for (j=0; j < numPCA; j++) {
			featureBase[i] += fCirc.invData[i*fCirc.numPCA + j] * 
				(curBase[j] - fCirc.data[j*(numF+1) + numF]);
		}
	}
	featureBase[i] = 1;
}

// make the current PCA values the base PCA values
void resetFromCur() {
	memcpy(curBase, curPCA, sizeof(float)*numPCA);
	updateFromBase();
}*/


void loadFeatures(char *fname) {
	ifstream fSet, in;
	char s[256], s2[256];
	int i, j, k, id, ind, numFeatures;
	double d;
	FILE *f;

	if (!openIFStream(&fSet, fname, "feature set")) {
		return;
	}

	fSet >> numFeatures;
	fCirc.init(numFeatures);

	vector<double> *featureVals = new vector<double>[dataSet.numCharacters];

	for (i=0; i < fCirc.numFeatures; i++) {
		fSet >> s;
		if (!openIFStream(&in, s, "feature file")) {
			return;
		}
		cout << "loading " << s << endl;

		in >> s; // name
		in >> s; // unit
		in >> d >> d >> d;

		fSet >> fCirc.mult[i] >> fCirc.iMult[i] >> 
				fCirc.power[i] >> fCirc.minV[i] >> fCirc.maxV[i];

		while (in.good()) {
			in >> s2 >> d;
			id = atoi(s2);

			for (j=0; j < dataSet.numCharacters; j++) {
				char *charName = dataSet.examples[dataSet.charIndex[j]].charName;
				if (charName[0] >= '0' && charName[0] <= '9' && atoi(charName) == id) {
					featureVals[j].push_back(pow(d, 1.0/fCirc.power[i]));
					break;
				}
			}
		}
		in.close();
	}

	int numComplete = 0;
	for (i = 0; i < dataSet.numCharacters; i++) {
		if (featureVals[i].size() == numFeatures)
			numComplete++;
	}

	cout << numComplete << " examples being used" << endl;
	if (numComplete == 0) return;

	VLMatd m(numComplete, numFeatures+1);
	ind = 0;
	int *pcaEx = new int[numComplete];
	for (i = 0; i < dataSet.numCharacters; i++) {
		if (featureVals[i].size() == numFeatures) {
			for (j = 0; j < numFeatures; j++)
				m[ind][j] = featureVals[i][j];
			m[ind][j] = 1;
			pcaEx[ind] = i;

			ind++;
		}
	}

	memset(featureBase, 0, sizeof(float)*maxFeatures);
	memset(curBase, 0, sizeof(float)*numPCA);
	for (i=0; i < numFeatures; i++) {
		for (j=0; j < numComplete; j++)
			featureBase[i] += m[j][i];
		featureBase[i] /= numComplete;
		cout << unconvertf(i, featureBase[i]) << " ";
	}
	cout << endl;

	// calculate pseudoinverse
	TMat U(numComplete, numFeatures+1), V(numFeatures+1, numFeatures+1);
	TVec diagonal(numFeatures+1);
	SVDFactorization(m, U, V, diagonal);
	TMat pinv(numFeatures+1, numComplete), diag(numFeatures+1, numFeatures+1);
	diag = vl_0;
	for (j = 0; j < numFeatures+1; j++) {
		if (fabs(diagonal[j]) > 0.00001)
			diag.Elt(j, j) = 1.0 / diagonal[j];
	}
	pinv = V * diag * trans(U);

	fCirc.data = new float[numPCA*(numFeatures+1)];
	for (i=0; i < numPCA; i++) {
		for (j=0; j < numFeatures+1; j++) {
			double sum = 0;
			for (k=0; k < numComplete; k++) {
				sum += dataSet.charMu[pcaEx[k]][i+1] * pinv[j][k];
			}
			fCirc.data[i * (numFeatures+1) + j] = sum;
		}
	}
}


void initFeatures(int iNumPCA) {
	numPCA = iNumPCA;

	featureBase = new float[maxFeatures];
	curPCA = new float[numPCA];
	curBase = new float[numPCA];
}
