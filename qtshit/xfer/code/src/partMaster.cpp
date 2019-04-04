#include "doppel2.h"
#include "cli.h"
#include "uMaster.h"
#include "uSkin.h"
#include "partMaster.h"
#include "partUI.h"
#include "skeleton.h"
#include "markers.h"
#include "trimesh.h"
#include "trimesh_util.h"
#include <float.h>
#include <GL/glu.h>
#include <assert.h>

extern PartUI *partUI;

// from uMaster
extern USkin skin;
extern UDataSet dataSet;
extern TriMesh *uMesh;

PCAPart *curPart = NULL;
bool showPart = false;
bool showTarget = false;
TriMesh *targetMesh = NULL;
MarkerSet *targetMkr = NULL;
static LBFGSSolver *solver = NULL;

// PCA code from Matt's project (anthropomorphic/pca.cpp)

void subtractMean(TVec& mean, std::vector<TVec>& data)
{
    assert(data.size() > 0);

    // resize output mean vector to match number of data vectors
    mean.SetSize(data[0].Elts());
    mean.MakeZero();
    for (int e=0; e<mean.Elts(); e++) {
        // find mean for each element from each data vector 
        for (int v=0; v<data.size(); v++) {
            assert(data[v].Elts() == mean.Elts());
            mean[e] += data[v][e];
        }
        mean[e] /= TVReal(data.size());

        // subtract mean from each data vector
        for (int v=0; v<data.size(); v++) {
            data[v][e] -= mean[e];
        }
    }
}

double calculateCovariance(TVec& a, TVec& b)
{
    assert(a.Elts() == b.Elts());
    assert(a.Elts() > 1);

    // calculate means
    double mean[2] = {0, 0};
    for (int e=0; e<a.Elts(); e++) {
        mean[0] += a[e];
        mean[1] += b[e];
    }
    mean[0] /= double(a.Elts());
    mean[1] /= double(a.Elts());
    
    // calculate covariance
    double covariance = 0;
    for (int e=0; e<a.Elts(); e++) {
        covariance += (a[e] - mean[0]) * (b[e] - mean[1]);
    }
    covariance /= double(a.Elts() - 1); // bias-corrected covariance

    return covariance;
}

void constructDataMatrix(TMat& dataMat, std::vector<TVec>& data)
{
    // resize data matrix
    dataMat.SetSize(data[0].Elts(), (int)data.size());

    // fill in data values 
    for (int r=0; r<data[0].Elts(); r++) {
        for (int c=0; c<data.size(); c++) {
            dataMat[r][c] = data[c][r];
        }
    }
}

void performPCA(std::vector<TVec>& data,  // input
                TVec& mean, TVec& variance, std::vector<TVec>& pcaVectors)  // output
{
    // find means and center data
    subtractMean(mean, data);

    // compute covariance matrix
    TMat dataMatrix;
    constructDataMatrix(dataMatrix, data);

    // perform SVD (and calculate variance)
    TMat U, V;
    assert(dataMatrix.Rows() >= dataMatrix.Cols()); // limitation of vl...
    SVDFactorization(dataMatrix, U, V, variance);

    // copy the resulting data to output format
    cout << "output..." << endl;
    pcaVectors.resize(U.Cols());
    for (int c=0; c<U.Cols(); c++) {
        pcaVectors[c].SetSize(U.Rows());
        for (int r=0; r<U.Rows(); r++) {
            pcaVectors[c][r] = U[r][c];
        }
        pcaVectors[c].Normalise();
    }
}



bool PCAPart::save(char *fname) {
	FILE *f;
	if (!openFile(&f, fname, "wb", "PCA part"))
		return false;

	fwrite(&numPts, sizeof(int), 1, f);
	fwrite(&numComponents, sizeof(int), 1, f);
	fwrite(weights, sizeof(double), numPts, f);
	fwrite(mean.Ref(), sizeof(double), numPts * 3, f);
	int comp;
	for (comp = 0; comp < numComponents; comp++) {
		fwrite(components[comp].Ref(), sizeof(double), numPts * 3, f);
	}
	fwrite(pcaW.n, sizeof(double), numComponents, f);
	fwrite(trans.n, sizeof(double), 16, f);

	fclose(f);

	char fname2[1024];
	strncpy(fname2, fname, 1023);
	fname2[strlen(fname2)-3] = 'p';
	fname2[strlen(fname2)-2] = 'l';
	fname2[strlen(fname2)-1] = 'y';
	mesh->savePly(fname2);
	return true;
}

bool PCAPart::load(char *fname) {
	FILE *f;
	if (!openFile(&f, fname, "rb", "PCA part"))
		return false;

	fread(&numPts, sizeof(int), 1, f);
	fread(&numComponents, sizeof(int), 1, f);
	pcaW.resize(numComponents, true);
	weights = new double[numPts];
	fread(weights, sizeof(double), numPts, f);
	mean.SetSize(numPts * 3);
	fread(mean.Ref(), sizeof(double), numPts * 3, f);
	int comp;
	components.resize(numComponents);
	for (comp = 0; comp < numComponents; comp++) {
		components[comp].SetSize(numPts * 3);
		fread(components[comp].Ref(), sizeof(double), numPts * 3, f);
	}
	fread(pcaW.n, sizeof(double), numComponents, f);
	fread(trans.n, sizeof(double), 16, f);

	fclose(f);

	char fname2[1024];
	strncpy(fname2, fname, 1023);
	fname2[strlen(fname2)-3] = 'p';
	fname2[strlen(fname2)-2] = 'l';
	fname2[strlen(fname2)-1] = 'y';
	mesh = new TriMesh();
	mesh->loadPly(fname2);
	return true;

}

void PCAPart::updateMesh() {
	int pt, comp;
	for (pt = 0; pt < mesh->numPts(); pt++) {
		Vec3d v(mean[pt*3+0], mean[pt*3+1], mean[pt*3+2]);

		for (comp = 0; comp < numComponents; comp++) {
			if (pcaW[comp] != 0)
				v += pcaW[comp] * 
					Vec3d(components[comp][pt*3+0], 
						components[comp][pt*3+1], 
						components[comp][pt*3+2]);
		}

		mesh->getPt(pt) = v;
	}

	mesh->calcNormals();
}



void PartUI::drawGL() {
	if (showPart && curPart) {
		renderTriMesh(curPart->mesh, dispMode);
	}
	if (showTarget && targetMesh) {
		renderTriMesh(targetMesh, dispMode);
	}
}

void partBuildPCA(const char *params) {
	char fname[80], transName[80];
	strcpy(fname, "-");
	params = extractString(params, fname, 80);

	PCAPart *newPart = new PCAPart();
	int pt, inf, ch, tri;
	bool hasPts = false;
	SkelTransform *primaryTransform = NULL;
	int *vertMap = new int[skin.numPts];
	newPart->weights = new double[skin.numPts];
	memset(newPart->weights, 0, sizeof(double) * skin.numPts);

	while (1) {
		transName[0] = 0;
		params = extractString(params, transName, 80);
		if (transName[0] == 0)
			break;

		int transID = skin.skel->transforms.lookupName(transName);
		if (transID < 0) {
			cout << "unknown transform: " << transName << endl;
			continue;
		}

		if (primaryTransform == NULL)
			primaryTransform = skin.skel->transforms.getT(transID);

		for (pt = 0; pt < skin.numPts; pt++) {
			for (inf = 0; inf < skin.maxInf; inf++) {
				if (skin.infJoints[pt * skin.maxInf + inf] == transID) {
					double w = skin.infWeights[pt * skin.maxInf + inf];
					if (w > 0) {
						hasPts = true;
						newPart->weights[pt] += w;
					}
				}
			}
		}
	}

	if (hasPts) {
		// make the mesh
		newPart->numPts = 0;
		newPart->mesh = new TriMesh();
		newPart->mesh->gsPtColors(1);
		for (pt = 0; pt < skin.numPts; pt++) {
			if (newPart->weights[pt] > 0) {
				newPart->mesh->addPoint(Vec3d());
				newPart->mesh->getPtColor(newPart->numPts) = Vec3d(0, 1, 0);
				vertMap[pt] = newPart->numPts;

				newPart->numPts++;
			}
			else {
				vertMap[pt] = -1;
			}
		}
		for (tri = 0; tri < uMesh->numTris(); tri++) {
			int v0 = vertMap[uMesh->getTri(tri, 0)];
			int v1 = vertMap[uMesh->getTri(tri, 1)];
			int v2 = vertMap[uMesh->getTri(tri, 2)];
			if (v0 >= 0 && v1 >= 0 && v2 >= 0)
				newPart->mesh->addTri(v0, v1, v2);
		}
		cout << "part has " << newPart->numPts << " points, and " << 
			newPart->mesh->numTris() << " triangles" << endl;

		vector<TVec> data;

		for (ch = 0; ch < dataSet.numCharacters; ch++) {
			Mat4d localTrans;

			data.push_back(TVec(newPart->numPts * 3));

			char s[1024];
			ifstream in;
			sprintf(s, "data/skels/%sa.po.txt", dataSet.examples[dataSet.charIndex[ch]].charName);
			if (openIFStream(&in, s, "pose file")) {
				skin.skel->loadPose(in);
				in.close();
			}
			localTrans = primaryTransform->globalCoord.mat;
			localTrans = localTrans.inverse();

			int curInd = 0;
			for (pt = 0; pt < skin.numPts; pt++) {
				if (newPart->weights[pt] > 0) {
					Vec3d v = dataSet.examples[dataSet.charIndex[ch]].pts[pt];
					v = localTrans * v;
					data[ch][curInd++] = v[0];
					data[ch][curInd++] = v[1];
					data[ch][curInd++] = v[2];
				}
			}
		}

		cout << "running PCA" << endl;
		baTimer();
		performPCA(data, newPart->mean, newPart->variance, newPart->components);
		cout << "took " << baTimer() << " ms" << endl;
		cout << newPart->variance << endl;

		newPart->numComponents = (int)newPart->components.size();
		newPart->pcaW.resize(newPart->numComponents, true);
	}
	else {
		cout << "no points found" << endl;
	}

	cout << "done" << endl;
	curPart = newPart;
	updatePartWeights(NULL, 0);
	if (strlen(fname) > 1)
		newPart->save(fname);
}

void partLoadPCA(const char *params) {
	char fname[80];
	params = extractString(params, fname, 80);
	curPart = new PCAPart();
	curPart->load(fname);
	updatePartWeights(NULL, 0);
}

void loadTargetMesh(const char *params) {
	char fname[256];
	params = extractString(params, fname, 256);
	char fname2[256];
	fname2[0] = 0;
	params = extractString(params, fname2, 256);

	targetMesh = new TriMesh();
	targetMesh->loadFile(fname);
	targetMesh->calcNormals();

	if (strlen(fname2) > 1) {
		targetMkr = new MarkerSet();
		targetMkr->load(fname2);
	}

	redrawV();
}

void updatePartWeights(double *w, int numW) {
	int comp;

	if (!curPart)
		return;

	curPart->pcaW.zeroElements();
	for (comp = 0; comp < numW; comp++)
		curPart->pcaW[comp] = w[comp];

	curPart->updateMesh();
	redrawV();
}

void partSolve(const char *params) {
	if (!curPart) {
		cout << "no active part" << endl;
		return;
	}
	if (solver) {
		cout << "solver already running" << endl;
		return;
	}

	char mode[80];
	int iter = 1000000;
	params = extractString(params, mode, 80);
	params = extractInt(params, &iter);

//	PartMatchFunction partMatchFn(curPart);

//	solver = new LBFGSSolver(partMatchFn);
}

void initPartMaster() {
	registerFunction(partBuildPCA, "partBuildPCA");
	registerFunction(partLoadPCA, "partLoadPCA");
	registerFunction(loadTargetMesh, "partLoadTarget");
	registerFunction(partSolve, "partSolve");
}