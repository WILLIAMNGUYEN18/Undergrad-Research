#include "doppel2.h"
#include "trimesh.h"
#include "ply.h"

struct PlyVertex {
	int face;
	float barya, baryb, baryc;
};

struct PlyFace {
  unsigned char nverts;    /* number of vertex indices in list */
  int *verts;              /* vertex index list */
};

static PlyProperty vert_props[] = { // list of property information for a vertex
  {"face", PLY_INT, PLY_INT, offsetof(PlyVertex,face), 0, 0, 0, 0},
  {"barya", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,barya), 0, 0, 0, 0},
  {"baryb", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,baryb), 0, 0, 0, 0},
  {"baryc", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,baryc), 0, 0, 0, 0},
};

static PlyProperty face_props[] = { /* list of property information for a face */
  {"vertex_indices", PLY_INT, PLY_INT, offsetof(PlyFace,verts),
   1, PLY_UCHAR, PLY_UCHAR, offsetof(PlyFace,nverts)},
};

struct PlyTristrip {
  int nverts;    /* number of vertex indices in list */
  int *verts;              /* vertex index list */
};

static PlyProperty tristrips_props[] = { /* list of property information for a tristrip */
  {"vertex_indices", PLY_INT, PLY_INT, offsetof(PlyTristrip,verts),
   1, PLY_INT, PLY_INT, offsetof(PlyTristrip,nverts)},
};


void createMapping(char *fnameS, char *fnameB, char *fnameM) {
	vector<int> smallTris;
	vector<PlyVertex> bigVerts;
	int numSmallVerts, numBigVerts;

	// load triangles in small mesh
	FILE *f;
	if (!openFile(&f, fnameS, "rb", "small ply"))
		return;

	int nelems;
	char **elemNames;
	PlyFile *ply = ply_read(f, &nelems, &elemNames);
	if (ply == NULL) {
		cerr << "couldn't read ply file " << fnameS << endl;
		return;
	}

	// go through each kind of element that we learned is in the file 
	// and read them
	char *elemName;
	PlyProperty **plist;
	int numElems, nprops;
	int i, j;
	for (i = 0; i < nelems; i++) {
		// get the description of the first element
		elemName = elemNames[i];
		plist = ply_get_element_description(ply, elemName, &numElems, &nprops);

		// if we're on vertex elements, read them in
		if (strcmp("vertex", elemName) == 0) {
			numSmallVerts = numElems;

			// set up for getting vertex elements
			ply_get_property(ply, elemName, &vert_props[0]);
			ply_get_property(ply, elemName, &vert_props[1]);
			ply_get_property(ply, elemName, &vert_props[2]);
			ply_get_property(ply, elemName, &vert_props[4]);

			// grab all the vertex elements
			for (j = 0; j < numElems; j++) {
				// grab an element from the file
				PlyVertex vert;
				ply_get_element(ply, &vert);
			}
		}

		// if we're on face elements, read them in
		if (strcmp("face", elemName) == 0) {
			// set up for getting face elements
			ply_get_property (ply, elemName, &face_props[0]);
      
			// grab all the face elements
			for (j = 0; j < numElems; j++) {
				// grab an element from the file
				PlyFace face;
				ply_get_element(ply, &face);

				smallTris.push_back(face.verts[0]);
				smallTris.push_back(face.verts[1]);
				smallTris.push_back(face.verts[2]);

				free(face.verts);
			}
		}

		// if we're on tristrips elements, read them in
		if (strcmp("tristrips", elemName) == 0) {
			// set up for getting tristrips
			ply_get_property (ply, elemName, &tristrips_props[0]);
      
			// grab all the face elements
			for (j = 0; j < numElems; j++) {
				// grab an element from the file
				PlyTristrip strip;
				ply_get_element(ply, &strip);

				bool clockwise = true;
				for (i=2; i < strip.nverts; i++) {
					if (strip.verts[i] == -1) {
						i += 2;
						clockwise = true;
					}
					else {
						if (clockwise) {
							smallTris.push_back(strip.verts[i-2]);
							smallTris.push_back(strip.verts[i-1]);
							smallTris.push_back(strip.verts[i-0]);
						}
						else {
							smallTris.push_back(strip.verts[i-2]);
							smallTris.push_back(strip.verts[i-0]);
							smallTris.push_back(strip.verts[i-1]);
						}
						clockwise = !clockwise;
					}
				}

				free(strip.verts);
			}
		}
	}

	ply_close(ply);


	// load vertices in big mesh
	if (!openFile(&f, fnameB, "rb", "big ply"))
		return;

	ply = ply_read(f, &nelems, &elemNames);
	if (ply == NULL) {
		cerr << "couldn't read ply file " << fnameB << endl;
		return;
	}

	// go through each kind of element that we learned is in the file 
	// and read them
	for (i = 0; i < nelems; i++) {
		// get the description of the first element
		elemName = elemNames[i];
		plist = ply_get_element_description(ply, elemName, &numElems, &nprops);

		// if we're on vertex elements, read them in
		if (strcmp("vertex", elemName) == 0) {
			numBigVerts = numElems;

			// set up for getting vertex elements
			ply_get_property(ply, elemName, &vert_props[0]);
			ply_get_property(ply, elemName, &vert_props[1]);
			ply_get_property(ply, elemName, &vert_props[2]);
			ply_get_property(ply, elemName, &vert_props[4]);

			// grab all the vertex elements
			for (j = 0; j < numElems; j++) {
				// grab an element from the file
				PlyVertex vert;
				ply_get_element(ply, &vert);
				vert.baryc = 1.0 - vert.barya - vert.baryb;
				bigVerts.push_back(vert);
			}
		}

		// if we're on face elements, read them in
		if (strcmp("face", elemName) == 0) {
			// set up for getting face elements
			ply_get_property (ply, elemName, &face_props[0]);
      
			// grab all the face elements
			for (j = 0; j < numElems; j++) {
				// grab an element from the file
				PlyFace face;
				ply_get_element(ply, &face);
				free(face.verts);
			}
		}

		// if we're on tristrips elements, read them in
		if (strcmp("tristrips", elemName) == 0) {
			// set up for getting tristrips
			ply_get_property (ply, elemName, &tristrips_props[0]);
      
			// grab all the face elements
			for (j = 0; j < numElems; j++) {
				// grab an element from the file
				PlyTristrip strip;
				ply_get_element(ply, &strip);
				free(strip.verts);
			}
		}
	}

	ply_close(ply);

	cout << "building mapping between " << numSmallVerts << " and " << numBigVerts << " vertices..." << endl;

	if (!openFile(&f, fnameM, "wb", "mapping file"))
		return;

	fwrite(&numSmallVerts, sizeof(int), 1, f);
	fwrite(&numBigVerts, sizeof(int), 1, f);

	for (i=0; i < numBigVerts; i++) {
		j = 3;
		fwrite(&j, sizeof(int), 1, f);

		j = smallTris[bigVerts[i].face * 3 + 0];
		fwrite(&j, sizeof(int), 1, f);
		fwrite(&bigVerts[i].barya, sizeof(float), 1, f);

		j = smallTris[bigVerts[i].face * 3 + 1];
		fwrite(&j, sizeof(int), 1, f);
		fwrite(&bigVerts[i].baryb, sizeof(float), 1, f);

		j = smallTris[bigVerts[i].face * 3 + 2];
		fwrite(&j, sizeof(int), 1, f);
		fwrite(&bigVerts[i].baryc, sizeof(float), 1, f);
	}

	fclose(f);
}

/* mirror-map version

void createMapping(char *fnameS, char *fnameB, char *fnameM) {
	vector<int> smallTris;
	vector<PlyVertex> bigVerts;
	int numSmallVerts, numBigVerts;

	// stupid hack: have to load mirror mappings
	FILE *mf;
	if (!openFile(&mf, "data/mmap-james-10k.dat", "rb", "small mirror map"))
		return;
	int mSmallN, *mSmallMap;
	fread(&mSmallN, sizeof(int), 1, mf);
	mSmallMap = new int[mSmallN];
	fread(mSmallMap, sizeof(int), mSmallN, mf);
	fclose(mf);
	if (!openFile(&mf, "data/mmap-james-50k.dat", "rb", "big mirror map"))
		return;
	int mBigN, *mBigMap;
	fread(&mBigN, sizeof(int), 1, mf);
	mBigMap = new int[mBigN];
	fread(mBigMap, sizeof(int), mBigN, mf);
	fclose(mf);

	// load triangles in small mesh
	FILE *f;
	if (!openFile(&f, fnameS, "rb", "small ply"))
		return;

	int nelems;
	char **elemNames;
	PlyFile *ply = ply_read(f, &nelems, &elemNames);
	if (ply == NULL) {
		cerr << "couldn't read ply file " << fnameS << endl;
		return;
	}

	// go through each kind of element that we learned is in the file 
	// and read them
	char *elemName;
	PlyProperty **plist;
	int numElems, nprops;
	int i, j;
	for (i = 0; i < nelems; i++) {
		// get the description of the first element
		elemName = elemNames[i];
		plist = ply_get_element_description(ply, elemName, &numElems, &nprops);

		// if we're on vertex elements, read them in
		if (strcmp("vertex", elemName) == 0) {
			numSmallVerts = numElems;

			// set up for getting vertex elements
			ply_get_property(ply, elemName, &vert_props[0]);
			ply_get_property(ply, elemName, &vert_props[1]);
			ply_get_property(ply, elemName, &vert_props[2]);
			ply_get_property(ply, elemName, &vert_props[4]);

			// grab all the vertex elements
			for (j = 0; j < numElems; j++) {
				// grab an element from the file
				PlyVertex vert;
				ply_get_element(ply, &vert);
			}
		}

		// if we're on face elements, read them in
		if (strcmp("face", elemName) == 0) {
			// set up for getting face elements
			ply_get_property (ply, elemName, &face_props[0]);
      
			// grab all the face elements
			for (j = 0; j < numElems; j++) {
				// grab an element from the file
				PlyFace face;
				ply_get_element(ply, &face);

				smallTris.push_back(face.verts[0]);
				smallTris.push_back(face.verts[1]);
				smallTris.push_back(face.verts[2]);

				free(face.verts);
			}
		}

		// if we're on tristrips elements, read them in
		if (strcmp("tristrips", elemName) == 0) {
			// set up for getting tristrips
			ply_get_property (ply, elemName, &tristrips_props[0]);
      
			// grab all the face elements
			for (j = 0; j < numElems; j++) {
				// grab an element from the file
				PlyTristrip strip;
				ply_get_element(ply, &strip);

				bool clockwise = true;
				for (i=2; i < strip.nverts; i++) {
					if (strip.verts[i] == -1) {
						i += 2;
						clockwise = true;
					}
					else {
						if (clockwise) {
							smallTris.push_back(strip.verts[i-2]);
							smallTris.push_back(strip.verts[i-1]);
							smallTris.push_back(strip.verts[i-0]);
						}
						else {
							smallTris.push_back(strip.verts[i-2]);
							smallTris.push_back(strip.verts[i-0]);
							smallTris.push_back(strip.verts[i-1]);
						}
						clockwise = !clockwise;
					}
				}

				free(strip.verts);
			}
		}
	}

	ply_close(ply);


	// load vertices in big mesh
	if (!openFile(&f, fnameB, "rb", "big ply"))
		return;

	ply = ply_read(f, &nelems, &elemNames);
	if (ply == NULL) {
		cerr << "couldn't read ply file " << fnameB << endl;
		return;
	}

	// go through each kind of element that we learned is in the file 
	// and read them
	for (i = 0; i < nelems; i++) {
		// get the description of the first element
		elemName = elemNames[i];
		plist = ply_get_element_description(ply, elemName, &numElems, &nprops);

		// if we're on vertex elements, read them in
		if (strcmp("vertex", elemName) == 0) {
			numBigVerts = numElems;

			// set up for getting vertex elements
			ply_get_property(ply, elemName, &vert_props[0]);
			ply_get_property(ply, elemName, &vert_props[1]);
			ply_get_property(ply, elemName, &vert_props[2]);
			ply_get_property(ply, elemName, &vert_props[4]);

			// grab all the vertex elements
			for (j = 0; j < numElems; j++) {
				// grab an element from the file
				PlyVertex vert;
				ply_get_element(ply, &vert);
				vert.baryc = 1.0 - vert.barya - vert.baryb;
				bigVerts.push_back(vert);
			}
		}

		// if we're on face elements, read them in
		if (strcmp("face", elemName) == 0) {
			// set up for getting face elements
			ply_get_property (ply, elemName, &face_props[0]);
      
			// grab all the face elements
			for (j = 0; j < numElems; j++) {
				// grab an element from the file
				PlyFace face;
				ply_get_element(ply, &face);
				free(face.verts);
			}
		}

		// if we're on tristrips elements, read them in
		if (strcmp("tristrips", elemName) == 0) {
			// set up for getting tristrips
			ply_get_property (ply, elemName, &tristrips_props[0]);
      
			// grab all the face elements
			for (j = 0; j < numElems; j++) {
				// grab an element from the file
				PlyTristrip strip;
				ply_get_element(ply, &strip);
				free(strip.verts);
			}
		}
	}

	ply_close(ply);

//	cout << "building mapping between " << numSmallVerts << " and " << numBigVerts << " vertices..." << endl;
	cout << "building mapping between " << mSmallN << " and " << mBigN << " vertices..." << endl;

	if (!openFile(&f, fnameM, "wb", "mapping file"))
		return;

	fwrite(&mSmallN, sizeof(int), 1, f);	// numSmallVerts
	fwrite(&mBigN, sizeof(int), 1, f);		// numBigVerts

	for (i=0; i < numBigVerts; i++) {
		j = 3;
		fwrite(&j, sizeof(int), 1, f);

		j = smallTris[bigVerts[i].face * 3 + 0];
		fwrite(&j, sizeof(int), 1, f);
		fwrite(&bigVerts[i].barya, sizeof(float), 1, f);

		j = smallTris[bigVerts[i].face * 3 + 1];
		fwrite(&j, sizeof(int), 1, f);
		fwrite(&bigVerts[i].baryb, sizeof(float), 1, f);

		j = smallTris[bigVerts[i].face * 3 + 2];
		fwrite(&j, sizeof(int), 1, f);
		fwrite(&bigVerts[i].baryc, sizeof(float), 1, f);
	}
	for (; i < mBigN; i++) {
		j = 3;
		fwrite(&j, sizeof(int), 1, f);

		j = mSmallMap[smallTris[bigVerts[mBigMap[i]].face * 3 + 0]];
		fwrite(&j, sizeof(int), 1, f);
		fwrite(&bigVerts[mBigMap[i]].barya, sizeof(float), 1, f);

		j = mSmallMap[smallTris[bigVerts[mBigMap[i]].face * 3 + 1]];
		fwrite(&j, sizeof(int), 1, f);
		fwrite(&bigVerts[mBigMap[i]].baryb, sizeof(float), 1, f);

		j = mSmallMap[smallTris[bigVerts[mBigMap[i]].face * 3 + 2]];
		fwrite(&j, sizeof(int), 1, f);
		fwrite(&bigVerts[mBigMap[i]].baryc, sizeof(float), 1, f);
	}

	fclose(f);
}
*/