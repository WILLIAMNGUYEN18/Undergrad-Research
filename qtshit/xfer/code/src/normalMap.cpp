#include "normalMap.h"
#include "surfaceDef.h"
#include "ba.h"

const bool gInTangentSpace = true;

void calcTexBary(TriMesh *srcMesh, TexBary *texBary) {
	int tri, i, x, y;

	// initialize to zero
	memset(texBary, 0, NORMAL_MAP_SIZE * sizeof(TexBary));

	for (tri = 0; tri < srcMesh->numTris(); tri++) {
		// get this triangle's texture coordinates
		double texX[3], texY[3];
		for (i=0; i < 3; i++) {
			texX[i] = srcMesh->getTriTexCoords(tri, i)[0] * NORMAL_MAP_W;
			texY[i] = srcMesh->getTriTexCoords(tri, i)[1] * NORMAL_MAP_H;
		}

		double bary0 = ((texX[1] - texX[0]) * (texY[2] - texY[0]) - 
			(texX[2] - texX[0]) * (texY[1] - texY[0]));

		// make sure bary0 is not zero
		if (fabs(bary0) < 1e-10)
			continue;

		// find minimum and maximum Y
		int minY = ceil(min(texY[0], min(texY[1], texY[2])));
		int maxY = floor(max(texY[0], max(texY[1], texY[2])));

		// pad outside the triangle by a couple of pixels
		int pad = 2;

		// loop over scanlines
		for (y = minY - pad; y <= maxY + pad; y++) {
			// make sure we're within the texture map
			if ((y < 0) || (y >= NORMAL_MAP_H)) {
				continue;
			}

			// now find the minimum and maximum X for this scanline
			// if we're above or below the triangle, use the bounds from the top or bottom
			int boundsY = y;
			if (boundsY < minY)
				boundsY = minY;
			if (boundsY > maxY)
				boundsY = maxY;
			int minX = NORMAL_MAP_W-1;
			int maxX = 0;
			double w, tempX;
			// figure out which triangle edges are on this scanline, and grow the bounds
			// (this could be done more efficiently)
			for (i=0; i < 3; i++) {
				int i2 = (i+1) % 3;
				if (texY[i] <= boundsY && texY[i2] >= boundsY) {
					w = (boundsY - texY[i]) / (texY[i2] - texY[i]);
					tempX = (1.0-w) * texX[i] + w * texX[i2];
					if (ceil(tempX) < minX)
						minX = ceil(tempX);
					if (floor(tempX) > maxX)
						maxX = floor(tempX);
				}
				else if (texY[i2] <= boundsY && texY[i] >= boundsY) {
					w = (boundsY - texY[i2]) / (texY[i] - texY[i2]);
					tempX = (1.0-w) * texX[i2] + w * texX[i];
					if (ceil(tempX) < minX)
						minX = ceil(tempX);
					if (floor(tempX) > maxX)
						maxX = floor(tempX);
				}
			}

			for (x = minX - pad; x <= maxX + pad; x++) {
				// make sure we're within the texture map
				if ((x < 0) || (x >= NORMAL_MAP_W)) {
					continue;
				}

				// find the barycentric coordinates
				double b1 = ((texX[1]-x) * (texY[2]-y) - (texX[2]-x) * (texY[1]-y)) / bary0;
				double b2 = ((texX[2]-x) * (texY[0]-y) - (texX[0]-x) * (texY[2]-y)) / bary0;
				double b3 = ((texX[0]-x) * (texY[1]-y) - (texX[1]-x) * (texY[0]-y)) / bary0;

				// if we're outside the triangle, then only update if there's no 'real'
				// record stored already
				if ((b1 > 1.0) || (b1 < 0.0) ||
					(b2 > 1.0) || (b2 < 0.0) ||
					(b3 > 1.0) || (b3 < 0.0)) {
					if (texBary[y*NORMAL_MAP_W + x].exists)
						continue;
				}

				// store the info for this sample
				texBary[y*NORMAL_MAP_W + x].exists = true;
				texBary[y*NORMAL_MAP_W + x].tri = tri;
				texBary[y*NORMAL_MAP_W + x].bary[0] = b1;
				texBary[y*NORMAL_MAP_W + x].bary[1] = b2;
				texBary[y*NORMAL_MAP_W + x].bary[2] = b3;
			}
 		}
 	}
}

void calcTangentSpace(TriMesh *mesh, Mat3d **mats) {
	if (*mats == NULL)
		*mats = new Mat3d[mesh->numTris() * 3];

	int tri, pt;
	for (tri = 0; tri < mesh->numTris(); tri++) {
		for (pt = 0; pt < 3; pt++) {
			int pt2 = (pt+1) % 3;
			int pt3 = (pt+2) % 3;
			Mat3d m;
			double det = 
				(mesh->getTriTexCoords(tri, pt2)[0] - mesh->getTriTexCoords(tri, pt)[0]) *
				(mesh->getTriTexCoords(tri, pt3)[1] - mesh->getTriTexCoords(tri, pt)[1]) -
				(mesh->getTriTexCoords(tri, pt3)[0] - mesh->getTriTexCoords(tri, pt)[0]) *
				(mesh->getTriTexCoords(tri, pt2)[1] - mesh->getTriTexCoords(tri, pt)[1]);
			Vec3d tangent = (1.0 / det) * (
				(mesh->getTriTexCoords(tri, pt3)[1] - mesh->getTriTexCoords(tri, pt)[1]) *
				(mesh->getPt(mesh->getTri(tri, pt2)) - mesh->getPt(mesh->getTri(tri, pt))) -
				(mesh->getTriTexCoords(tri, pt2)[1] - mesh->getTriTexCoords(tri, pt)[1]) *
				(mesh->getPt(mesh->getTri(tri, pt3)) - mesh->getPt(mesh->getTri(tri, pt))));
			Vec3d bitangent = (1.0 / det) * (
				-(mesh->getTriTexCoords(tri, pt3)[0] - mesh->getTriTexCoords(tri, pt)[0]) *
				(mesh->getPt(mesh->getTri(tri, pt2)) - mesh->getPt(mesh->getTri(tri, pt))) +
				(mesh->getTriTexCoords(tri, pt2)[0] - mesh->getTriTexCoords(tri, pt)[0]) *
				(mesh->getPt(mesh->getTri(tri, pt3)) - mesh->getPt(mesh->getTri(tri, pt))));
			Vec3d normal = -mesh->getPtNormal(mesh->getTri(tri, pt));
			tangent.normalize();
			bitangent.normalize();
//			Vec3d normal = tangent ^ binormal;

			(*mats)[tri*3 + pt] = Mat3d(
				tangent[0], bitangent[0], normal[0],
				tangent[1], bitangent[1], normal[1],
				tangent[2], bitangent[2], normal[2]);
		}
	}
}

static const double NORMAL_TOL = cos(70 * DEG_TO_RAD);

void buildNormalMap(TriMesh *srcMesh, TexBary *texBary, TriMesh *targetMesh, Vec3d *normalMap, 
				   Mat3d *tangentSpace) {
	int i;

	cout << "calculating HBBs..." << endl;
	targetMesh->calcNormals();
	targetMesh->calcHBB(32);

	EdgeList *edgeList;
	char *vertList;
	edgeList = new EdgeList();
	edgeList->buildFromTriMesh(*targetMesh);
	vertList = new char[edgeList->numVerts];
	edgeList->markVerts(vertList);

	// zero out texture
	memset(normalMap, 0, NORMAL_MAP_SIZE * sizeof(Vec3d));

	cout << "finding closest points..." << endl;
	for (i=0; i < NORMAL_MAP_SIZE; i++) {
		if (!texBary[i].exists)
			continue;

		// look up vertices
		int v[3];
		v[0] = srcMesh->getTri(texBary[i].tri, 0);
		v[1] = srcMesh->getTri(texBary[i].tri, 1);
		v[2] = srcMesh->getTri(texBary[i].tri, 2);

		// calculate interpolated position
		Vec3d pos = 
			srcMesh->getPt(v[0]) * texBary[i].bary[0] +
			srcMesh->getPt(v[1]) * texBary[i].bary[1] +
			srcMesh->getPt(v[2]) * texBary[i].bary[2];
		// caculate interpolated normal
		Vec3d norm = 
			srcMesh->getPtNormal(v[0]) * texBary[i].bary[0] +
			srcMesh->getPtNormal(v[1]) * texBary[i].bary[1] +
			srcMesh->getPtNormal(v[2]) * texBary[i].bary[2];
		norm.normalize();
		// calculate interpolate tangent space
		Mat3d ts = 
			tangentSpace[texBary[i].tri*3 + 0] * texBary[i].bary[0] + 
			tangentSpace[texBary[i].tri*3 + 1] * texBary[i].bary[1] +
			tangentSpace[texBary[i].tri*3 + 2] * texBary[i].bary[2];

		// find closest point
		targetMesh->closestRestrictNormal = true;
		targetMesh->closestNormalRestriction = norm;
		Vec3d hitNorm;
		if (targetMesh->calcClosestPoint(pos, 0.05)) {
			Vec3d curClosestPt = targetMesh->closestPt;
			Vec3d delta = targetMesh->closestPt - pos;
			delta.normalize();

			// is closest point on a vertex?
			if (targetMesh->closestTri[1] == -1) {
				// normal is vertex normal
				hitNorm = targetMesh->getPtNormal(targetMesh->closestTri[0]);
				// check if this vertex adjoins a hole
				if (vertList[targetMesh->closestTri[0]]) {
					continue;
				}
			}
			// is closest point on an edge?
			else if (targetMesh->closestTri[2] == -1) {
				// interpolate the normal
				hitNorm = 
					targetMesh->closestBary[0] * targetMesh->getPtNormal(targetMesh->closestTri[0]) +
					targetMesh->closestBary[1] * targetMesh->getPtNormal(targetMesh->closestTri[1]);
				// check if this edge adjoins a hole
				EdgeInfo *ei = edgeList->findEdge(targetMesh->closestTri[0], targetMesh->closestTri[1]);
				if (!ei) {
					cout << "error: unknown edge!!" << endl;
					return;
				}
				if (ei->count == 1) {
					continue;
				}
			}
			// is closest point on a triangle?
			else {
				// interpolate the normal
				hitNorm = 
					targetMesh->closestBary[0] * targetMesh->getPtNormal(targetMesh->closestTri[0]) +
					targetMesh->closestBary[1] * targetMesh->getPtNormal(targetMesh->closestTri[1]) +
					targetMesh->closestBary[2] * targetMesh->getPtNormal(targetMesh->closestTri[2]);
			}

			hitNorm.normalize();
			if (targetMesh->closestNormalRestriction * hitNorm > NORMAL_TOL) {
				hitNorm = ts.inverse() * hitNorm;
				hitNorm.normalize();
				normalMap[i] = hitNorm;
			}
		}
/*				if (targetMesh->calcRayIntersection(pos, norm)) {
			int pt = -1;
			if (targetMesh->hitPos && targetMesh->hitNeg) {
				pt = (fabs(targetMesh->tPos) < fabs(targetMesh->tNeg)) ? 
					targetMesh->ptPos : targetMesh->ptNeg;
			}
			else if (targetMesh->hitPos) {
				pt = targetMesh->ptPos;
			}
			else if (targetMesh->hitNeg) {
				pt = targetMesh->ptNeg;
			}
			if (pt >= 0) {
				norm = targetMesh->getPtNormal(pt);
//						norm = ts * norm;
				norm.normalize();
				normalMap[y*256 + x] = norm;
			}
		}*/
	}
	cout << "done building normal map" << endl;
}

void diffuseNormalMap(Vec3d *normalMap) {
	Vec3d *normalMap2 = new Vec3d[NORMAL_MAP_SIZE];
	unsigned char *mask = new unsigned char[NORMAL_MAP_SIZE];
	int i, j, iter;

	for (i=0; i < NORMAL_MAP_SIZE; i++) {
		mask[i] = normalMap[i].iszero() ? 1 : 0;
	}

	for (iter = 0; iter < 1000; iter++) {
		memcpy(normalMap2, normalMap, NORMAL_MAP_SIZE * sizeof(Vec3d));
		for (i=0; i < NORMAL_MAP_SIZE; i++) {
			if (mask[i]) {
				int x, y, count = 0;
				Vec3d sum;
				for (y = -1; y < 2; y++) {
					for (x = -1; x < 2; x++) {
						int ind = i + x + y * 256;
						if (ind < 0 || ind > NORMAL_MAP_SIZE)
							continue;
						if (normalMap2[ind].iszero())
							continue;
						sum += normalMap2[ind];
						count++;
					}
				}
				if (count > 0) {
					sum /= count;
					normalMap[i] = sum;
				}
			}
		}
	}

	delete []mask;
	delete []normalMap2;
}

void symmetrifyNormalMap(Vec3d *normalMap) {
	int x, y, ind, ind2;
	bool hasTop, hasBot;
	Vec3d sum;

	for (y=0; y < NORMAL_MAP_H / 2; y++) {
		for (x=0; x < NORMAL_MAP_W; x++) {
			ind = y*NORMAL_MAP_W + x;
			ind2 = (NORMAL_MAP_H-y-1)*NORMAL_MAP_W + x;
			hasTop = !normalMap[ind].iszero();
			hasBot = !normalMap[ind2].iszero();
			if (!hasTop && !hasBot)
				continue;

			sum.zeroElements();
			if (hasTop)
				sum += normalMap[ind];
			if (hasBot) {
				sum[0] += normalMap[ind2][0];
				sum[1] -= normalMap[ind2][1];
				sum[2] += normalMap[ind2][2];
			}
			if (hasTop && hasBot)
				sum /= 2;

			normalMap[ind] = sum;
			normalMap[ind2] = sum;
			normalMap[ind2][1] *= -1;
		}
	}
}

void renderNormalMap(Vec3d *normalMap, unsigned char *tex, 
					 TexBary *texBary, Mat3d *tangentSpace, QuatNorm lightRot, Mat4d modelMat) {
	int i;

	modelMat = modelMat.inverse();
	Vec3d lightPos = vec4to3(modelMat * lightRot.toMatrixD() * Vec4d(0,-1,0,0));
	lightPos.normalize();
	Vec3d viewPos = vec4to3(modelMat * Vec4d(0,0,-1,0));

	// diffusely illuminate the normal map
	for (i=0; i < NORMAL_MAP_SIZE; i++) {
		// interpolate tangent space
		Mat3d ts = 
			tangentSpace[texBary[i].tri*3 + 0] * texBary[i].bary[0] + 
			tangentSpace[texBary[i].tri*3 + 1] * texBary[i].bary[1] +
			tangentSpace[texBary[i].tri*3 + 2] * texBary[i].bary[2];
		// calculate global normal
		Vec3d v = ts * normalMap[i];
		v.normalize();
		double l = max(0, v * lightPos);
		double spec = pow(max(0, v * (0.5 * (lightPos + viewPos))), 60.0);
		int intensity = min(255, 16 + l * 180 + spec * 100);
		tex[i*3+0] = intensity;
		tex[i*3+1] = intensity;
		tex[i*3+2] = intensity;
/*		double l = max(0, v * Vec3d(0.707, 0.707, 0));
		double l2 = max(0, v * Vec3d(-0.707, 0, -0.707));
		double l3 = max(0, v * Vec3d(0, 0.707, -0.707));
		tex[i*3+0] = 32 + l * 84 + l2 * 96 + l3 * 96;
		tex[i*3+1] = 32 + l * 84 + l2 * 84 + l3 * 96;
		tex[i*3+2] = 32 + l * 96 + l2 * 84 + l3 * 106;*/
	}
}