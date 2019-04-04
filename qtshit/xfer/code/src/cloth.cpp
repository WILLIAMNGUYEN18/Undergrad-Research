#include "cloth.h"
#include "ba.h"

// ClothEdge ==============================================

ClothEdge::ClothEdge() {
	v0 = 0;
	v1 = 0;
	f0 = -1;
	f1 = -1;
}

void ClothEdge::update(TriMesh *tm) {
	Vec3d v;
	if (f1 >= 0)
		update(tm->getPt(opp0), tm->getPt(v0), tm->getPt(v1), tm->getPt(opp1));
	else
		update(tm->getPt(opp0), tm->getPt(v0), tm->getPt(v1), v);
	/*e = tm->getPt(v0) - tm->getPt(v1);
	eCurLen = e.length();
	e = (1.0 / eCurLen) * e;

	if (f1 >= 0) {
		n0 = (tm->getPt(v1) - tm->getPt(opp0)) ^ (tm->getPt(v0) - tm->getPt(opp0));
		n0Len = n0.length();
		n0 = (1.0/n0Len) * n0;

		n1 = (tm->getPt(v0) - tm->getPt(opp1)) ^ (tm->getPt(v1) - tm->getPt(opp1));
		n1Len = n1.length();
		n1 = (1.0/n1Len) * n1;

		eCurAngle = atan2((n0 ^ n1) * e, n0 * n1);

		he0 = n0Len / eCurLen;
		he1 = n1Len / eCurLen;
	}*/
}

void ClothEdge::update(Vec3d &x0, Vec3d &x1, Vec3d &x2, Vec3d &x3) {
	e = x1 - x2;
	eCurLen = e.length();
	e = (1.0 / eCurLen) * e;

	if (f1 >= 0) {
		n0 = (x2 - x0) ^ (x1 - x0);
		n0Len = n0.length();
		n0 = (1.0/n0Len) * n0;

		n1 = (x1 - x3) ^ (x0 - x3);
		n1Len = n1.length();
		n1 = (1.0/n1Len) * n1;

		eCurAngle = atan2((n0 ^ n1) * e, n0 * n1);

		he0 = n0Len / eCurLen;
		he1 = n1Len / eCurLen;
	}
}

void ClothEdge::init(TriMesh *tm) {
	update(tm);
	eRestLen = eCurLen;
	eRestAngle = eCurAngle;
	restBendFactor = 6.0 * sqr(eCurLen) / (n0Len + n1Len);
}

// Cloth ==================================================

Cloth::Cloth(TriMesh *mesh) {
	int i, j, k;

	tm = mesh;

	// build edge list
	vector<vector<ClothEdge> > tempEdges;
	numEdges = 0;

	tempEdges.resize(tm->numPts());
	for (i=0; i < tm->numTris(); i++) {
		for (j=0; j < 3; j++) {
			int v0 = tm->getTri(i, j);
			int v1 = tm->getTri(i, (j+1)%3);
			int opp = tm->getTri(i, (j+2)%3);

			if (v1 < v0)
				swap(v0, v1);

			int n = (int)tempEdges[v0].size();
			for (k=0; k < n; k++) {
				if (tempEdges[v0][k].v1 == v1) {
					tempEdges[v0][k].f1 = i;
					tempEdges[v0][k].opp1 = opp;
					break;
				}
			}
			if (k == n) {
				tempEdges[v0].push_back(ClothEdge());
				tempEdges[v0][n].v0 = v0;
				tempEdges[v0][n].v1 = v1;
				tempEdges[v0][n].f0 = i;
				tempEdges[v0][n].opp0 = opp;
				numEdges++;
			}
		}
	}

	edges = new ClothEdge[numEdges];
	int index = 0;
	for (i=0; i < tm->numPts(); i++) {
		for (j=0; j < tempEdges[i].size(); j++) {
			edges[index++] = tempEdges[i][j];
		}
	}

	numVerts = tm->numPts();
	pos = new DiagonalMatrix(numVerts);
	vel = new DiagonalMatrix(numVerts);
	mass = new DiagonalMatrix(numVerts);
	totalForces = new DiagonalMatrix(numVerts);

	mP = new DiagonalMatrix(numVerts);
	mPInv = new DiagonalMatrix(numVerts);
	mz = new DiagonalMatrix(numVerts);
	mb = new DiagonalMatrix(numVerts);
	mr = new DiagonalMatrix(numVerts);
	mc = new DiagonalMatrix(numVerts);
	mq = new DiagonalMatrix(numVerts);
	ms = new DiagonalMatrix(numVerts);
	my = new DiagonalMatrix(numVerts);
	mdv = new DiagonalMatrix(numVerts);
	mTemp = new DiagonalMatrix(numVerts);
	dmat = new DiagonalMatrix(numVerts);

	df_dx = new SymMatrixBlocks(numVerts);
	df_dv = new SymMatrixBlocks(numVerts);
	mA = new SymMatrixBlocks(numVerts);
	mM = new SymMatrixBlocks(numVerts);
	tm1 = new SymMatrixBlocks(numVerts);
	tm2 = new SymMatrixBlocks(numVerts);

	mS = new Mat3d[numVerts];

    for (i=0; i < numVerts; i++) {
		mass->m_block[i][0] = 1;
		mass->m_block[i][1] = 1;
		mass->m_block[i][2] = 1;

		pos->m_block[i] = tm->getPt(i);
		vel->m_block[i] = Vec3d();

		(*(*mM)(i,i)) = Mat3d();
	}

	for (i=0; i < numEdges; i++) {
		edges[i].init(tm);
		CStretch *pStretch = new CStretch(edges[i].v0, edges[i].v1, 100, 0.001, edges[i].eCurLen, edges[i].eCurLen*3);
		forces.push_back((CForce*)pStretch);
		if (edges[i].f1 > -1) {
			CBend *pBend = new CBend(&edges[i], 0.01, 0);
			forces.push_back((CForce*)pBend);
		}
	}

	cout << numVerts << " verts; " << numEdges << " edges" << endl;
}

//#define EULER_METHOD

void Cloth::simulate() {
	int i;

	for (i=0; i < numEdges; i++) {
		edges[i].update(tm);
	}

	for (i=0; i < numVerts; i++) {
		pos->m_block[i] = tm->getPt(i);
	}

	// based on http://www.people.fas.harvard.edu/~hchong/cs276r/index.html

	// the looping maximum value that worked for Dean Macri is:
	int maxLoops = (int)sqrt((double)tm->numPts())*3 + 3;

	// step size
	double fStep = 0.03;

	// the modified method requires a couple temporary variables
	double delta0, delta_new, alpha, delta_old;
	// pre-defined epsilon -- could allow user to modify later
	double epsilon = (double)1e-22;

	// Initialize forces information
	totalForces->Zero();
	df_dx->Zero();	
	df_dv->Zero();

	// Apply forces to vertices
	for(i=0; i < forces.size(); i++)
		forces[i]->Apply(fStep, *pos, *vel, *totalForces, *mass, *df_dx, *df_dv);

#ifdef EULER_METHOD
	(*vel) += (*totalForces) * fStep;  // euler method
#else
	// Calculate: A = (M - h*df/dv - h^2*df/dx)
	// We can eliminate one float multiply each cycle by regrouping into:
	//   A = M - (df/dv + df/dx * h) * h
	(*tm1) = (*df_dx);
	(*tm1) *= fStep;
	(*tm2) = (*df_dv);
	(*tm2) += (*tm1);
	(*tm2) *= fStep;
	(*mA) = (*mM);
	(*mA) -= (*tm2);

	// Calculate: b = h * (f(0) + h * df/dx * v(0) + df/dx*y)
	// again, eliminate some compute cycles by factoring:
	// b = (f(0) + df/dx * (h*v(0) + y)) * h
	(*mTemp) = (*vel);
	(*mTemp) *= fStep;
	(*dmat) = (*mTemp);
	(*dmat) += (*my);
	(*df_dx).SetDiag(*dmat, *mTemp);
	(*dmat) = (*totalForces);
	(*dmat) += (*mTemp);
	(*mb) = (*dmat);
	(*mb) *= fStep;

	//////////////////////////////////////////
	//  Modified Conjugate Gradient Method	//
	//////////////////////////////////////////

	// Create pre-conditioner based on diagonal of A
	for(i=0; i < numVerts; i++ ) {
		mPInv->m_block[i][0] = (double)(*mA)(i,i)->n[0];
		mPInv->m_block[i][1] = (double)(*mA)(i,i)->n[4];
		mPInv->m_block[i][2] = (double)(*mA)(i,i)->n[8];
	}
	if (!mPInv->Invert(*mP)) {
		cout << "warning: uninvertable matrix!" << endl;
	}

	// change in v = z
	(*mdv) = (*mz);

	// delta0 = filter(b)^t * P * filter(b)
	for(i=0; i < numVerts; i++)
		mTemp->m_block[i] = (mS[i] * (mb->m_block[i]));
	(*mP).Mult(*mTemp, *dmat);
	delta0 = mTemp->DotProductDiagonals( *dmat );

	// r = filter( b - A * dv )
	(*mA).SetDiag(*mdv, *dmat);
	(*mTemp) = (*mb);
	(*mTemp) -= (*dmat);
	for(i=0; i < numVerts; i++)
		mr->m_block[i] = (mS[i] * (mTemp->m_block[i]));

	// c = filter( PInv * r )
	(*mPInv).Mult(*mr, *mTemp);
	for(i=0; i < numVerts; i++)
		mc->m_block[i] = (mS[i] * (mTemp->m_block[i]));

	// delta_new = r^t * c
	delta_new = mr->DotProductDiagonals( (*mc) );

	int iLoops = 0;
	while( (delta_new > epsilon * delta0) && (iLoops < maxLoops) ) {
		// q = filter( A*c )
		(*mA).SetDiag(*mc, *mTemp);
		for(i=0; i < numVerts;i++)
			mq->m_block[i] = (mS[i] * (mTemp->m_block[i]));

		// alpha = delta_new / (c^t * q)
		alpha = delta_new / (mc->DotProductDiagonals( (*mq) ) );
		
		// deltav = deltav + alpha*c
		(*mTemp) = (*mc);
		(*mTemp) *= alpha;
		(*mdv) += (*mTemp);

		// r = r - alpha * q
		(*mTemp) = (*mq);
		(*mTemp) *= alpha;
		(*mr) -= (*mTemp);

		// s = PInv * r
		(*mPInv).Mult(*mr, *ms);

		delta_old = delta_new;
		// delta_new = r^t * s
		delta_new = mr->DotProductDiagonals((*ms));

		// c = filter( s + delta_new/delta_old * c )
		(*dmat) = (*mc);
		(*dmat) *= (delta_new/delta_old);
		(*mTemp) = (*dmat);
		(*mTemp) += (*ms);
		for(i=0; i < numVerts; i++)
			mc->m_block[i] = (mS[i] * (mTemp->m_block[i]));
		
		iLoops++;
	}

	// update vertex positions and velocities
	(*vel) += (*mdv);
#endif

	(*mTemp) = (*vel) * fStep;
	(*vel) *= 0.99;
//	(*pos) += (*mTemp);
	cout << totalForces->m_block[1] << " " << mdv->m_block[1] << endl;
	int lastVert = numVerts;
	for (i=0; i < lastVert; i++)
		pos->m_block[i] += mTemp->m_block[i];
	for (; i < numVerts; i++)
		vel->m_block[i] = Vec3d();

/*
	// do collision detections
	CCollisionObj *pCObj = m_context;
	while(pCObj!= NULL)
	{
		int i;
		switch(pCObj->m_iType)
		{
		case 0:
			// sphere
			for(i=0;i<m_iNumVerts;i++)
			{
				Vec3d *pos = &(m_Pos->m_block[i]);
				Vec3d diff = (*pos) - pCObj->m_center;
				float norm = sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
				if(norm==0)
				{
					// this really shouldn't be.  
				}
				else
				{
					if(norm<pCObj->m_size[0])
					{
						bDecStepSize = true;
						// have to push the vertex outside of sphere
						float oneover = 1.0/norm;
						diff *= oneover;  // normalize vector
						Vec3d outside = diff * pCObj->m_size[0];
						m_Pos->m_block[i] = outside + pCObj->m_center;
						m_Vel->m_block[i] -= diff * (diff*m_Vel->m_block[i]);
						
						if(m_vconstraint[i]>0)
						{
							
							// change the constrained direction
							m_Constraints[m_vconstraint[i]]->SetP(&diff);
							m_Constraints[m_vconstraint[i]]->SetOwner(pCObj);
							ReAttatchConstraint(m_vconstraint[i]);
						
						}
						else
						{
							
							// create a constraint
							CConstraint *cons = new CConstraint(i, 2, diff, diff, 0);
							m_vconstraint[i] = AddConstraint(cons);
							m_Constraints[m_vconstraint[i]]->SetOwner(pCObj);
							
						}
					}
					else
					{
						if(m_vconstraint[i]>0)
						{
							m_y->m_block[i][1] = 0;
							// check that the owner of a constraint is the current object 
							if(pCObj == m_Constraints[m_vconstraint[i]]->GetOwner())
							{
								// if the force is pointed away from the object, 
								//   release constraint because point is already outside
								float dp = (m_TotalForces->m_block[i] * diff);
								if( dp > 0 )
								{
									RelaxConstraint(m_vconstraint[i]);
									bDecStepSize = false;
								}
							}
						}
					}

				}
			}
			break;
		case 1:
			break;
		default:
			break;
		}
		pCObj = pCObj->m_next;
	}
	// restrict the cloth from falling through the floor
	for( i=0; i<m_iNumVerts; i++ )
	{
		if( m_Pos->m_block[i][1] < FLOOR_Y )
		{
			// Set Constraint such that the vertices may not go any farther in the y-dir
			m_S[i].Ident();
			m_S[i].n[4] = 0;
			// set change in height
			m_y->m_block[i][1] = (FLOOR_Y - m_Pos->m_block[i][1]);
			m_Pos->m_block[i][1] = FLOOR_Y;
			m_Vel->m_block[i][1] = 0;
		}
		else
		{
			// in other cases it might be worth releasing constraint after a while
			// here, we just undo the delta position so that it stays a certain
			//   level above the ground.  
			m_y->m_block[i][1] = 0;
		}
	}*/

    for (i=0; i < numVerts; i++) {
		tm->getPt(i) = pos->m_block[i];
	}
	tm->calcNormals();
}

// vector stuff ===========================================

Mat3d createMat(Vec3d v) {
	Mat3d ret;
	ret.n[0] = v[0]*v[0];
	ret.n[1] = v[0]*v[1];
	ret.n[2] = v[0]*v[2];
	ret.n[3] = ret.n[1];
	ret.n[4] = v[1]*v[1];
	ret.n[5] = v[1]*v[2];
	ret.n[6] = ret.n[2];
	ret.n[7] = ret.n[5];
	ret.n[8] = v[2]*v[2];
	return ret;
}

Mat3d createMat(Vec3d v, Vec3d src) {
	Mat3d ret;
	ret.n[0] = v[0]*src[0];
	ret.n[1] = v[0]*src[1];
	ret.n[2] = v[0]*src[2];
	ret.n[3] = v[1]*src[0];
	ret.n[4] = v[1]*src[1];
	ret.n[5] = v[1]*src[2];
	ret.n[6] = v[2]*src[0];
	ret.n[7] = v[2]*src[1];
	ret.n[8] = v[2]*src[2];
	return ret;
}

// CStretch ===============================================

CStretch::CStretch(int v1, int v2, double c1, double c2, double rest, double max) {
	if (v1 > v2)
		swap(v1, v2);

	m_bUse = true;
	m_vert[0] = v1;
	m_vert[1] = v2;
	m_k = c1;
	m_kDamp = c2;
	m_nat = rest;
	m_max = max;
}

void CStretch::Apply( float fTime, DiagonalMatrix &pos, DiagonalMatrix &vel, 
		DiagonalMatrix &totalforce, DiagonalMatrix &mass,
		SymMatrixBlocks &df_dx, SymMatrixBlocks &df_dv ) {
	Vec3d	force;							// force on vertices
	Vec3d	deltax;							// difference in positions between vertices
	Vec3d	dC_dx[2];						// dC/dx vectors for each vertex
	double		C;								// C <-- corresponds to positional stretch beyond rest length
	double		norm;							// norm of the change in position between vertices
	double		C_dot;							// dot-product
	Mat3d	dCdx_dCdy;						// dC/dx * (dC/dy)^t
	Mat3d	d2C_dxdy[2][2];					// 3D matrix of double derivatives
												// Each vector represends a diagonal 3x3 matrix
	int i, j;


	if(!m_bUse) return;

	//////////////////////////////////////////////////////////////////////////////////
	// A L G O R I T H M /////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	//
	//	This routine applies a stretch force between two vertices specified in
	//  initialization by indicies m_vert array.
	//  The formulation in Baraff and Witkin's large steps in cloth simulation
	//  paper uses the quadratic 1/2 * k * x^2 energy function for stretch 
	//  as opposed to the quartic function used by some because quadratic suffices
	//  and prevents "too much" stiffness.
	//  The vector quantity C is defined to be ||pos2-pos1|| - constant.
	//  More specifically:
	//    C = [(x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2]^(1/2) - natural_length
	//  The derivative of this is simply the unit vector that points in the direction
	//     of change.  C itself is the difference between positions minus the rest
	//     length of the "spring" between the vertices.  This gives the length by which
	//     the spring has been stretched/compressed.
	//  Note that since C does not depend on velocity, df/dv will be zero.
	//  The 2nd derivatives are stored with vector information, giving a 3D matrix
	//     d2C_dxdy that when multiplied by a standard vector will give a 3x3 matrix.
	//
	//  Damping is also handled here.  Damping is meant to resist cloth stretching/compression.
	//     Therefore, it is a function of the difference in velocities of each vertex.
	//     If the two vertices are moving the same direction with same magnitude, 
	//	   then there is no damping needed.
	//     If they move in opposite directions, then we need to apply a damping term.
	//	   This reduces to computing:
	//		dC/dp * (v2 - v1) for C_dot
	//
	/////////////////////////////////////////////////////////////////////////////////////////

	// calculate position data
	deltax = (pos.m_block[m_vert[1]] - pos.m_block[m_vert[0]]);
	norm = sqrt(deltax[0]*deltax[0]+deltax[1]*deltax[1]+deltax[2]*deltax[2]);
	C = norm - m_nat;

	// we can't have norm zero; if it's zero, give it some epsilon
	if(norm==0) norm = (double)1e-22;
	// divisions are more expensive than multiplies.
	double OneOverNorm = 1.0/norm;
	double OneOverCubed = OneOverNorm*OneOverNorm*OneOverNorm;
	// calculate dC/dx. one vertex's is negative the other's
	dC_dx[0] = (deltax*(-OneOverNorm));
	dC_dx[1] = (deltax*OneOverNorm);

	// remember dC_dx[0] has the negative sign, so this calculates
	//    dC_dx * (v2 - v1)
	C_dot = dC_dx[0] * vel.m_block[m_vert[0]] + 
			dC_dx[1] * vel.m_block[m_vert[1]];

	// spring force = - k * (displacement from rest) * (unit direction vector)
	// damping force = - k_Damp * direction * (velocity difference)

	// vertex 1
	force = -m_k * C * dC_dx[0] - m_kDamp * dC_dx[0] * C_dot;
	// add these to the cloth's total internal forces
	totalforce.m_block[m_vert[0]] += force;

	// vertex 2
	force = -m_k * C * dC_dx[1] - m_kDamp * dC_dx[1] * C_dot;
	// add contribution to internal forces
	totalforce.m_block[m_vert[1]] += force;


	// calculate second derivatives

	// d/dx1(d/dx1(C)) = -[(x2-x1)^2+(y2-y1)^2+(z2-z1)^2]^(-3/2) * (x2-x1)^2 + []^(-1/2) 
	// d/dx2(d/dx2(C)) = -[(x2-x1)^2+(y2-y1)^2+(z2-z1)^2]^(-3/2) * (x2-x1)^2 + []^(-1/2) 
	// d/dx2(d/dx1(C)) = [(x2-x1)^2+(y2-y1)^2+(z2-z1)^2]^(-3/2) * (x2-x1)^2 - []^(-1/2)
	// d/dx1(d/dx2(C)) = [(x2-x1)^2+(y2-y1)^2+(z2-z1)^2]^(-3/2) * (x2-x1)^2 - []^(-1/2)
	for (i=0; i < 3; i++)
		for (j=0; j < 3; j++) {
			if (i == j) {
				d2C_dxdy[0][0][i][j] = -OneOverCubed * ( deltax[i] * deltax[j] ) + OneOverNorm;
				d2C_dxdy[0][1][i][j] =  OneOverCubed * ( deltax[i] * deltax[j] ) - OneOverNorm;
			}
			else {
				d2C_dxdy[0][0][i][j] = -OneOverCubed * ( deltax[i] * deltax[j] );
				d2C_dxdy[0][1][i][j] =  OneOverCubed * ( deltax[i] * deltax[j] );
			}
		}
	d2C_dxdy[1][1] = d2C_dxdy[0][0]; 
	d2C_dxdy[1][0] = d2C_dxdy[0][1];


	Mat3d dfdx, dfdv, dCdC;

	for( i=0; i<2; i++ ) {
		// we're storing upper triangular only, so start at index i
		for( j=i; j<2; j++ ) {
			// Kij = -k (dC/dxi * (dC/dxj)^t + d2C/dxidxj * C)
			//  there's also the damping contribution: - k_Damp * (d2C/dxidxj * C_dot)
			//  Note: the damping term drops a term because that term is non-symmetric
			//   and would ruin our modified conjugate gradient method.  
			//   Baraff and Witkin assure us that there will be no "ill effects."
			// Making use of the fact that the 2nd derivatives are diagonal (very sparse):
			dCdC = createMat(dC_dx[i], dC_dx[j]);
			dfdx = -m_k * (dCdC + d2C_dxdy[i][j] * C) - m_kDamp * d2C_dxdy[i][j] * C_dot;

			// add contribution to overall force matrix
			(*(df_dx(m_vert[i], m_vert[j]))) +=  dfdx;

			// now compute the velocity derivatives.
			// remember, the spring energy function does not depend on velocity and therefore
			//   has no contribution here.  only the damping term factors in.
			// NOTE: we make use of the following derivation:
			//  C' = dC/dt = dC/dx * dx/dt = dC/dx * v
			//  dC'/dv = d/dv (dC/dx * v) = dC/dx
			// in our matrix notation, dC'/dv = d/dv ((dC/dx)^t * v) = dC/dx
			dfdv = -m_kDamp * dCdC;

			// add contribution to overall force derivative matrix
			(*(df_dv(m_vert[i], m_vert[j]))) += dfdv;
		}
	}
}


// CAttract ===============================================

CAttract::CAttract(int v1, Vec3d attractPt, double c1, double c2, double rest) {
	m_attract = attractPt;

	m_bUse = true;
	m_vert = v1;
	m_k = c1;
	m_kDamp = c2;
	m_nat = rest;
}

void   CAttract::Apply( float fTime, DiagonalMatrix &pos, DiagonalMatrix &vel, 
						DiagonalMatrix &totalforce, DiagonalMatrix &mass,
						SymMatrixBlocks &df_dx, SymMatrixBlocks &df_dv ) {
	Vec3d	*vel1;       					// velocities for vertices of interest
	Vec3d	force;							// force on vertices
	Vec3d	deltax;							// difference in positions between vertices
	Vec3d	dC_dx;							// dC/dx vectors for each vertex
	double		C;								// C <-- corresponds to positional stretch beyond rest length
	double		norm;							// norm of the change in position between vertices
	double		C_dot;							// dot-product
	Mat3d	dCdx_dCdy;						// dC/dx * (dC/dy)^t
	Mat3d	d2C_dxdy;					// 3D matrix of double derivatives
												// Each vector represends a diagonal 3x3 matrix
	int i, j;


	if(!m_bUse) return;

	//////////////////////////////////////////////////////////////////////////////////
	// A L G O R I T H M /////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	//
	//	This routine applies a stretch force between two vertices specified in
	//  initialization by indicies m_vert array.
	//  The formulation in Baraff and Witkin's large steps in cloth simulation
	//  paper uses the quadratic 1/2 * k * x^2 energy function for stretch 
	//  as opposed to the quartic function used by some because quadratic suffices
	//  and prevents "too much" stiffness.
	//  The vector quantity C is defined to be ||pos2-pos1|| - constant.
	//  More specifically:
	//    C = [(x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2]^(1/2) - natural_length
	//  The derivative of this is simply the unit vector that points in the direction
	//     of change.  C itself is the difference between positions minus the rest
	//     length of the "spring" between the vertices.  This gives the length by which
	//     the spring has been stretched/compressed.
	//  Note that since C does not depend on velocity, df/dv will be zero.
	//  The 2nd derivatives are stored with vector information, giving a 3D matrix
	//     d2C_dxdy that when multiplied by a standard vector will give a 3x3 matrix.
	//
	//  Damping is also handled here.  Damping is meant to resist cloth stretching/compression.
	//     Therefore, it is a function of the difference in velocities of each vertex.
	//     If the two vertices are moving the same direction with same magnitude, 
	//	   then there is no damping needed.
	//     If they move in opposite directions, then we need to apply a damping term.
	//	   This reduces to computing:
	//		dC/dp * (v2 - v1) for C_dot
	//
	/////////////////////////////////////////////////////////////////////////////////////////

	// point to relavant vertices
	vel1 = &vel.m_block[m_vert];

	// calculate position data
	deltax = (pos.m_block[m_vert] - m_attract);
	norm = sqrt(deltax[0]*deltax[0]+deltax[1]*deltax[1]+deltax[2]*deltax[2]);
	C = norm - m_nat;

	// we can't have norm zero; if it's zero, give it some epsilon
	if(norm==0) norm = (double)1e-22;
	// divisions are more expensive than multiplies.
	double OneOverNorm = 1.0/norm;
	double OneOverCubed = OneOverNorm*OneOverNorm*OneOverNorm;
	// calculate dC/dx. one vertex's is negative the other's
	dC_dx = (deltax*OneOverNorm);

	// remember dC_dx[0] has the negative sign, so this calculates
	//    dC_dx * (v2 - v1)
	C_dot = dC_dx[0] * (*vel1)[0] + 
			dC_dx[1] * (*vel1)[1] + 
			dC_dx[2] * (*vel1)[2];

	// spring force = - k * (displacement from rest) * (unit direction vector)
	// damping force = - k_Damp * direction * (velocity difference)

	// vertex 1
	force = -m_k * C * dC_dx - m_kDamp * dC_dx * C_dot;
	// add these to the cloth's total internal forces
	totalforce.m_block[m_vert] += force;

	// calculate second derivatives

	// d/dx1(d/dx1(C)) = -[(x2-x1)^2+(y2-y1)^2+(z2-z1)^2]^(-3/2) * (x2-x1)^2 + []^(-1/2) 
	// d/dx2(d/dx2(C)) = -[(x2-x1)^2+(y2-y1)^2+(z2-z1)^2]^(-3/2) * (x2-x1)^2 + []^(-1/2) 
	for (i=0; i < 3; i++)
		for (j=0; j < 3; j++) {
			if (i == j) {
				d2C_dxdy[i][j] = -OneOverCubed * ( deltax[i] * deltax[j] ) + OneOverNorm;
			}
			else {
				d2C_dxdy[i][j] = -OneOverCubed * ( deltax[i] * deltax[j] );
			}
		}
	
	Mat3d dfdx, dfdv, dCdC;

	// Kij = -k (dC/dxi * (dC/dxj)^t + d2C/dxidxj * C)
	//  there's also the damping contribution: - k_Damp * (d2C/dxidxj * C_dot)
	//  Note: the damping term drops a term because that term is non-symmetric
	//   and would ruin our modified conjugate gradient method.  
	//   Baraff and Witkin assure us that there will be no "ill effects."
	// Making use of the fact that the 2nd derivatives are diagonal (very sparse):
	dCdC = createMat(dC_dx, dC_dx);
	dfdx = -m_k * (dCdC + d2C_dxdy * C) - m_kDamp * d2C_dxdy * C_dot;
	
	// add contribution to overall force matrix
	(*(df_dx(m_vert, m_vert))) +=  dfdx;

	// now compute the velocity derivatives.
	// remember, the spring energy function does not depend on velocity and therefore
	//   has no contribution here.  only the damping term factors in.
	// NOTE: we make use of the following derivation:
	//  C' = dC/dt = dC/dx * dx/dt = dC/dx * v
	//  dC'/dv = d/dv (dC/dx * v) = dC/dx
	// in our matrix notation, dC'/dv = d/dv ((dC/dx)^t * v) = dC/dx
	dfdv = -m_kDamp * dCdC;

	// add contribution to overall force derivative matrix
	(*(df_dv(m_vert, m_vert))) += dfdv;
}


// CBend ==================================================

CBend::CBend(ClothEdge *myEdge, double c1, double c2) {
	edge = myEdge;
	m_k = c1;
	m_kDamp = c2;
}

void CBend::Apply( float fTime, DiagonalMatrix &pos, DiagonalMatrix &vel, 
						DiagonalMatrix &totalforce, DiagonalMatrix &mass,
						SymMatrixBlocks &df_dx, SymMatrixBlocks &df_dv ) {
	if(!m_bUse) return;

	int i, j, s, t;
	double C;			// difference from rest angle
	Vec3d	dC_dx[4];
	Vec3d   v[4];
	int     vert[4];
	Mat3d   dfdx;

	// store vertex positions and indices for easy reference
	vert[0] = edge->opp0;
	vert[1] = edge->v0;
	vert[2] = edge->v1;
	vert[3] = edge->opp1;
	for (i=0; i < 4; i++)
		v[i] = pos.m_block[vert[i]];

	// calculate delta
	C = (edge->eCurAngle - edge->eRestAngle);
	bool smallC = false;//(fabs(C) < 1e-4);

	// calculate dC/dx
	double cosTheta = (edge->n0 * edge->n1); //cos(edge->eCurAngle);
	double sinTheta = (edge->n0 ^ edge->n1) * edge->e; //sin(edge->eCurAngle);
	Vec3d qA[4], qB[4];
	qA[0] = v[2]-v[1]; qA[1] = v[0]-v[2]; qA[2] = v[1]-v[0]; qA[3] = 0;
	qB[0] = 0;         qB[1] = v[2]-v[3]; qB[2] = v[3]-v[1]; qB[3] = v[1]-v[2];
	double e[4];
	e[0] = 0; e[1] = 1; e[2] = -1; e[3] = 0;
	Vec3d nAGrad[4][3], nBGrad[4][3], eGrad[4][3];
	Vec3d cosGrad[4], sinGrad[4]; //, eLenGrad[4];
	Vec3d Sqa[4][3], Sqb[4][3], Sqe[4][3];
	Mat3d m0, m1, me;
	for (i=0; i < 3; i++)
		for (j=0; j < 3; j++) {
			m0[i][j] -= edge->n0[i] * edge->n0[j];
			m1[i][j] -= edge->n1[i] * edge->n1[j];
			me[i][j] -= edge->e[i] * edge->e[j];
		}
	for (i=0; i < 4; i++) {
		for (j=0; j < 3; j++) {
			switch (j) {
			case 0:
				Sqa[i][j] = Vec3d(0, -qA[i][2], qA[i][1]);
				Sqb[i][j] = Vec3d(0, -qB[i][2], qB[i][1]);
				Sqe[i][j] = Vec3d(e[i], 0, 0);
				break;
			case 1:
				Sqa[i][j] = Vec3d(qA[i][2], 0, -qA[i][0]);
				Sqb[i][j] = Vec3d(qB[i][2], 0, -qB[i][0]);
				Sqe[i][j] = Vec3d(0, e[i], 0);
				break;
			case 2:
				Sqa[i][j] = Vec3d(-qA[i][1], qA[i][0], 0);
				Sqb[i][j] = Vec3d(-qB[i][1], qB[i][0], 0);
				Sqe[i][j] = Vec3d(0, 0, e[i]);
			}
			nAGrad[i][j] = m0 * ((1.0/edge->n0Len) * Sqa[i][j]);
			nBGrad[i][j] = m1 * ((1.0/edge->n1Len) * Sqb[i][j]);
			eGrad[i][j] = me * ((1.0/edge->eCurLen) * Sqe[i][j]);
//			eLenGrad[i][j] = 2.0 * edge->e[j]*edge->eCurLen*e[j];
			cosGrad[i][j] = nAGrad[i][j] * edge->n1 + edge->n0 * nBGrad[i][j];
			sinGrad[i][j] = ((nAGrad[i][j] ^ edge->n1) + (edge->n0 ^ nBGrad[i][j])) * edge->e + (edge->n0 ^ edge->n1) * eGrad[i][j];
			if (smallC)
				dC_dx[i][j] = 0;
			else
				dC_dx[i][j] = cosTheta * sinGrad[i][j] - sinTheta * cosGrad[i][j];
		}
	}

	// calculate C_dot
	double C_dot = 0;
	for (i=0; i < 4; i++) {
		Vec3d nv = dC_dx[i];
		if (nv.length() < 1e6)
			nv = Vec3d();
		else
			nv.normalize();
		C_dot += nv * vel.m_block[vert[i]];
	}

	// accumulate forces
	for (i=0; i < 4; i++)
		totalforce.m_block[vert[i]] += -m_k * C * dC_dx[i] - m_kDamp * dC_dx[i] * C_dot;

	// calculate second derivatives
	for (i=0; i < 4; i++) {
		for (j=0; j < 4; j++) {
			// upper triangular part only
			if (vert[i] > vert[j])
				continue;

			Mat3d d2C_dxdy = Mat3d(0,0,0,0,0,0,0,0,0);

			double dQA[4][4] = {
				{0, -1, 1, 0},
				{1, 0, 0, -1},
				{-1, 1, 0, 0},
				{0, 0, 0, 0}
			};
			double dQB[4][4] = {
				{0, 0, 0, 0},
				{0, 0, 1, -1},
				{0, -1, 0, 1},
				{0, 1, -1, 0}
			};

			for (s=0; s < 3; s++) {
				for (t=0; t < 3; t++) {
					Vec3d d2nA_dxdy, d2nB_dxdy, d2e_dxdy;
					double d2sin_dxdy, d2cos_dxdy;

					d2nA_dxdy = Vec3d(0, 0, 0);
					d2nB_dxdy = Vec3d(0, 0, 0);
					Vec3d dqma_dxnt = Vec3d(0,0,0); dqma_dxnt[t] = dQA[i][j];
					Vec3d dqmb_dxnt = Vec3d(0,0,0); dqmb_dxnt[t] = dQB[i][j];
					Vec3d Sdqa_dx, Sdqb_dx;
					switch (s) {
					case 0:
						Sdqa_dx = Vec3d(0, -dqma_dxnt[2], dqma_dxnt[1]);
						Sdqb_dx = Vec3d(0, -dqmb_dxnt[2], dqmb_dxnt[1]);
						break;
					case 1:
						Sdqa_dx = Vec3d(dqma_dxnt[2], 0, -dqma_dxnt[0]);
						Sdqb_dx = Vec3d(dqmb_dxnt[2], 0, -dqmb_dxnt[0]);
						break;
					case 2:
						Sdqa_dx = Vec3d(-dqma_dxnt[1], dqma_dxnt[0], 0);
						Sdqb_dx = Vec3d(-dqmb_dxnt[1], dqmb_dxnt[0], 0);
						break;
					}
					d2nA_dxdy = 
						-(1.0 / edge->n0Len) *
							(createMat(nAGrad[j][t], edge->n0) + createMat(edge->n0, nAGrad[i][s])) 
							* Sqa[j][s] +
						(1.0/sqr(edge->n0Len)) * 
							m0 * 
							(edge->n0Len * Sdqa_dx - (Sqa[j][t]*edge->n0) * Sqa[i][s]);
					d2nB_dxdy = 
						-(1.0 / edge->n1Len) * 
							(createMat(nBGrad[j][t], edge->n1) + createMat(edge->n1, nBGrad[i][s])) 
							* Sqb[j][s] +
						(1.0/sqr(edge->n1Len)) * 
							m1 * 
							(edge->n1Len * Sdqb_dx - (Sqb[j][t]*edge->n1) * Sqb[i][s]);
/*
					// check w/finite diff:
					Vec3d oldNorm, newNorm, oldDeriv, newDeriv;
					oldNorm = (v[2] - v[0]) ^ (v[1] - v[0]);
					oldNorm.normalize();
					v[i][s] += 0.001;
					newNorm = (v[2] - v[0]) ^ (v[1] - v[0]);
					newNorm.normalize();
					v[i][s] -= 0.001;
					oldDeriv = (1.0/0.001)*(newNorm-oldNorm);
					v[j][t] += 0.001;
					oldNorm = (v[2] - v[0]) ^ (v[1] - v[0]);
					oldNorm.normalize();
					v[i][s] += 0.001;
					newNorm = (v[2] - v[0]) ^ (v[1] - v[0]);
					newNorm.normalize();
					v[i][s] -= 0.001;
					newDeriv = (1.0/0.001)*(newNorm-oldNorm);
					v[j][t] -= 0.001;
					d2nA_dxdy = (1.0/0.001)*(newDeriv-oldDeriv);
//					cout << i<<j<<s<<t<< ": " << d2nA_dxdy << " -- " << (1.0/0.001)*(newDeriv-oldDeriv) << endl;


					oldNorm = (v[1] - v[3]) ^ (v[2] - v[3]);
					oldNorm.normalize();
					v[i][s] += 0.001;
					newNorm = (v[1] - v[3]) ^ (v[2] - v[3]);
					newNorm.normalize();
					v[i][s] -= 0.001;
					oldDeriv = (1.0/0.001)*(newNorm-oldNorm);
					v[j][t] += 0.001;
					oldNorm = (v[1] - v[3]) ^ (v[2] - v[3]);
					oldNorm.normalize();
					v[i][s] += 0.001;
					newNorm = (v[1] - v[3]) ^ (v[2] - v[3]);
					newNorm.normalize();
					v[i][s] -= 0.001;
					newDeriv = (1.0/0.001)*(newNorm-oldNorm);
					v[j][t] -= 0.001;
					d2nB_dxdy = (1.0/0.001)*(newDeriv-oldDeriv);
*/
					/*
					switch (s) {
					case 0:
						if (t == 2) {
							d2nA_dxdy[1] = -(1.0/edge->n0Len)*dQA[i][j];
							d2nB_dxdy[1] = -(1.0/edge->n1Len)*dQB[i][j];
						}
						if (t == 1) {
							d2nA_dxdy[2] = (1.0/edge->n0Len)*dQA[i][j];
							d2nB_dxdy[2] = (1.0/edge->n1Len)*dQB[i][j];
						}
						break;
					case 1:
						if (t == 2) {
							d2nA_dxdy[0] = (1.0/edge->n0Len)*dQA[i][j];
							d2nB_dxdy[0] = (1.0/edge->n1Len)*dQB[i][j];
						}
						if (t == 0) {
							d2nA_dxdy[2] = -(1.0/edge->n0Len)*dQA[i][j];
							d2nB_dxdy[2] = -(1.0/edge->n1Len)*dQB[i][j];
						}
						break;
					case 2:
						if (t == 1) {
							d2nA_dxdy[0] = -(1.0/edge->n0Len)*dQA[i][j];
							d2nB_dxdy[0] = -(1.0/edge->n1Len)*dQB[i][j];
						}
						if (t == 0) {
							d2nA_dxdy[1] = (1.0/edge->n0Len)*dQA[i][j];
							d2nB_dxdy[1] = (1.0/edge->n1Len)*dQB[i][j];
						}
					}*/
					d2e_dxdy = Vec3d(0, 0, 0);

					d2cos_dxdy = 
						d2nA_dxdy * edge->n1 + nBGrad[j][t]*nAGrad[i][s] + 
						nAGrad[j][t]*nBGrad[i][s] + edge->n0*d2nB_dxdy;
					d2sin_dxdy = 
						(d2nA_dxdy ^ edge->n1 + nAGrad[i][s] ^ nBGrad[j][t] + 
						 nAGrad[j][t] ^ nBGrad[i][s] + edge->n0 ^ d2nB_dxdy) * edge->e +
						(nAGrad[i][s] ^ edge->n1 + edge->n0 ^ nBGrad[i][s]) * eGrad[j][t] +
						(nAGrad[j][t] ^ edge->n1 + edge->n0 ^ nBGrad[j][t]) * eGrad[i][s] +
						(edge->n0 ^ edge->n1) * d2e_dxdy;

					// check w/finite diff:
					double oldAngle, newAngle, oldDeriv, newDeriv;
					ClothEdge tempE = *edge;
					oldAngle = tempE.eCurAngle;
					v[i][s] += 0.001;
					tempE.update(v[0], v[1], v[2], v[3]);
					newAngle = tempE.eCurAngle;
					v[i][s] -= 0.001;
					oldDeriv = (1.0/0.001)*(newAngle-oldAngle);
					v[j][t] += 0.001;
					tempE.update(v[0], v[1], v[2], v[3]);
					oldAngle = tempE.eCurAngle;
					v[i][s] += 0.001;
					tempE.update(v[0], v[1], v[2], v[3]);
					newAngle = tempE.eCurAngle;
					v[i][s] -= 0.001;
					newDeriv = (1.0/0.001)*(newAngle-oldAngle);
					v[j][t] -= 0.001;
					d2C_dxdy[s][t] = (1.0/0.001)*(newDeriv-oldDeriv);
//					cout << newDeriv << " " << oldDeriv << endl;
//					cout << i<<j<<s<<t<< ": " << d2nA_dxdy << " -- " << (1.0/0.001)*(newDeriv-oldDeriv) << endl;


//					if (smallC)
//						d2C_dxdy[s][t] = 0;
//					else
//						d2C_dxdy[s][t] = cosGrad[j][t] * sinGrad[i][s] + cosTheta * d2sin_dxdy - sinGrad[j][t]*cosGrad[i][s] - sinTheta * d2cos_dxdy;
				}
			}

			// Kij = -k (dC/dxi * (dC/dxj)^t + d2C/dxidxj * C)
			dfdx = -m_k * (createMat(dC_dx[i], dC_dx[j]) + C * d2C_dxdy)
					-m_kDamp * C_dot * d2C_dxdy; //Mat3d(d2C_dxdy[0][0],0,0,0,d2C_dxdy[1][1],0,0,0,d2C_dxdy[2][2]);
//			cout << d2C_dxdy << endl;

			// add contribution to overall force matrix
			(*(df_dx(vert[i], vert[j]))) +=  dfdx;

			// damping
			Vec3d temp0 = dC_dx[i], temp1 = dC_dx[j];
			if (temp0.length() > 1e-6)
				temp0.normalize();
			if (temp1.length() > 1e-6)
				temp1.normalize();
			(*(df_dv(vert[i], vert[j]))) += -m_kDamp * createMat(temp0, temp1);
		}
	}
//	exit(0);
}


// DiagonalMatrix =========================================

DiagonalMatrix::DiagonalMatrix( int iSize ) {
	if( iSize > 0 ) {
		m_iSize = iSize;
		m_block = new Vec3d[m_iSize]; 
	}
	else {
		// this should never happen
		m_iSize = 0;
		m_block = NULL;
	}
}

DiagonalMatrix::DiagonalMatrix( const DiagonalMatrix &src ) {
	// initialize the array
	m_block = new Vec3d[src.m_iSize];
	m_iSize = src.m_iSize;

	// copy the values over
	if( src.m_iSize > 0 )
		memcpy( m_block, src.m_block, m_iSize * sizeof( Vec3d ) );
}

DiagonalMatrix::~DiagonalMatrix()
{
	if( m_block )
	{
		delete [] m_block;
		m_block = NULL;
	}
}

DiagonalMatrix& DiagonalMatrix::operator=( const DiagonalMatrix &src )
{
	if( m_iSize != src.m_iSize )
	{
		if( m_block )
			delete [] m_block;
		m_block = new Vec3d[src.m_iSize];
		m_iSize = src.m_iSize;
	}
	if( src.m_iSize > 0 )
	{
		memcpy( m_block, src.m_block, m_iSize * sizeof( Vec3d ) );
	}
	return *this;
}

void DiagonalMatrix::Zero()
{
	if( m_iSize && m_block )
		memset( m_block, 0, m_iSize * sizeof( Vec3d ) );
}

DiagonalMatrix& DiagonalMatrix::operator+=( const DiagonalMatrix &src )
{
	// ideally, these two diagonal matrices should be of the same size...
	//  but otherwise, there is a "natural" definition that will be allowed
	int minimum = m_iSize;
	if(src.m_iSize < minimum)
		minimum = src.m_iSize;
	for(int index=0;index<minimum;index++)
	{
		m_block[index] += src.m_block[index];
	}
	return *this;
}

DiagonalMatrix& DiagonalMatrix::operator-=( const DiagonalMatrix &src )
{
	// ideally, these two diagonal matrices should be of the same size...
	//  but otherwise, there is a "natural" definition that will be allowed
	int minimum = m_iSize;
	if(src.m_iSize < minimum)
		minimum = src.m_iSize;
	for(int index=0;index<minimum;index++)
	{
		m_block[index] -= src.m_block[index];
	}
	return *this;
}

DiagonalMatrix& DiagonalMatrix::operator*=( double scale )
{
	for(int index=0;index<m_iSize;index++)
		m_block[index] *= scale;
	return *this;
}

bool DiagonalMatrix::Mult( DiagonalMatrix &src , DiagonalMatrix &dst )
{
	if(dst.m_iSize != m_iSize) return false;
	for(int index=0;index<m_iSize;index++)
	{
		dst.m_block[index][0] = (m_block[index][0] * src.m_block[index][0]);
		dst.m_block[index][1] = (m_block[index][1] * src.m_block[index][1]);
		dst.m_block[index][2] = (m_block[index][2] * src.m_block[index][2]);
	}
	return true;
}

// helper function for checking whether a matrix-matrix op would be 
//   well-defined in the usual sense; just checks sizes of matrices
bool well_defined_op(DiagonalMatrix &m1, DiagonalMatrix &m2)
{
	if(m1.m_iSize != m2.m_iSize )
		return false;
	return true;
}

bool DiagonalMatrix::Invert(DiagonalMatrix &ret) {
	// inversion of a diagonal matrix is simply the inversion of each
	//   element. however, zero element entries make for a singular
	//   matrix and is therefore inversion is not well-defined. 
	// singular matrices will just return themselves as inverses
	// (user's responsibility to not invert singular matrices)
	if(ret.m_iSize != m_iSize) return false;
	bool nonzero = true;
	for( int i=0; i<m_iSize; i++ ) {
		double bx = m_block[i][0];
		double by = m_block[i][1];
		double bz = m_block[i][2];
		if((bx!=0)&&(by!=0)&&(bz!=0)) {
			ret.m_block[i][0] = 1.0 / bx;
			ret.m_block[i][1] = 1.0 / by;
			ret.m_block[i][2] = 1.0 / bz;
		}
		else
			nonzero = false;
	}
	if(!nonzero) {
		// make singular matrices return themselves
		return false;
	}
	return true;
}

// this treates each diagonal matrix as a vector and performs the dot product
double DiagonalMatrix::DotProductDiagonals( const DiagonalMatrix &src )
{
	// should be equal size, but if not, do best we can.
	int minimum = m_iSize;
	if(src.m_iSize<minimum) minimum = src.m_iSize;
	double ret = 0.0;

	// perform dot product on each vector element
	for( int i=0; i<minimum; i++ )
		ret += (m_block[i] * src.m_block[i]); 

	return ret;
}


// SymMatrixBlocks ========================================

SymMatrixBlocks::SymMatrixBlocks(int size) {
	if (size > 0) {
		m_iSize = size;
		// allocate one more than size so we can store the 
		//  "length" of each row the same way
		m_rowstart = new int[size+1];
		for(int i=0;i<=size;i++)
			m_rowstart[i] = 0;
	}
	else {
		m_iSize = 0;
		m_rowstart = NULL;
	}
	m_bAllZero = true;
}

SymMatrixBlocks::SymMatrixBlocks( const SymMatrixBlocks &src ) {
	// some dummy init just in case...
	int i;
	m_iSize = 0;
	m_rowstart = NULL;
	m_col_lookup.clear();
	m_matBlock.clear();
	m_bAllZero = true;

	// allocate for new information
	if (src.m_iSize > 0) {
		m_iSize = src.m_iSize;
		m_rowstart = new int[m_iSize+1];
		if(src.m_matBlock.size()>0) m_bAllZero = false;

		// now copy over new information
		for(i=0;i<=m_iSize;i++)
			m_rowstart[i] = src.m_rowstart[i];
		for(i=0;i<src.m_col_lookup.size();i++)
			m_col_lookup.push_back(src.m_col_lookup[i]);
		m_bAllZero = src.m_bAllZero;
		for(i=0;i<src.m_matBlock.size();i++) {
			Mat3d *temp = new Mat3d;
			(*temp) = *(src.m_matBlock[i]);
			m_matBlock.push_back(temp);
		}
	}
	else {
		m_iSize = 0;
		m_rowstart = NULL;
	}
}

SymMatrixBlocks::~SymMatrixBlocks() {
	m_iSize = 0;
	if(m_rowstart) {
		delete [] m_rowstart;
		m_rowstart = NULL;
	}
	m_col_lookup.erase(m_col_lookup.begin(),m_col_lookup.end());
	for(int i=0; i < m_matBlock.size(); i++) {
		if(m_matBlock[i])
			delete m_matBlock[i];
		m_matBlock.erase(m_matBlock.begin(),m_matBlock.end());
	}
	m_bAllZero = true;
}

SymMatrixBlocks& SymMatrixBlocks::operator=(const SymMatrixBlocks &src) {
	int i;

	if (src.m_iSize != m_iSize) {
		// cleanup what's currently stored here
		int i;
		m_iSize = 0;
		if(m_rowstart) {
			delete [] m_rowstart;
			m_rowstart = NULL;
		}
		m_col_lookup.erase(m_col_lookup.begin(),m_col_lookup.end());
		for(i=0; i<m_matBlock.size(); i++) {
			if(m_matBlock[i])
				delete m_matBlock[i];
			m_matBlock.erase(m_matBlock.begin(),m_matBlock.end());
		}
		m_bAllZero = true;

		// allocate for new information
		if (src.m_iSize > 0) {
			m_iSize = src.m_iSize;
			m_rowstart = new int[m_iSize+1];
			if(src.m_matBlock.size()>0) m_bAllZero = false;

			// now copy over new information
			for(i=0;i<=m_iSize;i++)
				m_rowstart[i] = src.m_rowstart[i];
			for(i=0;i<src.m_col_lookup.size();i++)
				m_col_lookup.push_back(src.m_col_lookup[i]);
			m_bAllZero = src.m_bAllZero;
			for(i=0;i<src.m_matBlock.size();i++) {
				Mat3d *temp = new Mat3d;
				(*temp) = *(src.m_matBlock[i]);
				m_matBlock.push_back(temp);
			}

		}
		else {
			m_iSize = 0;
			m_rowstart = NULL;
		}
	}
	else {
		// non-destructive copy
		int j, start, length;

		for(i=0; i < m_matBlock.size(); i++) {
			if(m_matBlock[i])
				m_matBlock[i]->zero();
		}

		for(i=0; i < src.m_iSize; i++) {
			start = src.m_rowstart[i];
			length = src.m_rowstart[i+1] - start;
			for(j = 0; j<length; j++)
				(*((*this)(i, src.m_col_lookup[start+j] ))) = (*(src.m_matBlock[start+j]));
		}
	}

	return *this;

}

void SymMatrixBlocks::Zero(void) {
	// modified to keep sparsity arrangement
	int i;
	for(i=0;i<m_matBlock.size();i++) {
		if(m_matBlock[i])
			m_matBlock[i]->zero();
	}
	m_bAllZero = true;
	/*
	int i;
	for(i=0;i<=m_iSize;i++)
		m_rowstart[i] = 0;
	m_col_lookup.erase(m_col_lookup.begin(),m_col_lookup.end());
	for(i=0;i<m_matBlock.size();i++) {
		if(m_matBlock[i])
			delete m_matBlock[i];
	}
	m_matBlock.erase(m_matBlock.begin(),m_matBlock.end());
	m_bAllZero = true;
	*/
}

SymMatrixBlocks& SymMatrixBlocks::operator+=( const SymMatrixBlocks &src ) {
	// check that matrices have compatible sizes; if not, simply return this matrix
	if(m_iSize != src.m_iSize)
		return *this;

	// don't waste time looping through a zero matrix
	if(src.m_bAllZero)
		return *this;

	int i, j, start, length;

	for( i=0; i<src.m_iSize; i++ ) {
		start = src.m_rowstart[i];
		length = src.m_rowstart[i+1] - start;
		for( j = 0; j<length; j++ )
			(*((*this)(i, src.m_col_lookup[start+j] ))) += (*(src.m_matBlock[start+j]));
	}

	return *this;
}

SymMatrixBlocks& SymMatrixBlocks::operator-=(const SymMatrixBlocks &src) {
	// check that matrices have compatible sizes; if not, simply return this matrix
	if(m_iSize != src.m_iSize)
		return *this;

	// don't waste time looping through a zero matrix
	if(src.m_bAllZero)
		return *this;

	int i, j, start, length;

	for( i=0; i<src.m_iSize; i++ ) {
		start = src.m_rowstart[i];
		length = src.m_rowstart[i+1] - start;
		for( j = 0; j<length; j++ )
			(*((*this)(i, src.m_col_lookup[start+j] ))) -= (*(src.m_matBlock[start+j]));
	}

	return *this;
}

SymMatrixBlocks& SymMatrixBlocks::operator*=( double scale ) {
	int i,j;
	int start,length;

	if(!m_bAllZero) {
		for( i=0; i<m_iSize; i++ )
		{
			start = m_rowstart[i];
			length = m_rowstart[i+1] - start;
			for( j = 0; j<length; j++ )
			{
				(*((*this)(i,m_col_lookup[start+j]))) *= scale;
			}
		}
	}

	return *this;
}

void SymMatrixBlocks::SetDiag(DiagonalMatrix &src, DiagonalMatrix &ret) {
	// this function treats DiagonalMatrix src as a column vector with
	// the diagonal entries as the vector components.
	// the result is another diagonal matrix with the diagonal entries
	// treates as the components of the resulting vector from the 
	// matrix-vector multiply.

	ret.Zero();
	// these really must be of the same size to make sense.
	if( m_iSize != src.m_iSize)
		return;

	// if all entries are still zero, we're already done
	if(m_bAllZero) return;
	
	int i, start, length, j , iCol;

	for( i=0; i<m_iSize; i++ )
	{
		start = m_rowstart[i];
		length = m_rowstart[i+1] - start;
		for( j = 0; j<length; j++ )
		{
			// perform multiplication for nonzero matrix blocks
			iCol = m_col_lookup[start+j];
			Mat3d *temp = m_matBlock[start+j];
			ret.m_block[i][0] += temp->n[0] * src.m_block[iCol][0] + temp->n[1] * src.m_block[iCol][1] + temp->n[2] * src.m_block[iCol][2];
			ret.m_block[i][1] += temp->n[3] * src.m_block[iCol][0] + temp->n[4] * src.m_block[iCol][1] + temp->n[5] * src.m_block[iCol][2];
			ret.m_block[i][2] += temp->n[6] * src.m_block[iCol][0] + temp->n[7] * src.m_block[iCol][1] + temp->n[8] * src.m_block[iCol][2];
			// since we stored stuff in upper triangular form, we need to use the transpose of the appropriate matrix block for the complementary matrix-vector mult
			// NOTE: symmetry does not hold for the supermatrix or the submatrices necessarily, but for the matrix obtained by embedding the submatrices into the supermatrix structure
			if( i != iCol )
			{
				ret.m_block[iCol][0] += temp->n[0] * src.m_block[i][0] + temp->n[3] * src.m_block[i][1] + temp->n[6] * src.m_block[i][2];
				ret.m_block[iCol][1] += temp->n[1] * src.m_block[i][0] + temp->n[4] * src.m_block[i][1] + temp->n[7] * src.m_block[i][2];
				ret.m_block[iCol][2] += temp->n[2] * src.m_block[i][0] + temp->n[5] * src.m_block[i][1] + temp->n[8] * src.m_block[i][2];
			}
		}
	}
}

Mat3d* SymMatrixBlocks::operator() (int row, int col) {
	int i, start, l, length;

	m_bAllZero = false;

	// store using upper triangular matrix.
	// so if indexing lower blocks, reverse indexing to get transpose.
	if( row > col  )
	{
		start = row;
		row = col;
		col = start;
	}

	start = m_rowstart[row];
	//if(row+1<m_iSize)
	//	length = m_rowstart[row+1] - start;
	//else
	//	length = 0;
	length = m_rowstart[row+1] - start;
	// if entree already exists, find it.
	for( l=0; l<length; l++ )
		if( m_col_lookup[start+l] == col )
			break;
	// if found a pre-existing entree, return it
	if(  (length != 0) && (l != length) )
	{
		// return the sub-matrix pointer
		return m_matBlock[start+l];
	}
	else
	{
		// find appropriate place in row (keep in order although not necessary)
		for( l=0; l<length; l++ )
			if( m_col_lookup[start+l] >= col )
				break;

		m_col_lookup.insert(m_col_lookup.begin()+start+l, col);
		Mat3d *tempmat = new Mat3d(0,0,0,0,0,0,0,0,0);
		m_matBlock.insert(m_matBlock.begin()+start+l, tempmat);

		for( i=row+1; i<=m_iSize; i++ )
			m_rowstart[i]++;

		return m_matBlock[start+l];
	}
}