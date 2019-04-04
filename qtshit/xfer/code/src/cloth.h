#ifndef CLOTH_H
#define CLOTH_H

#include "trimesh.h"

class DiagonalMatrix;
class SymMatrixBlocks;
class CForce;

class ClothEdge {
public:
	int v0, v1;
	int f0, f1;
	int opp0, opp1;
	Vec3d n0, n1, e;
	double n0Len, n1Len;
	double eCurLen, eCurAngle;
	double eRestLen, eRestAngle;
	double he0, he1;
	double restBendFactor;

	ClothEdge();
	void update(TriMesh *tm);
	void update(Vec3d &x0, Vec3d &x1, Vec3d &x2, Vec3d &x3);
	void init(TriMesh *tm);
};

class Cloth { 
public:
	TriMesh *tm;
	int numEdges;
	ClothEdge *edges;

	vector<CForce *> forces;

	int numVerts;
	DiagonalMatrix *pos, *vel, *totalForces, *mass;
	DiagonalMatrix *mP, *mPInv, *mz, *mb, *mr, *mc, *mq, *ms, *my, *mdv;
	DiagonalMatrix *mTemp, *dmat;
	SymMatrixBlocks *df_dx, *df_dv;
	SymMatrixBlocks *mA, *mM, *tm1, *tm2;
	Mat3d *mS;

	Cloth(TriMesh *mesh);
	
	void simulate();
};


// Hamilton Chong's force classes:

class CForce  {
public:
	bool		m_bUse;
public:
	CForce(){m_bUse = true;};
	virtual ~CForce(){};
	// Apply force to all vertices
	virtual void Apply( float fTime, DiagonalMatrix &pos, DiagonalMatrix &vel, 
						DiagonalMatrix &totalforce, DiagonalMatrix &mass, 
						SymMatrixBlocks &df_dx, SymMatrixBlocks &df_dv ) {};
};

class CStretch:public CForce {
protected:
	int				m_vert[2];			// vertices to measure stretch between
	double			m_k, m_kDamp;		// physical constants
	double			m_nat;				// natural length
	double			m_max;				// maximum length
public:
	CStretch(int v1, int v2, double c1, double c2, double rest, double max);
	~CStretch(){};
	void Apply( float fTime, DiagonalMatrix &pos, DiagonalMatrix &vel, 
				DiagonalMatrix &totalforce, DiagonalMatrix &mass,
				SymMatrixBlocks &df_dp, SymMatrixBlocks &df_dv );
};

class CAttract:public CForce {
protected:
	int				m_vert;			// vertices to measure stretch between
	double			m_k, m_kDamp;		// physical constants
	double			m_nat;				// natural length
	Vec3d           m_attract;
public:
	CAttract(int v1, Vec3d attractPt, double c1, double c2, double rest);
	~CAttract(){};
	void Apply( float fTime, DiagonalMatrix &pos, DiagonalMatrix &vel, 
				DiagonalMatrix &totalforce, DiagonalMatrix &mass,
				SymMatrixBlocks &df_dp, SymMatrixBlocks &df_dv );
};

class CBend : public CForce {
protected:
	ClothEdge		*edge;
	double			m_k, m_kDamp;		// physical constants
public:
	CBend(ClothEdge *myEdge, double c1, double c2);
	~CBend(){};
	void Apply( float fTime, DiagonalMatrix &pos, DiagonalMatrix &vel, 
				DiagonalMatrix &totalforce, DiagonalMatrix &mass,
				SymMatrixBlocks &df_dp, SymMatrixBlocks &df_dv );
	void Apply( int ind, float fTime, DiagonalMatrix &pos, DiagonalMatrix &vel, 
						Vec3d &tforce, DiagonalMatrix &mass );
};


// Hamilton Chong's sparse matrix classes:

////////////////////////////////////////////////////////
// Diagonal Matrix -- special case of a sparse matrix //
////////////////////////////////////////////////////////
class DiagonalMatrix
{
public:
	int m_iSize;				// number of blocks
	Vec3d*	m_block;		// 3 entries per block

public:
	DiagonalMatrix( int iSize=1 );
	DiagonalMatrix( const DiagonalMatrix &src );
	~DiagonalMatrix();

	// again, for convenience. remember that certain operations may be 
	//   expensive albeit seemingly innocuous to code (not such a big problem
	//   because matrices are diagonal)
	DiagonalMatrix& operator=( const DiagonalMatrix &src );
	void Zero();
	DiagonalMatrix& operator+=( const DiagonalMatrix &src );
	DiagonalMatrix& operator-=( const DiagonalMatrix &src );
	DiagonalMatrix& operator*=( double scale );

	bool Mult( DiagonalMatrix &src , DiagonalMatrix &dst );
	bool well_defined_op(DiagonalMatrix &m1, DiagonalMatrix &m2);

	bool Invert( DiagonalMatrix &ret );
	double DotProductDiagonals( const DiagonalMatrix &src );

};

inline const DiagonalMatrix operator +(const DiagonalMatrix &lhs, const DiagonalMatrix &rhs)
{
	return DiagonalMatrix(lhs) += rhs;
}

inline const DiagonalMatrix operator -(const DiagonalMatrix &lhs, const DiagonalMatrix &rhs)
{
	return DiagonalMatrix(lhs) -= rhs;
}

inline const DiagonalMatrix operator *(const DiagonalMatrix &lhs, double rhs)
{
	return DiagonalMatrix(lhs) *= rhs;
}

inline const DiagonalMatrix operator *(double lhs, const DiagonalMatrix &rhs)
{
	return DiagonalMatrix(rhs) *= lhs;
}

//////////////////////////////////////////////////////////////////////////
// SymMatrixBlocks -- sparse matrix with symmetric 3x3 matrix blocks	//
//////////////////////////////////////////////////////////////////////////
// usually such matrices are composed of sub-blocks of 3x3 matrices.
// sparcity refers to the fact that many of these sub-blocks are all zero.
// however, in many cases, the nonzero sub-blocks are not themselves 
//   sparse (as in the case of cloth simulation).
// Furthermore, the matrix is symmetric, so we only need to store upper triangular
//   supermatrix
// Inspiration for using such an implementation comes from Dean Macri's implementation.
// Dean Macri is Senior Technical Marketing Engineer with Intel’s Developer Relations Division.
class SymMatrixBlocks {
public:
	int							m_iSize;		// size of supermatrix w/ sub-blocks as unit
	int							*m_rowstart;	// beginning of row
	std::vector<int>			m_col_lookup;	// indicies to nonzero blocks (ordered)
	std::vector<Mat3d *>	m_matBlock;		// nonzero blocks of data
	bool						m_bAllZero;		// all entries zero still?
public:
	SymMatrixBlocks(int size);
	SymMatrixBlocks( const SymMatrixBlocks &src );
	~SymMatrixBlocks();

	void Zero(void);							// zero the matrix

	SymMatrixBlocks& operator=( const SymMatrixBlocks &src );
	SymMatrixBlocks& operator+=( const SymMatrixBlocks &src );
	SymMatrixBlocks& operator-=( const SymMatrixBlocks &src );
	SymMatrixBlocks& operator*=( double scale );

	void SetDiag(DiagonalMatrix &src, DiagonalMatrix &ret);

	Mat3d* operator() (int row, int col);	// index overloading -- also allocates space
};

inline const SymMatrixBlocks operator +(const SymMatrixBlocks &lhs, const SymMatrixBlocks &rhs)
{
	return SymMatrixBlocks(lhs) += rhs;
}

inline const SymMatrixBlocks operator -(const SymMatrixBlocks &lhs, const SymMatrixBlocks &rhs)
{
	return SymMatrixBlocks(lhs) -= rhs;
}

inline const SymMatrixBlocks operator *(const SymMatrixBlocks &lhs, double rhs)
{
	return SymMatrixBlocks(lhs) *= rhs;
}

inline const SymMatrixBlocks operator *(double lhs, const SymMatrixBlocks &rhs)
{
	return SymMatrixBlocks(rhs) *= lhs;
}



#endif
