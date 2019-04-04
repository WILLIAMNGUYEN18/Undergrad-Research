#include "expmap.h"

// eponential map calculations
// based on Grassia's C code:
//   http://www.acm.org/jgt/papers/Grassia98/code.html

/* crossover point to Taylor Series approximation.  Figuring 16
 * decimal digits of precision for doubles, the Taylor approximation
 * should be indistinguishable (to machine precision) for angles
 * of magnitude less than 1e-4. To be conservative, we add on three
 * additional orders of magnitude.  */
const double MIN_ANGLE = 1e-7;

/* Angle beyond which we perform dynamic reparameterization of a 3 DOF EM */
const double CUTOFF_ANGLE = M_PI;


int normalizeExpMap(Vec3d &em) {
	int rep = 0;
	double theta = em.length();
	if (theta > CUTOFF_ANGLE) {
		double scl = theta;
		if (theta > 2*M_PI) {	/* first get theta into range 0..2PI */
			theta = fmod(theta, 2*M_PI);
			scl = theta/scl;
			em *= sc1;
			rep = 1;
		}
		if (theta > CUTOFF_ANGLE) {
			scl = theta;
			theta = 2*M_PI - theta;
			scl = 1.0 - 2*M_PI/scl;
			em *= sc1;
			rep = 1;
		}
	}

	return rep;	
}

QuatNorm expMapToQuatNorm(Vec3d em, bool normalize) {
	double cosp, sinp, theta;
	QuatNorm q;

	if (normalize)
		normalizeExpMap(em);
		
	theta = em.length();
	cosp = cos(.5*theta);
	sinp = sin(.5*theta);

	q.w = cosp;
	if (theta < MIN_ANGLE)
		em *= 0.5 - theta*theta/48.0;	/* Taylor Series for sinc */
	else
		em *= sinp/theta;
	q.x = em[0];
	q.y = em[1];
	q.z = em[2];
		
	return q;
}

Mat4d expMapToMat(Vec3d em) {
	QuatNorm q = expMapToQuatNorm(ex);
	return q.toMatrixD();
}

/* -----------------------------------------------------------------
 * 'Partial_Q_Partial_3V' Partial derivative of quaternion wrt i'th
 * component of EM vector 'v'
 * -----------------------------------------------------------------*/
void Partial_Q_Partial_3V(Vec3d em[3], QuatNorm *dqdx) {
	double   theta = em.length();
	double   cosp = cos(.5*theta), sinp = sin(.5*theta);
    
	int i;
	for (i=0; i < 3; i++) {
		/* This is an efficient implementation of the derivatives given
		* in Appendix A of the paper with common subexpressions factored out */
		if (theta < MIN_ANGLE) {
			const int i2 = (i+1)%3, i3 = (i+2)%3;
			double Tsinc = 0.5 - theta*theta/48.0;
			double vTerm = v[i] * (theta*theta/40.0 - 1.0) / 24.0;
			
			dqdx[W] = -.5*v[i]*Tsinc;
			dqdx[i]  = v[i]* vTerm + Tsinc;
			dqdx[i2] = v[i2]*vTerm;
			dqdx[i3] = v[i3]*vTerm;
		}
		else {
			const int i2 = (i+1)%3, i3 = (i+2)%3;
			const double  ang = 1.0/theta, ang2 = ang*ang*v[i], sang = sinp*ang;
			const double  cterm = ang2*(.5*cosp - sang);
			
			dqdx[i]  = cterm*v[i] + sang;
			dqdx[i2] = cterm*v[i2];
			dqdx[i3] = cterm*v[i3];
			dqdx[W] = -.5*v[i]*sang;
		}
	}
}