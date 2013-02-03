/**@brief      Useful functions for operating on triads in C.
 * @author     Rob Simpson
 * @contact    r.simpson11@imperial.ac.uk
 * @version    0.0
 * @date       31/01/2013
 * @pre        None
 * @warning    None
 */

#include <cmath>
#include <stdio.h>
#include <cstdlib>

void AddTriad(const double* x1, const double* x2, double* xOut) {
	xOut[0] = x1[0] + x2[0];
	xOut[1] = x1[1] + x2[1];
	xOut[2] = x1[2] + x2[2];
}

void SubTriad(const double* x1, const double* x2, double* xOut) {
	xOut[0] = x1[0] - x2[0];
	xOut[1] = x1[1] - x2[1];
	xOut[2] = x1[2] - x2[2];
}

void MulTriad(const double* x1, const double Factor, double* xOut) {
	xOut[0] = Factor*x1[0];
	xOut[1] = Factor*x1[1];
	xOut[2] = Factor*x1[2];
}

void DivTriad(const double* x1, const double Factor, double* xOut) {
	if (Factor == 0.0) {
		fprintf(stderr, "Division by zero! Aborting...\n");
		exit(EXIT_FAILURE);
	}
	xOut[0] = x1[0]/Factor;
	xOut[1] = x1[1]/Factor;
	xOut[2] = x1[2]/Factor;
}

double DotTriad(const double* x1, const double* x2) {
	return (x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2]);
}

double NormTriad(const double* x1) {
	return (sqrt(x1[0] * x1[0] + x1[1] * x1[1] + x1[2] * x1[2]));
}

void NormaliseTriad(const double* x1, double* xOut) {
	double Factor = NormTriad(x1);
	DivTriad(x1,Factor,xOut);
}

void CrossTriad(const double* x1, const double* x2, double* xOut) {
	xOut[0] = x1[1]*x2[2] - x1[2]*x2[1];
	xOut[1] = x1[2]*x2[0] - x1[0]*x2[2];
	xOut[2] = x1[0]*x2[1] - x1[1]*x2[0];
}

void BilinearMapTriad(const double* p1, const double* p2, \
					  const double* p3, const double* p4, \
					  double* pOut) {
	/** @brief Bilinear map on aero surface.
	 * @details Maps to centre point only just now.
	 */
	pOut[0] = (p1[0] +p2[0] +p3[0] + p4[0])/4.0;
	pOut[1] = (p1[1] +p2[1] +p3[1] + p4[1])/4.0;
	pOut[2] = (p1[2] +p2[2] +p3[2] + p4[2])/4.0;
}

void CopyTriad(double* pTarget, double* pSrc) {
	pTarget[0] = pSrc[0];
	pTarget[1] = pSrc[1];
	pTarget[2] = pSrc[2];
}

void PrintTriad(double* x) {
	printf("(%f %f %f)",x[0],x[1],x[2]);
}
