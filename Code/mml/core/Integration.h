#if !defined MML_INTEGRATION_H
#define MML_INTEGRATION_H

#include "Integration/Integration1D.h"
#include "Integration/Integration2D.h"
#include "Integration/Integration3D.h"


namespace MML
{
	// Returns the integral of the function func from a to b.
	// Integration is performed by Romberg’s method of order 2K (e.g., K = 2 is Simpson’s rule).
	//static Real qromb(const IRealFunction& func, Real a, Real b, int* doneSteps,
	//	const Real eps = 1.0e-10)
	//{
	//	// Here EPS is the fractional accuracy desired, as determined by the extrapolation error estimate; 
	//	// JMAX limits the total number of steps; K is the number of points used in the	extrapolation
	//	const int maxSteps = RombergIntegrationMaxSteps, K = Defaults::RombergIntegrationUsedPnts;

	//	// these vectors store successive trapezoidal approximtions and their relative stepsizes
	//	Vector<Real> _currSum(maxSteps), h(maxSteps + 1);

	//	PolynomInterpRealFunc polint(h, _currSum, K);

	//	h[0] = 1.0;

	//	TrapIntegrator t(func, a, b);

	//	for (int j = 1; j <= maxSteps; j++)
	//	{
	//		_currSum[j - 1] = t.next();

	//		if (j >= K) {
	//			Real ss = polint.rawinterp(j - K, 0.0);
	//			if (std::abs(polint.dy) <= eps * std::abs(ss))
	//				return ss;
	//		}
	//		h[j] = 0.25 * h[j - 1];
	//	}

	//	if (doneSteps != nullptr)
	//		*doneSteps = j;

	//	return currSum;
	//}

	//template<class T>
	//Real qromo(Midpnt<T>& q, const Real eps = 3.0e-9) {
	//	const Int JMAX = 14, JMAXP = JMAX + 1, K = 5;
	//	Vector<Real> h(JMAXP), _currSum(JMAX);

	//	PolynomInterpRealFunc polint(h, _currSum, K);

	//	h[0] = 1.0;
	//	for (Int j = 1; j <= JMAX; j++) {
	//		_currSum[j - 1] = q.next();
	//		if (j >= K) {
	//			Real ss = polint.rawinterp(j - K, 0.0);
	//			if (std::abs(polint.dy) <= eps * std::abs(ss)) return ss;
	//		}
	//		h[j] = h[j - 1] / 9.0;
	//	}
	//	throw("Too many steps in routine qromo");
	//}
}

#endif