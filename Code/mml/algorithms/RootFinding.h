#if !defined MML_ROOTFINDING_H
#define MML_ROOTFINDING_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "core/Derivation.h"

namespace MML
{
	namespace RootFinding
	{
		// For a given real function 'func' and an initial guessed range x1 to x2, the routine expands
		// the range geometrically until a root is bracketed by the returned values x1 and x2 (in which
		// case function returns true) or until the range becomes unacceptably large (in which case it
		// returns false).
		static bool BracketRoot(const IRealFunction& func, Real& x1, Real& x2, int MaxTry = 50)
		{
			const Real scaleFactor = 1.6;

			if (x1 == x2)
				throw RootFindingError("Bad initial range in BracketRoot");

			Real f1 = func(x1);
			Real f2 = func(x2);

			for (int j = 0; j < MaxTry; j++)
			{
				if (f1 * f2 < 0.0)
					return true;

				if (std::abs(f1) < std::abs(f2))
					f1 = func(x1 += scaleFactor * (x1 - x2));
				else
					f2 = func(x2 += scaleFactor * (x2 - x1));
			}
			return false;
		}

		// For a given a function 'func' defined on the interval[x1, x2], subdivides the interval into
		// numPoints equally spaced segments, and searches for zero crossings of the function.
		// numRoots will be set to the number of bracketing pairs found.
		// If it is positive, the vectors xb1[0..numRoots - 1] and xb2[0..numRoots - 1] will be filled 
		// sequentially with any bracketing pairs that are found.
		static void FindRootBrackets(const IRealFunction& func, const Real x1, const Real x2, const int numPoints, 
																 Vector<Real>& xb1, Vector<Real>& xb2, int& numRoots)
		{
			int numBrackets = 20;
			xb1.Resize(numBrackets);
			xb2.Resize(numBrackets);
			numRoots = 0;
			Real dx = (x2 - x1) / numPoints;
			Real x = x1;
			Real fp = func(x1);

			for (int i = 0; i < numPoints; i++)
			{
				x += dx;
				Real fc = func(x);

				if (fc * fp <= 0.0)
				{
					xb1[numRoots]   = x - dx;
					xb2[numRoots++] = x;
					if (numRoots == numBrackets)
					{
						xb1.Resize(numBrackets * 2, true);
						xb2.Resize(numBrackets * 2, true);
						numBrackets *= 2;
					}
				}
				fp = fc;
			}
		}

		// Using bisection, return the root of a function 'func' known to lie between x1 and x2.
		// The root will be refined until its absolute accuracy is xacc.
		static Real FindRootBisection(const IRealFunction& func, Real x1, Real x2, Real xacc)
		{
			Real dx, xmid, rtb;

			Real f = func(x1);
			Real fmid = func(x2);

			if (f * fmid >= 0.0)
				throw RootFindingError("Root must be bracketed for bisection in FindRootBisection");

			if( f < 0.0 ) {
				dx = x2 - x1;
				rtb = x1;
			}
			else {
				dx = x1 - x2;
				rtb = x2;
			}

			for (int j = 0; j < Defaults::BisectionMaxSteps; j++)
			{
				dx  *= 0.5;
				xmid = rtb + dx;
				fmid = func(xmid);

				if (fmid <= 0.0)
					rtb = xmid;

				if (std::abs(dx) < xacc || fmid == 0.0)
					return rtb;
			}
			throw RootFindingError("Too many bisections in FindRootBisection");
		}

		// Using the Newton-Raphson method, return the root of a function known to lie in the interval
		// [x1, x2]. The root will be refined until its accuracy is known within ˙xacc.
		static Real FindRootNewton(const IRealFunction& func, Real x1, Real x2, Real xacc) 
		{
			Real rtn = 0.5 * (x1 + x2);
			for (int j = 0; j < Defaults::NewtonRaphsonMaxSteps; j++)
			{
				Real f  = func(rtn);
				Real df = Derivation::NDer4(func, rtn);
				
				Real dx = f / df;
				
				rtn -= dx;

				if ((x1 - rtn) * (rtn - x2) < 0.0)
					throw RootFindingError("Jumped out of brackets in FindRootNewton");

				if (std::abs(dx) < xacc)
					return rtn;
			}
			throw RootFindingError("Maximum number of iterations exceeded in FindRootNewton");
		}
	};
}
#endif


