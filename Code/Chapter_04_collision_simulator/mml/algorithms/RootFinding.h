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
		// Given a function or functor func and an initial guessed range x1 to x2, the routine expands
		// the range geometrically until a root is bracketed by the returned values x1 and x2(in which
		// case zbrac returns true) or until the range becomes unacceptably large(in which case zbrac
		// returns false).
		static bool BracketRoot(const IRealFunction& func, double& x1, double& x2, int MaxTry = 50)
		{
			const double scaleFactor = 1.6;

			if (x1 == x2)
				throw("Bad initial range in zbrac");

			double f1 = func(x1);
			double f2 = func(x2);

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

		// Given a function or functor fx defined on the interval[x1, x2], subdivide the interval into
		// n equally spaced segments, and search for zero crossings of the function.nroot will be set
		// to the number of bracketing pairs found.If it is positive, the arrays xb1[0..nroot - 1] and
		// xb2[0..nroot - 1] will be filled sequentially with any bracketing pairs that are found.On input,
		// these vectors may have any size, including zero; they will be resized to   nroot.
		static void zbrak(const IRealFunction& fx, const Real x1, const Real x2, const int n, Vector<Real>& xb1,
			Vector<Real>& xb2, int& nroot)
		{
			int nb = 20;
			xb1.resize(nb);
			xb2.resize(nb);
			nroot = 0;
			Real dx = (x2 - x1) / n;
			Real x = x1;
			Real fp = fx(x1);

			for (int i = 0; i < n; i++)
			{
				Real fc = fx(x += dx);

				if (fc * fp <= 0.0)
				{
					xb1[nroot] = x - dx;
					xb2[nroot++] = x;
					if (nroot == nb)
					{
						Vector<Real> tempvec1(xb1), tempvec2(xb2);
						xb1.resize(2 * nb);
						xb2.resize(2 * nb);
						for (int j = 0; j < nb; j++)
						{
							xb1[j] = tempvec1[j];
							xb2[j] = tempvec2[j];
						}
						nb *= 2;
					}
				}
				fp = fc;
			}
		}

		// Using bisection, return the root of a function or functor func known to lie between x1 and x2.
		// The root will be refined until its accuracy is ˙xacc.
		static Real FindRootBisection(const IRealFunction& func, const Real x1, const Real x2, const Real xacc)
		{
			const int JMAX = 50;
			Real dx, xmid, rtb;

			Real f = func(x1);
			Real fmid = func(x2);

			if (f * fmid >= 0.0)
				throw("Root must be bracketed for bisection in FindRootBisection");

			rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
			for (int j = 0; j < JMAX; j++)
			{
				fmid = func(xmid = rtb + (dx *= 0.5));

				if (fmid <= 0.0)
					rtb = xmid;

				if (std::abs(dx) < xacc || fmid == 0.0)
					return rtb;
			}
			throw("Too many bisections in FindRootBisection");
		}

		// Using the false - position method, return the root of a function or functor func known to lie
		// between x1 and x2.The root is refined until its accuracy is ˙xacc
		static Real FindRootFalsePosition(const IRealFunction& func, const Real x1, const Real x2, const Real xacc)
		{
			const int MAXIT = 30;

			Real xl, xh, del;

			Real fl = func(x1);
			Real fh = func(x2);

			if (fl * fh > 0.0)
				throw("Root must be bracketed in FindRootFalsePosition");

			if (fl < 0.0) {
				xl = x1;
				xh = x2;
			}
			else {
				xl = x2;
				xh = x1;
				std::swap(fl, fh);
			}
			Real dx = xh - xl;
			for (int j = 0; j < MAXIT; j++)
			{
				Real rtf = xl + dx * fl / (fl - fh);
				Real f = func(rtf);

				if (f < 0.0) {
					del = xl - rtf;
					xl = rtf;
					fl = f;
				}
				else {
					del = xh - rtf;
					xh = rtf;
					fh = f;
				}
				dx = xh - xl;
				if (std::abs(del) < xacc || f == 0.0)
					return rtf;
			}
			throw("Maximum number of iterations exceeded in FindRootFalsePosition");
		}

		// Using the secant method, return the root of a function or functor func thought to lie between
		// x1 and x2.The root is refined until its accuracy is ˙xacc.
		static Real FindRootSecant(const IRealFunction& func, const Real x1, const Real x2, const Real xacc) {
			const int MAXIT = 30;
			Real xl, rts;
			Real fl = func(x1);
			Real f = func(x2);
			if (std::abs(fl) < std::abs(f)) {
				rts = x1;
				xl = x2;
				std::swap(fl, f);
			}
			else {
				xl = x1;
				rts = x2;
			}
			for (int j = 0; j < MAXIT; j++)
			{
				Real dx = (xl - rts) * f / (f - fl);
				xl = rts;
				fl = f;
				rts += dx;
				f = func(rts);

				if (std::abs(dx) < xacc || f == 0.0)
					return rts;
			}
			throw("Maximum number of iterations exceeded in FindRootSecant");
		}

		// Using Ridders’ method, return the root of a function or functor func known to lie between x1
		// and x2.The root will be refined to an approximate accuracy xacc.
		static Real FindRootRidders(const IRealFunction& func, const Real x1, const Real x2, const Real xacc) {
			const int MAXIT = 60;
			Real fl = func(x1);
			Real fh = func(x2);
			if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0))
			{
				Real xl = x1;
				Real xh = x2;
				Real ans = -9.99e99;

				for (int j = 0; j < MAXIT; j++)
				{
					Real xm = 0.5 * (xl + xh);
					Real fm = func(xm);
					Real s = sqrt(fm * fm - fl * fh);

					if (s == 0.0)
						return ans;

					Real xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);

					if (std::abs(xnew - ans) <= xacc)
						return ans;

					ans = xnew;
					Real fnew = func(ans);

					if (fnew == 0.0)
						return ans;

					if (SIGN(fm, fnew) != fm) {
						xl = xm;
						fl = fm;
						xh = ans;
						fh = fnew;
					}
					else if (SIGN(fl, fnew) != fl) {
						xh = ans;
						fh = fnew;
					}
					else if (SIGN(fh, fnew) != fh) {
						xl = ans;
						fl = fnew;
					}
					else throw("never get here.");

					if (std::abs(xh - xl) <= xacc)
						return ans;
				}
				throw("FindRootRidders exceed maximum iterations");
			}
			else {
				if (fl == 0.0)
					return x1;
				if (fh == 0.0)
					return x2;

				throw("root must be bracketed in FindRootRidders.");
			}
		}

		// Using Brent’s method, return the root of a function or functor func known to lie between x1
		// and x2.The root will be refined until its accuracy is tol.
		static Real FindRootBrent(IRealFunction& func, const Real x1, const Real x2, const Real tol)
		{
			const int ITMAX = 100;
			const Real EPS = std::numeric_limits<Real>::epsilon();
			Real a = x1, b = x2, c = x2, d, e, fa = func(a), fb = func(b), fc, p, q, r, s, tol1, xm;

			if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
				throw("Root must be bracketed in FindRootBrent");

			fc = fb;
			for (int iter = 0; iter < ITMAX; iter++)
			{
				if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
					c = a;
					fc = fa;
					e = d = b - a;
				}
				if (std::abs(fc) < std::abs(fb)) {
					a = b;
					b = c;
					c = a;
					fa = fb;
					fb = fc;
					fc = fa;
				}

				tol1 = 2.0 * EPS * std::abs(b) + 0.5 * tol;
				xm = 0.5 * (c - b);

				if (std::abs(xm) <= tol1 || fb == 0.0)
					return b;

				if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
					s = fb / fa;
					if (a == c) {
						p = 2.0 * xm * s;
						q = 1.0 - s;
					}
					else {
						q = fa / fc;
						r = fb / fc;
						p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
						q = (q - 1.0) * (r - 1.0) * (s - 1.0);
					}

					if (p > 0.0)
						q = -q;
					p = std::abs(p);

					Real min1 = 3.0 * xm * q - std::abs(tol1 * q);
					Real min2 = std::abs(e * q);

					if (2.0 * p < (min1 < min2 ? min1 : min2)) {
						e = d;
						d = p / q;
					}
					else {
						d = xm;
						e = d;
					}
				}
				else {
					d = xm;
					e = d;
				}
				a = b;
				fa = fb;
				if (std::abs(d) > tol1)
					b += d;
				else
					b += SIGN(tol1, xm);
				fb = func(b);
			}
			throw("Maximum number of iterations exceeded in FindRootBrent");
		}

		// Using the Newton-Raphson method, return the root of a function known to lie in the interval
		// x1; x2.The root will be refined until its accuracy is known within ˙xacc.
		static Real FindRootNewton(const IRealFunction& funcd, const Real x1, const Real x2, const Real xacc) {
			const int JMAX = 20;
			Real rtn = 0.5 * (x1 + x2);
			for (int j = 0; j < JMAX; j++)
			{
				Real f = funcd(rtn);
				Real df = Derivation::NDer4(funcd, rtn);
				Real dx = f / df;
				rtn -= dx;
				if ((x1 - rtn) * (rtn - x2) < 0.0)
					throw("Jumped out of brackets in rtnewt");
				if (std::abs(dx) < xacc)
					return rtn;
			}
			throw("Maximum number of iterations exceeded in rtnewt");
		}

		//Using a combination of Newton - Raphson and bisection, return the root of a function bracketed
		//between x1 and x2.The root will be refined until its accuracy is known within ˙xacc.funcd
		//is a user - supplied struct that returns the function value as a functor and the first derivative of
		//the function at the point x as the function df(see text).
		static Real FindRootSafe(const IRealFunction& funcd, const Real x1, const Real x2, const Real xacc) {
			const int MAXIT = 100;
			Real xh, xl;
			Real fl = funcd(x1);
			Real fh = funcd(x2);

			if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
				throw("Root must be bracketed in rtsafe");

			if (fl == 0.0) return x1;
			if (fh == 0.0) return x2;

			if (fl < 0.0) {
				xl = x1;
				xh = x2;
			}
			else {
				xh = x1;
				xl = x2;
			}
			Real rts = 0.5 * (x1 + x2);
			Real dxold = std::abs(x2 - x1);
			Real dx = dxold;

			Real f = funcd(rts);
			Real df = Derivation::NDer4(funcd, rts);

			for (int j = 0; j < MAXIT; j++)
			{
				if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0)
					|| (std::abs(2.0 * f) > std::abs(dxold * df))) {
					dxold = dx;
					dx = 0.5 * (xh - xl);
					rts = xl + dx;
					if (xl == rts) return rts;
				}
				else {
					dxold = dx;
					dx = f / df;
					Real temp = rts;
					rts -= dx;
					if (temp == rts) return rts;
				}
				if (std::abs(dx) < xacc)
					return rts;

				f = funcd(rts);
				df = Derivation::NDer4(funcd, rts);

				if (f < 0.0)
					xl = rts;
				else
					xh = rts;
			}
			throw("Maximum number of iterations exceeded in rtsafe");
		}
	};
}
#endif


