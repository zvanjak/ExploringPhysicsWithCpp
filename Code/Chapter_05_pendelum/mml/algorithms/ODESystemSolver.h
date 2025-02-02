#if !defined MML_ODE_SYSTEM_SOLVERS_H
#define MML_ODE_SYSTEM_SOLVERS_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "core/ODESystem.h"

namespace MML
{
	class RungeKuttaSolverDumb
	{
	public:
		void rk4(Vector<Real>& y, Vector<Real>& dydx, const Real x, const Real h, Vector<Real>& yout, IODESystem& sys)
		{
			int i;
			Real xh, hh, h6;

			int n = (int)y.size();
			Vector<Real> dym(n), dyt(n), yt(n);
			hh = h * 0.5;
			h6 = h / 6.0;
			xh = x + hh;

			for (i = 0; i < n; i++) 
				yt[i] = y[i] + hh * dydx[i];
			
			sys.derivs(xh, yt, dyt);
			
			for (i = 0; i < n; i++) 
				yt[i] = y[i] + hh * dyt[i];
			
			sys.derivs(xh, yt, dym);
			
			for (i = 0; i < n; i++) {
				yt[i] = y[i] + h * dym[i];
				dym[i] += dyt[i];
			}
			
			sys.derivs(x + h, yt, dyt);
			
			for (i = 0; i < n; i++)
				yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
		}

		ODESystemSolutionEqualSpacing integrate(IODESystem& sys, const Vector<Real>& vstart, const Real x1, const Real x2, int numSteps)
		{
			int i, k;
			Real x, h;
			int dim = sys.getDim();

			ODESystemSolutionEqualSpacing sol(dim, numSteps);

			Vector<Real> v(vstart), vout(dim), dv(dim);

			for (i = 0; i < dim; i++) {
				sol.yval[i][0] = v[i];
			}
			sol.xval[0] = x1;
			x = x1;
			h = (x2 - x1) / numSteps;
			
			for (k = 0; k < numSteps; k++) 
			{
				sys.derivs(x, v, dv);

				rk4(v, dv, x, h, vout, sys);

				x += h;

				sol.xval[k + 1] = x;

				for (i = 0; i < dim; i++) {
					v[i] = vout[i];
					sol.yval[i][k + 1] = v[i];
				}
			}

			return sol;
		}
	};

	class RungeKuttaSolverSimple
	{
		void rkck(Vector<Real>& y, Vector<Real>& dydx, const Real x,
			const Real h, Vector<Real>& yout, Vector<Real>& yerr,
			IODESystem& sys)
		{
			static const Real a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875,
				b21 = 0.2, b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9,
				b43 = 1.2, b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0,
				b54 = 35.0 / 27.0, b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0,
				b63 = 575.0 / 13824.0, b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0,
				c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0,
				dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
				dc4 = c4 - 13525.0 / 55296.0, dc5 = -277.00 / 14336.0, dc6 = c6 - 0.25;
			int i;

			int n = y.size();
			Vector<Real> ak2(n), ak3(n), ak4(n), ak5(n), ak6(n), ytemp(n);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + b21 * h * dydx[i];
			sys.derivs(x + a2 * h, ytemp, ak2);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
			sys.derivs(x + a3 * h, ytemp, ak3);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
			sys.derivs(x + a4 * h, ytemp, ak4);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
			sys.derivs(x + a5 * h, ytemp, ak5);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
			sys.derivs(x + a6 * h, ytemp, ak6);
			for (i = 0; i < n; i++)
				yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
			for (i = 0; i < n; i++)
				yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
		}

		void rkqs(Vector<Real>& y, Vector<Real>& dydx, Real& x, const Real htry,
			const Real eps, Vector<Real>& yscal, Real& hdid, Real& hnext,
			IODESystem& sys)
		{
			const Real SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;
			int i;
			Real errmax, h, htemp, xnew;

			int n = (int)y.size();
			h = htry;
			Vector<Real> yerr(n), ytemp(n);

			for (;;) {
				rkck(y, dydx, x, h, ytemp, yerr, sys);

				errmax = 0.0;
				for (i = 0; i < n; i++)
					errmax = std::max(errmax, fabs(yerr[i] / yscal[i]));
				errmax /= eps;
				if (errmax <= 1.0)
					break;

				htemp = SAFETY * h * pow(errmax, PSHRNK);
				h = (h >= Real{ 0 } ? std::max<Real>(htemp, 0.1 * h) : std::min<Real>(htemp, 0.1 * h));
				xnew = x + h;

				if (xnew == x)
					throw("stepsize underflow in rkqs");
			}
			if (errmax > ERRCON)
				hnext = SAFETY * h * pow(errmax, PGROW);
			else
				hnext = 5.0 * h;

			x += (hdid = h);

			for (i = 0; i < n; i++)
				y[i] = ytemp[i];
		}

	public:
		ODESystemSolution integrate(IODESystem& sys, const Vector<Real>& ystart, const Real x1, const Real x2, int maxSteps, Real minSaveInterval,
																const Real eps, const Real h1, const Real hmin, int& nok, int& nbad)
		{
			const int MAXSTP = 10000;
			const Real TINY = 1.0e-30;

			int dim = sys.getDim();
			ODESystemSolution sol(x1, x2, dim, maxSteps);

			int kount = 0;
			int i, nstp;
			Real xsav, x, hnext, hdid, h;
			Vector<Real> yscal(dim), y(ystart), dydx(dim);

			x = x1;
			h = SIGN(h1, x2 - x1);
			nok = nbad = kount = 0;

			if (maxSteps > 0) xsav = x - minSaveInterval * 2.0;
			
			for (nstp = 0; nstp < MAXSTP; nstp++) 
			{
				sys.derivs(x, y, dydx);
			
				for (i = 0; i < dim; i++)
					yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;
				
				if (maxSteps > 0 && kount < maxSteps - 1 && fabs(x - xsav) > fabs(minSaveInterval)) {
					for (i = 0; i < dim; i++) sol._yval[i][kount] = y[i];
					sol._xval[kount++] = x;
					xsav = x;
				}
				
				if ((x + h - x2) * (x + h - x1) > 0.0) 
					h = x2 - x;
				
				rkqs(y, dydx, x, h, eps, yscal, hdid, hnext, sys);
				
				if (hdid == h) 
					++nok; 
				else 
					++nbad;
				
				if ((x - x2) * (x2 - x1) >= 0.0) {
					if (maxSteps != 0) {
						for (i = 0; i < dim; i++) sol._yval[i][kount] = y[i];
						sol._xval[kount++] = x;
					}
					return sol;
				}
				
				if (fabs(hnext) <= hmin) 
					throw("Step size too small in integrate");
				
				h = hnext;
			}
			throw("Too many steps in routine integrate");
		}
	};
}
#endif