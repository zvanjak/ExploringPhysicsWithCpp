#if !defined MML_ODE_SYSTEM_STEPPERS_H
#define MML_ODE_SYSTEM_STEPPERS_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"
#include "interfaces/IODESystemStepCalculator.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "core/ODESystem.h"

#include "algorithms/ODESystemStepCalculators.h"

namespace MML
{
	// Base stepper class
	class StepperBase {
	public:
		// references that stepper gets from the solver
		IODESystem&		_sys;

		Real&					_t;
		Vector<Real>& _x;
		Vector<Real>& _dxdt;

		Real	_absTol, _relTol;
		Real	_eps;

		Real	_tOld;

		Real	_hDone, _hNext;

		Vector<Real> _xout, _xerr;

	public:
		StepperBase(IODESystem& sys, Real& t, Vector<Real>& x, Vector<Real>& dxdt)
			: _sys(sys), _t(t), _x(x), _dxdt(dxdt) {
		}

		virtual void doStep(Real htry, Real eps) = 0;

		void setHDone(Real h) { _hDone = h; }	
		void setHNext(Real h) { _hNext = h; }
	};

	class RK4_CashKarp_Stepper : public StepperBase
	{
	private:
		RK4_CashKarp_Calculator	_stepCalc;

	public:
		RK4_CashKarp_Stepper(IODESystem& sys,	Real& t, Vector<Real>& x, Vector<Real>& dxdt)
			: StepperBase(sys, t, x, dxdt)
		{	}

		void doStep(Real htry, Real eps) override
		{
			/*
			Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and
			adjust stepsize. Input are the dependent variable vector y[1.._numPoints] and its derivative dydx[1.._numPoints]
			at the starting value of the independent variable x. Also input are the stepsize to be attempted
			htry, the required accuracy eps, and the vector xscale[1.._numPoints] against which the error is
			scaled. On output, y and x are replaced by their new values, hdid is the stepsize that was
			actually accomplished, and hnext is the estimated next stepsize. derivs is the user-supplied
			routine that computes the right-hand side derivatives.
			*/
			const Real SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;
			Real errmax, h, htemp, tnew;

			int i, n = _sys.getDim();
			Vector<Real> xerr(n), xscale(n), xtemp(n);

			h = htry;						// Set stepsize to the initial trial value.
			
			// Scaling used to monitor accuracy. This general-purpose choice can be modified if need be
			for (i = 0; i < n; i++)
				xscale[i] = fabs(_x[i]) + fabs(_dxdt[i] * h) + 1e-30;

			for (;;) 
			{
				_stepCalc.calcStep(_sys, _t, _x, _dxdt, h, xtemp, xerr);

				// Evaluating accuracy.
				errmax = 0.0;
				for (int i = 0; i < n; i++)
					errmax = std::max(errmax, fabs(xerr[i] / xscale[i]));
				errmax /= eps;					// Scale relative to required tolerance
				if (errmax <= 1.0)			// Step succeeded. Compute size of next step
					break;

				htemp = SAFETY * h * pow(errmax, PSHRNK);
				// Truncation error too large, reduce stepsize.
				// No more than a factor of 10.
				if (h >= Real{ 0 })
					h = std::max(htemp, 0.1 * h);
				else
					h = std::min(htemp, 0.1 * h);

				tnew = _t + h;

				if (tnew == _t)
					throw("stepsize underflow in rkqs");
			}
			// computing size of the next step
			if (errmax > ERRCON)
				_hNext = SAFETY * h * pow(errmax, PGROW);
			else
				_hNext = 5.0 * h;				// No more than a factor of 5 increase

			_hDone = h;
			_t += h;

			for (i = 0; i < n; i++)
				_x[i] = xtemp[i];
		}
	};
}

#endif