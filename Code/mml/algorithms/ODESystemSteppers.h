#if !defined MML_ODE_SYSTEM_STEPPERS_H
#define MML_ODE_SYSTEM_STEPPERS_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"
#include "interfaces/IODESystemStepCalculator.h"

#include "core/ODESystem.h"

#include "algorithms/ODESystemStepCalculators.h"

namespace MML
{
	// Base stepper class
	class StepperBase {
	protected:
		// references that stepper gets from the solver
		const IODESystem&		_sys;

		Real&					_t;
		Vector<Real>& _x;
		Vector<Real>& _dxdt;

		Real	_absTol, _relTol;
		Real	_eps;

		Real	_tOld;
		Real	_hDone, _hNext;

		Vector<Real> _xout, _xerr;

	public:
		StepperBase(const IODESystem& sys, Real& t, Vector<Real>& x, Vector<Real>& dxdt)
			: _sys(sys), _t(t), _x(x), _dxdt(dxdt) {
		}

		virtual void doStep(Real htry, Real eps) = 0;

		Real	hDone() const { return _hDone; }
		Real& hDone()				{ return _hDone; }

		Real	hNext() const { return _hNext; }
		Real& hNext()				{ return _hNext; }
	};

	class RK5_CashKarp_Stepper : public StepperBase
	{
	private:
		RK5_CashKarp_Calculator	_stepCalc;

	public:
		// Initializing stepper with references to the solver's system, time, state and derivative vectors.
		RK5_CashKarp_Stepper(const IODESystem& sys,	Real& t, Vector<Real>& x, Vector<Real>& dxdt)
			: StepperBase(sys, t, x, dxdt) { }

		/*
		Implements fifth-order Runge-Kutta step with monitoring of local truncation error 
		to ensure accuracy and	adjust stepsize. 
		Inputs are the stepsize to be attempted htry and the required accuracy eps.
		On output, references to x and t are replaced by their new values, hdid is the stepsize that was
		actually accomplished, and hnext is the estimated next stepsize. 
		*/
		void doStep(Real htry, Real eps) override
		{
			const Real SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;
			Real	errmax, h, htemp, tnew;
			int		i, n = _sys.getDim();
			Vector<Real> xerr(n), xscale(n), xtemp(n);

			h = htry;						// Set stepsize to the initial trial value.
			
			// Scaling used to monitor accuracy. 
			for (i = 0; i < n; i++)
				xscale[i] = std::abs(_x[i]) + std::abs(_dxdt[i] * h) + 1e-30;

			for (;;) {
				_stepCalc.calcStep(_sys, _t, _x, _dxdt, h, xtemp, xerr);

				errmax = 0.0;								// Evaluating accuracy.
				for (int i = 0; i < n; i++)
					errmax = std::max(errmax, fabs(xerr[i] / xscale[i]));
				errmax /= eps;							// Scale relative to required tolerance
				if (errmax <= 1.0)			
					break;										// Step succeeded. Compute size of next step.	

				htemp = SAFETY * h * pow(errmax, PSHRNK);
				// Truncation error too large, reduce stepsize, but no more than a factor of 10.
				if (h >= Real{ 0 })
					h = std::max(htemp, 0.1 * h);
				else
					h = std::min(htemp, 0.1 * h);

				tnew = _t + h;
				if (tnew == _t)	
					throw ODESolverError("Stepsize underflow in RK5_CashKarp_Stepper::doStep()");
			}
			// computing size of the next step
			if (errmax > ERRCON)
				_hNext = SAFETY * h * pow(errmax, PGROW);
			else
				_hNext = 5.0 * h;				// No more than a factor of 5 increase

			_hDone = h;
			_t += h;
			_x = xtemp;
		}
	};
}

#endif