#if !defined MML_ODE_SYSTEM_SOLVERS_H
#define MML_ODE_SYSTEM_SOLVERS_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"
#include "interfaces/IODESystemStepCalculator.h"
#include "interfaces/IODESystemStepper.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "core/ODESystem.h"
#include "core/ODESystemSolution.h"


namespace MML
{
	class ODESystemFixedStepSolver
	{
		IODESystem& _odeSys;
		IODESystemStepCalculator& _stepCalc;

	public:
		ODESystemFixedStepSolver(IODESystem& inOdeSys, IODESystemStepCalculator& inStepCalc)
			: _odeSys(inOdeSys), _stepCalc(inStepCalc) { }
		
		// For given numSteps, ODESystemSolution will have numSteps+1 values!
		ODESystemSolution integrate(const Vector<Real>& initCond, Real x1, Real x2, int numSteps)
		{
			int		dim = _odeSys.getDim();

			ODESystemSolution sol(x1, x2, dim, numSteps);
			Vector<Real>			y(initCond), y_out(dim), dydx(dim), y_err(dim);

			sol.fillValues(0, x1, y);

			Real x = x1;
			Real h = (x2 - x1) / numSteps;

			for (int k = 1; k <= numSteps; k++) {
				_odeSys.derivs(x, y, dydx);

				_stepCalc.calcStep(_odeSys, x, y, dydx, h, y_out, y_err);

				x += h;

				y = y_out;
				sol.fillValues(k, x, y);
			}

			return sol;
		}
	};

	template<class Stepper>	class ODESystemSolver 
	{
		IODESystem&		_sys;
		Stepper				_stepper;

		// used as reference by stepper!
		Real					_curr_t;               
		Vector<Real>	_curr_x;
		Vector<Real>	_curr_dxdt;

	public:
		ODESystemSolver(IODESystem& sys)
			: _sys(sys), _stepper(sys, _curr_t, _curr_x, _curr_dxdt) 
		{
			_curr_x.Resize(_sys.getDim());
			_curr_dxdt.Resize(_sys.getDim());
		}

		int getDim() { return _sys.getDim(); }

		// ODE system integrator using given stepper for adaptive stepsize control. 
		// Integrates starting values ystart[1..nvar] from t1 to t2 with accuracy eps 
		ODESystemSolution integrate(const Vector<Real>& initCond,
																Real t1, Real t2, Real minSaveInterval,
																Real eps, Real h1, Real hmin = 0)
		{
			Real xsav, h;
			int i, stepNum, numSavedSteps=0, sysDim = _sys.getDim();
			int expectedSteps = (int)((t2 - t1) / h1) + 1;
			ODESystemSolution sol(t1, t2, sysDim, expectedSteps);

			_curr_t = t1;
			_curr_x = initCond;
			h = SIGN(h1, t2 - t1);

			if (expectedSteps > 0)	xsav = _curr_t - minSaveInterval * 2.0;

			for (stepNum = 0; stepNum < Defaults::ODESolverMaxSteps; stepNum++)	{
				_sys.derivs(_curr_t, _curr_x, _curr_dxdt);

				// storing intermediate results, if we have advanced enough
				if (expectedSteps > 0 && fabs(_curr_t - xsav) > fabs(minSaveInterval))
				{
					sol.fillValues(numSavedSteps, _curr_t, _curr_x);
					numSavedSteps++;
					xsav = _curr_t;
				}
				// If stepsize overshoots end of interval, decrease to adjust
				if ((_curr_t + h - t2) * (_curr_t + h - t1) > 0.0)
					h = t2 - _curr_t;

				_stepper.doStep(h, eps);

				if (_stepper._hDone == h) 
					sol.incrementNumStepsOK();
				else 
					sol.incrementNumStepsBad();

				// check if we are done
				if ((_curr_t - t2) * (t2 - t1) >= 0.0) {
					if (expectedSteps != 0) {
						sol.fillValues(numSavedSteps, _curr_t, _curr_x);
					}
					sol.setFinalSize(numSavedSteps);
					return sol;
				}
				if (fabs(_stepper._hNext) <= hmin) throw ODESolverError("Step size too small in integrate");

				h = _stepper._hNext;
			}
			throw ODESolverError("Too many steps in routine integrate");
		}
	};
}
#endif