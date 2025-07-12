#if !defined MML_ODE_SYSTEM_SOLVERS_H
#define MML_ODE_SYSTEM_SOLVERS_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"
#include "interfaces/IODESystemStepCalculator.h"
#include "interfaces/IODESystemStepper.h"

#include "base/ODESystem.h"
#include "base/ODESystemSolution.h"


namespace MML
{
	class ODESystemFixedStepSolver
	{
		const IODESystem& _odeSys;
		const IODESystemStepCalculator& _stepCalc;

	public:
		ODESystemFixedStepSolver(const IODESystem& inOdeSys, const IODESystemStepCalculator& inStepCalc)
			: _odeSys(inOdeSys), _stepCalc(inStepCalc) {
		}

		// For given numSteps, ODESystemSolution will have numSteps+1 values!
		ODESystemSolution integrate(const Vector<Real>& initCond, Real t1, Real t2, int numSteps)
		{
			int		dim = _odeSys.getDim();

			ODESystemSolution sol(t1, t2, dim, numSteps);
			Vector<Real>			x(initCond), x_out(dim), dxdt(dim), x_err(dim);

			sol.fillValues(0, t1, x);		// store initial values

			Real t = t1;
			Real h = (t2 - t1) / numSteps;

			for (int k = 1; k <= numSteps; k++) {
				_odeSys.derivs(t, x, dxdt);

				_stepCalc.calcStep(_odeSys, t, x, dxdt, h, x_out, x_err);

				t += h;

				x = x_out;
				sol.fillValues(k, t, x);
			}

			return sol;
		}
	};

	// CHECK!!!
	class ODESystemLeapfrogSolver
	{
		IODESystem& _odeSys;

	public:
		ODESystemLeapfrogSolver(IODESystem& inOdeSys)	: _odeSys(inOdeSys) {
		}
		// Leapfrog integrator for ODE systems.
		ODESystemSolution integrate(const Vector<Real>& initCond, Real t1, Real t2, int numSteps)
		{
			int dim = _odeSys.getDim();
			int half_n = dim / 2;
			ODESystemSolution sol(t1, t2, dim, numSteps);

			Vector<Real> x(initCond); // current state
			Vector<Real> dxdt(dim);   // derivatives
			Vector<Real> x_out(dim);  // next state

			sol.fillValues(0, t1, x); // store initial values

			Real t = t1;
			Real h = (t2 - t1) / numSteps;

			// Initial acceleration
			_odeSys.derivs(t, x, dxdt);

			// Split x into positions and velocities
			Vector<Real> pos(half_n), vel(half_n), acc(half_n);
			for (int i = 0; i < half_n; ++i) {
				pos[i] = x[i];
				vel[i] = x[half_n + i];
				acc[i] = dxdt[half_n + i];
			}

			for (int k = 1; k <= numSteps; ++k) {
				// 1. Advance velocities by half step
				for (int i = 0; i < half_n; ++i)
					vel[i] += 0.5 * h * acc[i];

				// 2. Advance positions by full step
				for (int i = 0; i < half_n; ++i)
					pos[i] += h * vel[i];

				t += h;

				// 3. Compute new acceleration at new position
				for (int i = 0; i < half_n; ++i) {
					x[i] = pos[i];
					x[half_n + i] = vel[i]; // temporary, needed for acceleration calculation
				}
				_odeSys.derivs(t, x, dxdt);
				for (int i = 0; i < half_n; ++i)
					acc[i] = dxdt[half_n + i];

				// 4. Advance velocities by another half step
				for (int i = 0; i < half_n; ++i)
					vel[i] += 0.5 * h * acc[i];

				// 5. Store new state in x_out and save
				for (int i = 0; i < half_n; ++i) {
					x_out[i] = pos[i];
					x_out[half_n + i] = vel[i];
				}
				sol.fillValues(k, t, x_out);

				// Prepare for next step
				x = x_out;
			}

			return sol;
		}
	};

	template<class Stepper>	class ODESystemSolver
	{
		IODESystem& _sys;
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
		// Integrates starting values initCond[1..nvar] from t1 to t2 with accuracy eps 
		ODESystemSolution integrate(const Vector<Real>& initCond,
																Real t1, Real t2, Real minSaveInterval,
																Real eps, Real h1, Real hmin = 0)
		{
			Real xsav, h;
			int	i, stepNum, numSavedSteps = 0, sysDim = _sys.getDim();
			int expectedSteps = (int)((t2 - t1) / minSaveInterval) + 1;
			ODESystemSolution sol(t1, t2, sysDim, expectedSteps);

			_curr_t = t1;
			_curr_x = initCond;

			// making sure our steps are pointing the right way
			if (t2 - t1 > 0.0) h = std::abs(h1);
			else h = -std::abs(h1);

			if (expectedSteps > 0)	xsav = _curr_t - minSaveInterval * 2.0;

			for (stepNum = 0; stepNum < Defaults::ODESolverMaxSteps; stepNum++) {
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

				if (_stepper.hDone() == h)
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
				if (fabs(_stepper.hNext()) <= hmin) throw ODESolverError("Step size too small in integrate");

				h = _stepper.hNext();
			}
			throw ODESolverError("Too many steps in routine integrate");
		}
	};
}
#endif