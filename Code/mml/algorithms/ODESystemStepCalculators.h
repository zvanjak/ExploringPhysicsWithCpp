#if !defined MML_ODE_SYSTEM_STEP_CALCULATORS_H
#define MML_ODE_SYSTEM_STEP_CALCULATORS_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"
#include "interfaces/IODESystemStepCalculator.h"

#include "core/ODESystem.h"

namespace MML
{
	class EulerStep_Calculator : public IODESystemStepCalculator
	{
	public:
		// For a given IODESystem of dimension n, and given initial values for the variables x_start[0..n-1] 
		// and their derivatives dxdt[0..n-1] known at t, uses the Euler method to advance the solution 
		// over an interval h and return the incremented variables as x_out[0..n-1].
		void calcStep(const IODESystem& odeSystem,
									const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
									const Real h, Vector<Real>& x_out, Vector<Real>& x_err_out) const override
		{
			int i, n = odeSystem.getDim();
			for (i = 0; i < n; i++)
				x_out[i] = x_start[i] + h * dxdt[i];
		}
	};

	class MidpointStep_Calculator : public IODESystemStepCalculator
	{
	public:
		// For a given ODESystem, of dimension n, and given initial values for the variables x_start[0..n-1] 
		// and their derivatives dxdt[0..n-1] known at t, uses the midpoint method to advance the solution 
		// over an interval h and return the incremented variables as xout[0..n-1].
		void calcStep(const IODESystem& odeSystem,
									const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
									const Real h, Vector<Real>& x_out, Vector<Real>& x_err_out) const override
		{
			int i, n = odeSystem.getDim();
			Vector<Real> x_mid(n);

			for (i = 0; i < n; i++)
				x_mid[i] = x_start[i] + 0.5 * h * dxdt[i];
			
			odeSystem.derivs(t + 0.5 * h, x_mid, x_out);
		}
	};

	class RK4_FixedStep_Calculator : public IODESystemStepCalculator
	{
	public:
		// For a given ODESystem, of dimension n, and given initial values for the variables x_start[0..n-1] 
		// and their derivatives dxdt[0..n-1] known at t, uses the fourth-order Runge-Kutta method 
		// to advance the solution over an interval h and return the incremented variables as xout[0..n-1].
		void calcStep(const IODESystem& odeSystem, 
									const Real t, const Vector<Real> &x_start, const Vector<Real> &dxdt,
									const Real h, Vector<Real> &x_out, Vector<Real> &x_err_out) const override
		{
			int i, n = odeSystem.getDim();
			Vector<Real> dx_mid(n), dx_temp(n), x_temp(n);

			Real xh, hh, h6;
			hh = h * 0.5;
			h6 = h / 6.0;
			xh = t + hh;

			for (i = 0; i < n; i++)												// First step
				x_temp[i] = x_start[i] + hh * dxdt[i];

			odeSystem.derivs(xh, x_temp, dx_temp);				// Second step

			for (i = 0; i < n; i++)
				x_temp[i] = x_start[i] + hh * dx_temp[i];

			odeSystem.derivs(xh, x_temp, dx_mid);					// Third step

			for (i = 0; i < n; i++) {
				x_temp[i] = x_start[i] + h * dx_mid[i];
				dx_mid[i] += dx_temp[i];
			}

			odeSystem.derivs(t + h, x_temp, dx_temp);			// Fourth step	

			for (i = 0; i < n; i++)
				x_out[i] = x_start[i] + h6 * (dxdt[i] + dx_temp[i] + 2.0 * dx_mid[i]);
		}
	};

	class RK5_CashKarp_Calculator : public IODESystemStepCalculator
	{
	public:
		/*
		Given initial values for n variables x[0..n-1] and their derivatives dxdt[0..n-1] known at t, 
		usees	the fifth-order Cash-Karp Runge-Kutta method to advance the solution over an interval h
		and return the incremented variables as x_out[0..n-1].
		Also returns an estimate of the local truncation error in x_out using the embedded fourth-order method.
		*/
		void calcStep(const IODESystem& sys, 
									Real t, const Vector<Real>& x, const Vector<Real>& dxdt, 
									Real h, Vector<Real>& x_out, Vector<Real>& x_err) const override
		{
			static const Real a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875,
				b21 = 0.2, b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9,
				b43 = 1.2, b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0,
				b54 = 35.0 / 27.0, b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0,
				b63 = 575.0 / 13824.0, b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0,
				c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0,
				dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
				dc4 = c4 - 13525.0 / 55296.0, dc5 = -277.00 / 14336.0, dc6 = c6 - 0.25;

			int i, n = x.size();
			Vector<Real> ak2(n), ak3(n), ak4(n), ak5(n), ak6(n), xtemp(n);

			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + b21 * h * dxdt[i];

			sys.derivs(t + a2 * h, xtemp, ak2);
			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + h * (b31 * dxdt[i] + b32 * ak2[i]);

			sys.derivs(t + a3 * h, xtemp, ak3);
			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + h * (b41 * dxdt[i] + b42 * ak2[i] + b43 * ak3[i]);

			sys.derivs(t + a4 * h, xtemp, ak4);
			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + h * (b51 * dxdt[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);

			sys.derivs(t + a5 * h, xtemp, ak5);
			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + h * (b61 * dxdt[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);

			sys.derivs(t + a6 * h, xtemp, ak6);

			for (i = 0; i < n; i++)
				x_out[i] = x[i] + h * (c1 * dxdt[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
			for (i = 0; i < n; i++)
				x_err[i] = h * (dc1 * dxdt[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
		}
	};

	class StepCalculators
	{
	public:
		static inline EulerStep_Calculator			EulerStepCalc;
		static inline MidpointStep_Calculator		MidpointStepCalc;
		static inline RK4_FixedStep_Calculator	RK4_Basic;
		static inline RK5_CashKarp_Calculator		RK5_CashKarp;
	};
}

#endif // MML_ODE_SYSTEM_STEP_CALCULATORS_H
