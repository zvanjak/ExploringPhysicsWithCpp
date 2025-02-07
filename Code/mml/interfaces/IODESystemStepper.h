#if !defined MML_IODE_SYSTEM_STEPPER_H
#define MML_IODE_SYSTEM_STEPPER_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Vector.h"
#include "base/Matrix.h"

namespace MML
{
	class IODESystemStepperCalculator	{
	public:
		virtual void doStep(const IODESystem& odeSystem, 
												const Real x, const Vector<Real>& y_start, const Vector<Real>& dydx,
												const Real h, Vector<Real>& yout, Vector<Real>& yerr) = 0;
	};
}
#endif
