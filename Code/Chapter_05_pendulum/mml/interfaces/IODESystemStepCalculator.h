#if !defined MML_IODE_SYSTEM_STEP_CALC_H
#define MML_IODE_SYSTEM_STEP_CALC_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Vector.h"
#include "base/Matrix.h"

namespace MML
{
	class IODESystemStepCalculator
	{
	public:
		virtual void calcStep(const IODESystem& odeSystem, 
													const Real x, const Vector<Real>& y_start, const Vector<Real>& dydx,
													const Real h, Vector<Real>& yout, Vector<Real>& yerr) const = 0;
	};
}
#endif
