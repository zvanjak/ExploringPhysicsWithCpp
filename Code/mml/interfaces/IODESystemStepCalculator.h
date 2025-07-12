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
													const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
													const Real h, Vector<Real>& xout, Vector<Real>& xerr) const = 0;
	};
}
#endif
