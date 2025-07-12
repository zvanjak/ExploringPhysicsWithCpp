#if !defined MML_IODE_SYSTEM_STEPPER_H
#define MML_IODE_SYSTEM_STEPPER_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Vector.h"
#include "base/Matrix.h"

namespace MML
{
	class IODESystemStepper	{
	public:
		virtual void doStep(const IODESystem& odeSystem, 
												const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
												const Real h, Vector<Real>& xout, Vector<Real>& xerr) = 0;
	};
}
#endif
