#if !defined  MML_ODE_SYSTEM_H
#define MML_ODE_SYSTEM_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

//#include "core/Derivation.h"

namespace MML
{
	class ODESystem : public IODESystem
	{
	protected:
		int _dim;
		void (*_func)(Real, const Vector<Real>&, Vector<Real>&);

	public:
		ODESystem() 
				: _dim(0), _func(nullptr) { }
		ODESystem(int n, void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&)) 
				: _dim(n), _func(inFunc) { }

		int		getDim() const { return _dim; }
		void	derivs(const Real t, const Vector<Real> &x, Vector<Real> &dxdt) const
		{
			_func(t, x, dxdt);
		}

		void  operator()(const Real t, const Vector<Real>& x, Vector<Real>& dxdt)
		{
			derivs(t, x, dxdt);
		}
	};

	class ODESystemWithJacobian : public ODESystem
	{
	private:
		void (*_funcJac)(const Real, const Vector<Real>&, Vector<Real>&, Matrix<Real>&);

	public:
		ODESystemWithJacobian() : _funcJac(nullptr) { }
		ODESystemWithJacobian(int n,
			void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&),
			void (*inFuncJac)(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
		) : ODESystem(n, inFunc), _funcJac(inFuncJac) { }

		void jacobian(const Real t, Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
		{
			if (_funcJac != nullptr)
				_funcJac(t, x, dxdt, dydx);
		}
	};

	// calculates needed Jacobian numerically
//	class ODESystemWithNumJacobian : public ODESystem
//	{
//	public:
//		ODESystemWithNumJacobian() { }
//		ODESystemWithNumJacobian(int n, void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&) ) 
//			: ODESystem(n, inFunc) { }
//
//		void jacobian(const Real t, Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
//		{
//			// formirati lokalnu vektorsku funkciju na bazi derivs()
////			_funcJac(t, x, dxdt, dydx);
//		}
//	};

}
#endif