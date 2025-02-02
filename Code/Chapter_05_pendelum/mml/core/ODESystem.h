#if !defined  MML_ODE_SYSTEM_H
#define MML_ODE_SYSTEM_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

	#include "base/InterpolatedFunction.h"

namespace MML
{
	class ODESystem : public IODESystem
	{
	protected:
		int _dim;
		void (*_func)(Real, const Vector<Real>&, Vector<Real>&);

	public:
		ODESystem() : _dim(0), _func(nullptr) { }
		ODESystem(int n, void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&)) : _dim(n), _func(inFunc) { }

		int getDim() const { return _dim; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const
		{
			_func(t, x, dxdt);
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
	class ODESystemWithNumJacobian : public ODESystem
	{
	public:
		ODESystemWithNumJacobian() { }
		ODESystemWithNumJacobian(int n, void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&) ) 
			: ODESystem(n, inFunc) { }

		void jacobian(const Real t, Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
		{
			// formirati lokalnu vektorsku funkciju na bazi derivs()
//			_funcJac(t, x, dxdt, dydx);
		}
	};

	class ODESystemSolution
	{
		// first values are initial conditions
	public:
		int  _sys_dim;
		int  _count;
		Real _x1, _x2;
		Vector<Real> _xval;
		Matrix<Real> _yval;

		ODESystemSolution() {}
		ODESystemSolution(Real x1, Real x2, int dim, int maxSteps) : _sys_dim(dim), _count(maxSteps + 1), _x1(x1), _x2(x2)
		{
			_xval.resize(maxSteps + 1);
			_yval.Resize(dim, maxSteps + 1);
		}

		template<int N>
		ParametricCurveInterpolated<N> getSolutionAsParametricCurve() const
		{
			ParametricCurveInterpolated<N> curve(_xval, _yval);
			return curve;
		}

		LinearInterpRealFunc getSolutionAsLinearInterp(int component) const
		{
			Vector<Real> xsave = _xval;
			Vector<Real> ysave = _yval.VectorFromRow(component);

			return LinearInterpRealFunc(xsave, ysave);
		}

		bool Serialize(std::string fileName, std::string title) const
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "MULTI_REAL_FUNCTION_VARIABLE_SPACED" << std::endl;

			file << title << std::endl;
			file << _sys_dim << std::endl;
			file << _count << std::endl;
			file << _x1 << std::endl;
			file << _x2 << std::endl;

			for (int i = 0; i < _count; i++)
			{
				file << _xval[i] << " ";
				for (int j = 0; j < _sys_dim; j++)
				{
					file << _yval[j][i] << " ";
				}
				file << std::endl;
			}
			file.close();
			return true;
		}
		bool SerializeAsParametricCurve3D(std::string fileName, std::string title) const
		{
			if (_sys_dim != 3)
				return false;

			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "PARAMETRIC_CURVE_CARTESIAN_3D" << std::endl;
			file << "t1: " << _x1 << std::endl;
			file << "t2: " << _x2 << std::endl;
			file << "NumPoints: " << _count << std::endl;

			for (int i = 0; i < _count; i++)
			{
				file << _xval[i] << " ";
				for (int j = 0; j < _sys_dim; j++)
				{
					file << _yval[j][i] << " ";
				}
				file << std::endl;
			}

			file.close();
			return true;
		}
	};

	class ODESystemSolutionEqualSpacing
	{
		// first values are initial conditions
	public:
		int  _sys_dim;
		int  _count;
		Real x1, x2;
		Vector<Real> xval;
		Matrix<Real> yval;

		ODESystemSolutionEqualSpacing() {}
		ODESystemSolutionEqualSpacing(int dim, int numSteps) : _sys_dim(dim), _count(numSteps + 1)
		{
			xval.resize(numSteps + 1);
			yval.Resize(dim, numSteps + 1);
		}
	};
}
#endif