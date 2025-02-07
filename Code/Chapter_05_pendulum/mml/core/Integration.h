#if !defined MML_INTEGRATION_H
#define MML_INTEGRATION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/Matrix.h"

namespace MML
{

	enum IntegrationMethod { TRAP, SIMPSON, ROMBERG, GAUSS10 };

	class IQuadrature {
	protected:
		int _currStep = 0;

	public:
		// Returns the nth stage of refinement of the extended trapezoidal rule.
		virtual Real next() = 0;
	};

	class TrapIntegrator : IQuadrature
	{
	public:
		Real _a, _b, _currSum;
		const IRealFunction& _func;

		TrapIntegrator(const IRealFunction& func, Real aa, Real bb) 
			:	_func(func), _a(aa), _b(bb), _currSum(0.0)	{	}

		Real next() {
			Real x, sum, del;
			int subDivNum, j;

			_currStep++;
			if (_currStep == 1) {
				return (_currSum = 0.5 * (_b - _a) * (_func(_a) + _func(_b)));
			}
			else {
				for (subDivNum = 1, j = 1; j < _currStep - 1; j++)
					subDivNum *= 2;

				del = (_b - _a) / subDivNum;
				x = _a + 0.5 * del;

				for (sum = 0.0, j = 0; j < subDivNum; j++, x += del)
					sum += _func(x);

				_currSum = 0.5 * (_currSum + (_b - _a) * sum / subDivNum);

				return _currSum;
			}
		}
	};

	// Using extended trapezoidal rule, returns the integral of the function func from a to b.
	// Parameter eps can be	set to the desired fractional accuracy.
	static Real IntegrateTrap(const IRealFunction& func, Real a, Real b,
														int* doneSteps, Real* achievedPrec,
														const Real eps = Defaults::TrapezoidIntegrationEPS)	{
		int		j;
		Real	currSum, oldSum = 0.0;

		TrapIntegrator t(func, a, b);

		for (j = 0; j < Defaults::TrapezoidIntegrationMaxSteps; j++)
		{
			currSum = t.next();

			if (j > 5) {
				if ( std::abs(currSum - oldSum) < eps * std::abs(oldSum) ||
					   (currSum == 0.0 && oldSum == 0.0) )
				{
					if (doneSteps != nullptr)			*doneSteps = j;
					if (achievedPrec != nullptr)  *achievedPrec = std::abs(currSum - oldSum);

					return currSum;
				}
			}

			oldSum = currSum;
		}
		if (doneSteps != nullptr)			*doneSteps = j;
		if (achievedPrec != nullptr)	*achievedPrec = std::abs(currSum - oldSum);

		return currSum;
	}
	static Real IntegrateTrap(const IRealFunction& func, Real a, Real b, 
														Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		return IntegrateTrap(func, a, b, nullptr, nullptr, eps);
	}

	// Using Simpson's rule, returns the integral of the function func from a to b.
	// Parameter eps can be set to the desired fractional accuracy.
	static Real IntegrateSimpson(const IRealFunction& func, Real a, Real b, 
															 int* doneSteps, Real* achievedPrec,
															 const Real eps = Defaults::SimpsonIntegrationEPS)
	{
		int j;
		Real currSum, st, ost = 0.0, oldSum = 0.0;

		TrapIntegrator t(func, a, b);

		for (j = 0; j < Defaults::SimpsonIntegrationMaxSteps; j++)
		{
			st = t.next();

			currSum = (4.0 * st - ost) / 3.0;

			if (j > 5) {
				if (std::abs(currSum - oldSum) < eps * std::abs(oldSum) ||
					(currSum == 0.0 && oldSum == 0.0))
				{
					if (doneSteps != nullptr)			*doneSteps = j;
					if (achievedPrec != nullptr)	*achievedPrec = std::abs(currSum - oldSum);

					return currSum;
				}
			}

			oldSum = currSum;
			ost = st;
		}
		if (doneSteps != nullptr)			*doneSteps = j;
		if (achievedPrec != nullptr)	*achievedPrec = std::abs(currSum - oldSum);

		return currSum;
	}
	static Real IntegrateSimpson(const IRealFunction& func, Real a, Real b,
															 const Real eps = Defaults::SimpsonIntegrationEPS)
	{
		return IntegrateSimpson(func, a, b, nullptr, nullptr, eps);
	}

	// Using ten-point GaussLegendre integration, returns the integral of the function func 
	// between a and b (the function is evaluated exactly ten times at interior points=
	static Real IntegrateGauss10(const IRealFunction& func, const Real a, const Real b)
	{
		static const Real x[] = { 0.1488743389816312,0.4333953941292472,
														0.6794095682990244,0.8650633666889845,0.9739065285171717 };
		static const Real w[] = { 0.2955242247147529,0.2692667193099963,
														0.2190863625159821,0.1494513491505806,0.0666713443086881 };
		
		Real xm = 0.5 * (b + a);
		Real xr = 0.5 * (b - a);
		
		Real s = 0;
		for (int j = 0; j < 5; j++) 
		{
			Real dx = xr * x[j];
			s += w[j] * (func(xm + dx) + func(xm - dx));
		}
		return s *= xr;
	}

	struct SurfaceIntegralInner : public IRealFunction
	{
		mutable Real				_currX;
		IScalarFunction<2>& _funcToIntegrate;

		SurfaceIntegralInner(IScalarFunction<2>& func) : _funcToIntegrate(func) {}

		Real operator()(const Real y) const
		{
			VectorN<Real, 2> v{ _currX, y };
			return _funcToIntegrate(v);
		}
	};

	struct SurfaceIntegralOuter : public IRealFunction
	{
		mutable SurfaceIntegralInner _fInner;

		IntegrationMethod			_integrMethod;
		IScalarFunction<2>&		_funcToIntegrate;

		Real(*_yRangeLow)(Real);
		Real(*_yRangeUpp)(Real);

		SurfaceIntegralOuter(IScalarFunction<2>& func, IntegrationMethod inMethod, 
												 Real yy1(Real), Real yy2(Real)) 
				: _yRangeLow(yy1), _yRangeUpp(yy2), 
					_fInner(func), _funcToIntegrate(func), _integrMethod(inMethod)	{	}

		// for given x, will return (ie. integrate function over Y-range)
		Real operator()(const Real x) const
		{
			_fInner._currX = x;
			switch (_integrMethod)
			{
				case SIMPSON:
					return IntegrateSimpson(_fInner, _yRangeLow(x), _yRangeUpp(x));
				//case ROMBERG:
				//	return IntegrateRomberg(f2, y1(x), _secDerY(x));
				case GAUSS10:
					return IntegrateGauss10(_fInner, _yRangeLow(x), _yRangeUpp(x));
				default:
					return IntegrateTrap(_fInner, _yRangeLow(x), _yRangeUpp(x));
			}
		}
	};

	static Real IntegrateSurface(IScalarFunction<2>& func, IntegrationMethod method, 
															 const Real x1, const Real x2, 
															 Real y1(Real), Real y2(Real))
	{
		SurfaceIntegralOuter f1(func, method, y1, y2);

		switch (method)
		{
			case SIMPSON:
				return IntegrateSimpson(f1, x1, x2);
			//case ROMBERG:
			//	return IntegrateRomberg(f1, x1, x2);
			case GAUSS10:
				return IntegrateGauss10(f1, x1, x2);
			default:
				return IntegrateTrap(f1, x1, x2);
		}
	}

	struct VolumeIntegralInnermost : public IRealFunction
	{
		mutable Real					_currX, _currY;
		IScalarFunction<3>&		_funcToIntegrate;

		VolumeIntegralInnermost(IScalarFunction<3>& func) 
			: _funcToIntegrate(func), _currX{ 0 }, _currY{ 0 } {}
		
		Real operator()(const Real z) const
		{
			VectorN<Real, 3> v{ _currX, _currY, z };

			return _funcToIntegrate(v);
		}
	};

	struct VolumeIntegralInner : public IRealFunction
	{
		mutable VolumeIntegralInnermost _fInnermost;

		IScalarFunction<3>& _funcToIntegrate;

		Real(*_zRangeLow)(Real, Real);
		Real(*_zRangeUpp)(Real, Real);

		VolumeIntegralInner(IScalarFunction<3>& func, 
												Real zz1(Real, Real), Real zz2(Real, Real)) 
			: _zRangeLow(zz1), _zRangeUpp(zz2), _funcToIntegrate(func), _fInnermost(func) {}

		Real operator()(const Real y) const
		{
			_fInnermost._currY = y;

			return IntegrateGauss10(_fInnermost, 
															_zRangeLow(_fInnermost._currX, y), 
															_zRangeUpp(_fInnermost._currX, y));
		}
	};

	struct VolumeIntegralOuter : public IRealFunction
	{
		mutable VolumeIntegralInner _fInner;

		IScalarFunction<3>& _funcToIntegrate;
		Real(*_yRangeLow)(Real);
		Real(*_yRangeUpp)(Real);

		VolumeIntegralOuter(IScalarFunction<3>& func, Real yy1(Real), Real yy2(Real), 
												Real z1(Real, Real), Real z2(Real, Real)) 
			: _yRangeLow(yy1), _yRangeUpp(yy2), _funcToIntegrate(func), _fInner(func, z1, z2)
		{
		}

		Real operator()(const Real x) const
		{
			_fInner._fInnermost._currX = x;

			return IntegrateGauss10(_fInner, _yRangeLow(x), _yRangeUpp(x));
		}
	};

	static Real IntegrateVolume(IScalarFunction<3>& func, 
															const Real x1, const Real x2, 
															Real y1(Real), Real y2(Real),
															Real z1(Real, Real), Real z2(Real, Real))
	{
		VolumeIntegralOuter f1(func, y1, y2, z1, z2);

		return IntegrateGauss10(f1, x1, x2);
	}

	static inline Real(*Integrate)(const MML::IRealFunction& f, Real a, Real b, int* done_steps, Real* achieved_precision, Real req_eps) = IntegrateTrap;

} // end namespace
#endif