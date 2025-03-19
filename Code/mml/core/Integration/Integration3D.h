#if !defined MML_INTEGRATION_3D_H
#define MML_INTEGRATION_3D_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"


namespace MML
{
	struct Integral3DInnermost : public IRealFunction
	{
		mutable Real					_currX, _currY;
		IScalarFunction<3>&		_funcToIntegrate;

		Integral3DInnermost(IScalarFunction<3>& func) 
			: _funcToIntegrate(func), _currX{ 0 }, _currY{ 0 } {}
		
		Real operator()(const Real z) const
		{
			VectorN<Real, 3> v{ _currX, _currY, z };

			return _funcToIntegrate(v);
		}
	};

	struct Integral3DInner : public IRealFunction
	{
		mutable Integral3DInnermost _fInnermost;

		IScalarFunction<3>& _funcToIntegrate;

		Real(*_zRangeLow)(Real, Real);
		Real(*_zRangeUpp)(Real, Real);

		Integral3DInner(IScalarFunction<3>& func, 
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

	struct Integral3DOuter : public IRealFunction
	{
		mutable Integral3DInner _fInner;

		IScalarFunction<3>& _funcToIntegrate;
		Real(*_yRangeLow)(Real);
		Real(*_yRangeUpp)(Real);

		Integral3DOuter(IScalarFunction<3>& func, Real yy1(Real), Real yy2(Real), 
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

	static Real Integrate3D(IScalarFunction<3>& func, 
													const Real x1, const Real x2, 
													Real y1(Real), Real y2(Real),
													Real z1(Real, Real), Real z2(Real, Real))
	{
		Integral3DOuter f1(func, y1, y2, z1, z2);

		return IntegrateGauss10(f1, x1, x2);
	}	
}

#endif