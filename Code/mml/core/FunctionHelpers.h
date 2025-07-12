#if !defined  MML_FUNCTION_HELPERS_H
#define MML_FUNCTION_HELPERS_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/VectorN.h"
#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/Polynom.h"

#include "core/Derivation.h"

namespace MML
{

	static PolynomRealFunc TaylorSeries2(IRealFunction& f, Real a)
	{
		PolynomRealFunc ret(2);

		Real val = f(a);
		Real coef1 = Derivation::NDer6(f, a);
		Real coef2 = Derivation::NSecDer4(f, a) / 2.0;

		ret[0] = val - coef1 * a + coef2 * POW2(a);
		ret[1] = coef1 - 2.0 * coef2 * a;
		ret[2] = coef2;

		return ret;
	}
	static PolynomRealFunc TaylorSeries3(IRealFunction& f, Real a)
	{
		PolynomRealFunc ret(3);

		Real val = f(a);
		Real coef1 = Derivation::NDer6(f, a);
		Real coef2 = Derivation::NSecDer4(f, a) / 2.0;
		Real coef3 = Derivation::NThirdDer2(f, a) / 6.0;

		ret[0] = val - coef1 * a + coef2 * POW2(a) - coef3 * pow(a, 3);
		ret[1] = coef1 - 2.0 * coef2 * a + 3.0 * coef3 * POW2(a);
		ret[2] = coef2 - 3.0 * coef3 * a;
		ret[3] = coef3;

		return ret;
	}

	/////////////////////       FUNCTION HELPERS         //////////////////////////////////
	class RealFuncDerived1 : public IRealFunction
	{
		IRealFunction& _f1;
		Real _deriv_step;
	public:
		RealFuncDerived1(IRealFunction& f1) : _f1(f1), _deriv_step(0.0) {}
		RealFuncDerived1(IRealFunction& f1, Real deriv_step) : _f1(f1), _deriv_step(deriv_step) {}

		Real operator()(Real x) const {
			if (_deriv_step != 0.0)
				return Derivation::NDer1(_f1, x, _deriv_step);
			else
				return Derivation::NDer1(_f1, x);
		}
	};
	class RealFuncDerived2 : public IRealFunction
	{
		IRealFunction& _f1;
		Real _deriv_step;
	public:
		RealFuncDerived2(IRealFunction& f1) : _f1(f1), _deriv_step(0.0) {}
		RealFuncDerived2(IRealFunction& f1, Real deriv_step) : _f1(f1), _deriv_step(deriv_step) {}

		Real operator()(Real x) const {
			if (_deriv_step != 0.0)
				return Derivation::NDer2(_f1, x, _deriv_step);
			else
				return Derivation::NDer2(_f1, x);
		}
	};
	class RealFuncDerived4 : public IRealFunction
	{
		IRealFunction& _f1;
		Real _deriv_step;
	public:
		RealFuncDerived4(IRealFunction& f1) : _f1(f1), _deriv_step(0.0) {}
		RealFuncDerived4(IRealFunction& f1, Real deriv_step) : _f1(f1), _deriv_step(deriv_step) {}

		Real operator()(Real x) const {
			if (_deriv_step != 0.0)
				return Derivation::NDer4(_f1, x, _deriv_step);
			else
				return Derivation::NDer4(_f1, x);
		}
	};
	class RealFuncDerived6 : public IRealFunction
	{
		IRealFunction& _f1;
		Real _deriv_step;
	public:
		RealFuncDerived6(IRealFunction& f1) : _f1(f1), _deriv_step(0.0) {}
		RealFuncDerived6(IRealFunction& f1, Real deriv_step) : _f1(f1), _deriv_step(deriv_step) {}

		Real operator()(Real x) const {
			if (_deriv_step != 0.0)
				return Derivation::NDer6(_f1, x, _deriv_step);
			else
				return Derivation::NDer6(_f1, x);
		}
	};
	class RealFuncDerived8 : public IRealFunction
	{
		IRealFunction& _f1;
		Real _deriv_step;
	public:
		RealFuncDerived8(IRealFunction& f1) : _f1(f1), _deriv_step(0.0) {}
		RealFuncDerived8(IRealFunction& f1, Real deriv_step) : _f1(f1), _deriv_step(deriv_step) {}

		Real operator()(Real x) const {
			if (_deriv_step != 0.0)
				return Derivation::NDer8(_f1, x, _deriv_step);
			else
				return Derivation::NDer8(_f1, x);
		}
	};

	// TODO - integrate helper - ima i x1, a x2 je parametar funkcije

	class RealFuncDiffHelper : public IRealFunction
	{
		IRealFunction& _f1, & _f2;
	public:
		RealFuncDiffHelper(IRealFunction& f1, IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return _f1(x) - _f2(x); }
	};

	class RealFuncDiffAbsHelper : public IRealFunction
	{
		IRealFunction& _f1, & _f2;
	public:
		RealFuncDiffAbsHelper(IRealFunction& f1, IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return std::abs(_f1(x) - _f2(x)); }
	};

	class RealFuncDiffSqrHelper : public IRealFunction
	{
		IRealFunction& _f1, & _f2;
	public:
		RealFuncDiffSqrHelper(IRealFunction& f1, IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return POW2(_f1(x) - _f2(x)); }
	};

} // end namespace

#endif