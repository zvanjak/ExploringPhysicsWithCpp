#if !defined  MML_FUNCTION_H
#define MML_FUNCTION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/VectorN.h"
#include "base/Vector.h"
#include "base/Matrix.h"

namespace MML
{
	///////////////////////////     REAL FUNCTION      ////////////////////////////////////
	class RealFunction : public IRealFunction
	{
		Real(*_func)(const Real);
	public:
		RealFunction(Real(*inFunc)(const Real)) : _func(inFunc) {}

		Real operator()(const Real x) const { return _func(x); }
	};
	class RealFunctionFromStdFunc : public IRealFunction
	{
		std::function<Real(const Real)> _func;
	public:
		RealFunctionFromStdFunc(std::function<Real(const Real)> inFunc) : _func(inFunc) {}

		Real operator()(const Real x) const { return _func(x); }
	};

	///////////////////////////     SCALAR FUNCTION       //////////////////////////////////
	template<int N>
	class ScalarFunction : public IScalarFunction<N>
	{
		Real(*_func)(const VectorN<Real, N>&);
	public:
		ScalarFunction(Real(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		Real operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	template<int N>
	class ScalarFunctionFromStdFunc : public IScalarFunction<N>
	{
		std::function<Real(const VectorN<Real, N>&)> _func;
	public:
		ScalarFunctionFromStdFunc(std::function<Real(const VectorN<Real, N>&)> inFunc) : _func(inFunc) {}

		Real operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	/////////////////////////    VECTOR FUNCTION N -> N      ///////////////////////////////////
	template<int N>
	class VectorFunction : public IVectorFunction<N>
	{
		VectorN<Real, N>(*_func)(const VectorN<Real, N>&);
	public:
		VectorFunction(VectorN<Real, N>(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		VectorN<Real, N>     operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	template<int N>
	class VectorFunctionFromStdFunc : public IVectorFunction<N>
	{
		std::function<VectorN<Real, N>(const VectorN<Real, N>&)> _func;
	public:
		VectorFunctionFromStdFunc(std::function<VectorN<Real, N>(const VectorN<Real, N>&)>& inFunc) : _func(inFunc) {}

		VectorN<Real, N>     operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	/////////////////////////    VECTOR FUNCTION N -> M      ///////////////////////////////////
	template<int N, int M>
	class VectorFunctionNM : public IVectorFunctionNM<N, M>
	{
		VectorN<Real, M>(*_func)(const VectorN<Real, N>&);
	public:
		VectorFunctionNM(VectorN<Real, M>(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		VectorN<Real, M>     operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	template<int N, int M>
	class VectorFunctionNMFromStdFunc : public IVectorFunctionNM<N, M>
	{
		std::function<VectorN<Real, M>(const VectorN<Real, N>&)> _func;
	public:
		VectorFunctionNMFromStdFunc(std::function<VectorN<Real, M>(const VectorN<Real, N>&)>& inFunc) : _func(inFunc) {}

		VectorN<Real, M>     operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	//////////////////////     PARAMETRIC CURVE             ///////////////////////////////////
	template<int N>
	class   ParametricCurve : public IParametricCurve<N>
	{
		Real _minT;
		Real _maxT;
		VectorN<Real, N>(*_func)(Real);
	public:
		ParametricCurve(VectorN<Real, N>(*inFunc)(Real)) : _func(inFunc), _minT(Constants::NegativeInf), _maxT(Constants::PositiveInf) {}
		ParametricCurve(Real minT, Real maxT, VectorN<Real, N>(*inFunc)(Real)) : _func(inFunc), _minT(minT), _maxT(maxT) {}

		Real getMinT() const { return _minT; }
		Real getMaxT() const { return _maxT; }

		virtual VectorN<Real, N> operator()(Real x) const { return _func(x); }
	};

	template<int N>
	class ParametricCurveFromStdFunc : public IParametricCurve<N>
	{
		Real _minT;
		Real _maxT;
		std::function<VectorN<Real, N>(Real)> _func;
	public:
		ParametricCurveFromStdFunc(std::function<VectorN<Real, N>(Real)>& inFunc) : _func(inFunc), _minT(Constants::NegativeInf), _maxT(Constants::PositiveInf) {}
		ParametricCurveFromStdFunc(Real minT, Real maxT, std::function<VectorN<Real, N>(Real)>& inFunc) : _func(inFunc), _minT(minT), _maxT(maxT) {}

		Real getMinT() const { return _minT; }
		Real getMaxT() const { return _maxT; }

		VectorN<Real, N> operator()(Real x) const { return _func(x); }
	};

	template<int N>
	class   ParametricCurveParametrized : public IParametricCurveParametrized<N>
	{
		Real _minT;
		Real _maxT;
		Vector<Real> _params;
		VectorN<Real, N>(*_func)(Real);
	public:
		ParametricCurveParametrized(VectorN<Real, N>(*inFunc)(Real), int numParams) : _func(inFunc), _minT(Constants::NegativeInf), _maxT(Constants::PositiveInf) 
		{
			_params.resize(numParams);
		}
		ParametricCurveParametrized(Real minT, Real maxT, VectorN<Real, N>(*inFunc)(Real), int numParams) : _func(inFunc), _params(numParams), _minT(minT), _maxT(maxT) {}

		Real getMinT() const { return _minT; }
		Real getMaxT() const { return _maxT; }

		virtual VectorN<Real, N> operator()(Real x) const { return _func(x); }

		virtual int		getNumParam() const { return _params.size(); }
		virtual Real	getParam(int i) const { return _params[i]; }
		virtual void	setParam(int i, Real val) { _params[i] = val; }

		virtual VectorN<Real, N> operator()(Real t, const Vector<Real>& params)
		{
			_params = params;
			return _func(t);
		}

	};

	/////////////////////       PARAMETRIC SURFACE         //////////////////////////////////
	template<int N>
	class ParametricSurface : public IParametricSurfaceRect<N>
	{
		// TODO - ensure that N is at least 3!!!
		Real _minX;
		Real _maxX;
		Real _minY;
		Real _maxY;
		VectorN<Real, N>(*_func)(Real u, Real w);

	public:
		ParametricSurface(VectorN<Real, N>(*inFunc)(Real u, Real w)) : _func(inFunc), _minX(Constants::NegativeInf), _maxX(Constants::PositiveInf), _minY(Constants::NegativeInf), _maxY(Constants::PositiveInf) {}
		ParametricSurface(VectorN<Real, N>(*inFunc)(Real u, Real w), Real minX, Real maxX, Real minY, Real maxY) : _func(inFunc), _minX(minX), _maxX(maxX), _minY(minY), _maxY(maxY) {}

		VectorN<Real, N> operator()(Real u, Real w) const { return _func(u, w); }

		virtual Real getMinU() const { return _minX; }
		virtual Real getMaxU() const { return _maxX; }
		virtual Real getMinW() const { return _minY; }
		virtual Real getMaxW() const { return _maxY; }
	};

	// imati cemo i surface Discrete, kreiran od triangles?

	template<int N>
	class ParametricSurfaceFromStdFunc : public IParametricSurfaceRect<N>
	{
		Real _minX;
		Real _maxX;
		Real _minY;
		Real _maxY;
		std::function<VectorN<Real, N>(Real u, Real w)> _func;
	public:
		ParametricSurfaceFromStdFunc(std::function<VectorN<Real, N>(Real u, Real w)>& inFunc) : _func(inFunc), _minX(Constants::NegativeInf), _maxX(Constants::PositiveInf), _minY(Constants::NegativeInf), _maxY(Constants::PositiveInf) {}
		ParametricSurfaceFromStdFunc(std::function<VectorN<Real, N>(Real u, Real w)>& inFunc, Real minX, Real maxX, Real minY, Real maxY) : _func(inFunc), _minX(minX), _maxX(maxX), _minY(minY), _maxY(maxY) {}

		VectorN<Real, N> operator()(Real u, Real w) const { return _func(u, w); }

		virtual Real getMinU() const { return _minX; }
		virtual Real getMaxU() const { return _maxX; }
		virtual Real getMinW() const { return _minY; }
		virtual Real getMaxW() const { return _maxY; }
	};
} // end namespace

#endif