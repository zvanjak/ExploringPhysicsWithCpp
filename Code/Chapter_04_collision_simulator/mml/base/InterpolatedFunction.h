#if !defined MML_INTERPOLATEDFUNCTION_H
#define MML_INTERPOLATEDFUNCTION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "base/Function.h"

namespace MML
{
	class RealFunctionInterpolated : public IRealFunction
	{
		double _xmin, _xmax;

	public:
		Real virtual rawinterp(int jlo, Real x) const = 0;
	};

	class RealFunctionInterpolatedBase : public RealFunctionInterpolated
	{
		// TODO - HIGH, LAKO, sve privatizirati!
	public:
		mutable int jsav, cor;

		int n, mm, dj;
		Vector<Real> xx, yy;

		RealFunctionInterpolatedBase(const Vector<Real>& x, const Vector<Real>& y, int m)
			: n((int)x.size()), mm(m), jsav(0), cor(0), xx(x), yy(y)
		{
			dj = std::min(1, (int)pow((Real)n, 0.25));
		}

		Real operator()(Real x) const
		{
			int jlo = cor ? hunt(x) : locate(x);
			return rawinterp(jlo, x);
		}

		// Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
		// xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
		// increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
		int locate(const Real x) const
		{
			int ju, jm, jl;
			if (n < 2 || mm < 2 || mm > n) throw("locate size error");
			bool ascnd = (xx[n - 1] >= xx[0]);
			jl = 0;
			ju = n - 1;
			while (ju - jl > 1) {
				jm = (ju + jl) >> 1;
				if (x >= xx[jm] == ascnd)
					jl = jm;
				else
					ju = jm;
			}
			cor = std::abs(jl - jsav) > dj ? 0 : 1;
			jsav = jl;
			return std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1)));
		}

		// Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
		// xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
		// increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
		int hunt(const Real x) const
		{
			int jl = jsav, jm, ju, inc = 1;
			if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
			bool ascnd = (xx[n - 1] >= xx[0]);
			if (jl < 0 || jl > n - 1) {
				jl = 0;
				ju = n - 1;
			}
			else {
				if (x >= xx[jl] == ascnd) {
					for (;;) {
						ju = jl + inc;
						if (ju >= n - 1) { ju = n - 1; break; }
						else if (x < xx[ju] == ascnd) break;
						else {
							jl = ju;
							inc += inc;
						}
					}
				}
				else {
					ju = jl;
					for (;;) {
						jl = jl - inc;
						if (jl <= 0) { jl = 0; break; }
						else if (x >= xx[jl] == ascnd) break;
						else {
							ju = jl;
							inc += inc;
						}
					}
				}
			}
			while (ju - jl > 1) {
				jm = (ju + jl) >> 1;
				if (x >= xx[jm] == ascnd)
					jl = jm;
				else
					ju = jm;
			}
			cor = std::abs(jl - jsav) > dj ? 0 : 1;
			jsav = jl;
			return std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1)));
		}
	};

	struct LinearInterpRealFunc : RealFunctionInterpolatedBase
	{
		LinearInterpRealFunc(Vector<Real>& xv, Vector<Real>& yv) : RealFunctionInterpolatedBase(xv, yv, 2) {}

		Real rawinterp(int j, Real x) const {
			if (xx[j] == xx[j + 1]) return yy[j];
			else return yy[j] + ((x - xx[j]) / (xx[j + 1] - xx[j])) * (yy[j + 1] - yy[j]);
		}
	};


	template<int N>
	class ParametricCurveInterpolated : public IParametricCurve<N>
	{
		// TODO 0.9 - HIGH, verify
		Real _minT;
		Real _maxT;
		std::shared_ptr<Vector<Real>> _xvals;
		std::shared_ptr<Matrix<Real>> _yvals;

	public:
		ParametricCurveInterpolated() {}
		ParametricCurveInterpolated(std::shared_ptr<Vector<Real>> xsave, std::shared_ptr<Matrix<Real>> ysave) : _xvals(xsave), _yvals(ysave) {}
		ParametricCurveInterpolated(int inPoints, Real* xsave, Real* ysave)
		{
			_minT = xsave[0];
			_maxT = xsave[inPoints - 1];

			_xvals = std::make_shared<Vector<Real>>(inPoints);
			_yvals = std::make_shared<Matrix<Real>>(N, inPoints);

			for (int i = 0; i < inPoints; i++)
			{
				(*_xvals)[i] = xsave[i];
				for (int j = 0; j < N; j++)
					(*_yvals)(j, i) = ysave[i * N + j];
			}
		}
		ParametricCurveInterpolated(const Vector<Real>& xsave, const Matrix<Real>& ysave)
		{
			_minT = xsave[0];
			_maxT = xsave[xsave.size() - 1];

			_xvals = std::make_shared<Vector<Real>>(xsave);
			_yvals = std::make_shared<Matrix<Real>>(ysave);
		}

		Real getMinT() const { return _minT; }
		Real getMaxT() const { return _maxT; }

		virtual VectorN<Real, N> operator()(Real x) const { return VectorN<Real, N>{0}; }
	};

	template<int N>
	class InterpolatedSurface : public IParametricSurfaceRect<N>
	{
	public:
		InterpolatedSurface() {}

		VectorN<Real, N> operator()(const VectorN<Real, 2>& x) const { return VectorN<Real, N>{}; }
	};
}

#endif