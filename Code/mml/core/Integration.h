#if !defined MML_INTEGRATION_H
#define MML_INTEGRATION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/Matrix.h"

namespace MML
{
	enum IntegrationMethod { TRAP, SIMPSON, ROMBERG, GAUSS10 };


	static Real TrapRefine(const IRealFunction& func, const Real a, const Real b, const int n)
	{
		// This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input
		// as a pointer to the function to be integrated between limits a and b, also input. When called with
		// n=1, the routine returns the crudest estimate of Rab f(x)dx. Subsequent calls with n=2,3,...
		// (in that sequential order) will improve the accuracy by adding 2n-2 additional interior points.
		Real x, tnm, sum, del;
		static Real s;
		int it, j;

		if (n == 1) {
			return (s = 0.5 * (b - a) * (func(a) + func(b)));
		}
		else
		{
			for (it = 1, j = 1; j < n - 1; j++)
				it <<= 1;

			tnm = it;
			del = (b - a) / tnm;
			x = a + 0.5 * del;

			for (sum = 0.0, j = 0; j < it; j++, x += del)
				sum += func(x);

			s = 0.5 * (s + (b - a) * sum / tnm);

			return s;
		}
	}

	static Real IntegrateTrap(const IRealFunction& func, const Real a, const Real b, Real req_eps, Real* achieved_precision)
	{
		// Returns the integral of the function func from a to b. The parameters EPS can be set to the
		// desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
		// number of steps. Integration is performed by the trapezoidal rule.

		// Unsophisticated as it is, routine qtrap is in fact a fairly robust way of doing
		// integrals of functions that are not very smooth. Increased sophistication will usually
		// translate into a higher-order method whose efficiency will be greater only for
		// sufficiently smooth integrands. qtrap is the method of choice, e.g., for an integrand
		// which is a function of a variable that is linearly interpolated between measured data
		// points. Be sure that you do not require too stringent an EPS, however: If qtrap takes
		// too many steps in trying to achieve your required accuracy, accumulated roundoff
		// errors may start increasing, and the routine may never converge. 
		// Value 1e-6 is just on the edge of trouble for most 32-bit machines; it is achievable when the
		// convergence is moderately rapid, but not otherwise.
		int j;
		Real s, olds = 0.0;
		Real diff = 0.0, threshold = 0.0;

		for (j = 0; j < Defaults::IntegrateTrapMaxSteps; j++)
		{
			s = TrapRefine(func, a, b, j + 1);

			if (j > 5)
			{
				diff = s - olds;
				threshold = req_eps * std::abs(olds);
				if (std::abs(diff) < threshold || (s == 0.0 && olds == 0.0))
				{
					// std::cout << "\ns : " << s << " olds : " << olds <<  " diff : " << diff << " threshold : " << threshold << std::endl;
					if (achieved_precision != nullptr)
						*achieved_precision = std::abs(diff);

					return s;
				}
			}
			olds = s;
		}
		if (achieved_precision != nullptr)
			*achieved_precision = std::abs(diff);

		return s;
	}
	static Real IntegrateTrap(const IRealFunction& func, const Real a, const Real b)
	{
		return IntegrateTrap(func, a, b, Defaults::IntegrateTrapEPS, nullptr);
	}


	static Real IntegrateGauss10(const IRealFunction& func, const Real a, const Real b)
	{
		// Returns the integral of the function func between a and b, by ten-point GaussLegendre integration: 
		// the function is evaluated exactly ten times at interior points in the range of integration.        
		static const Real x[] = { 0.1488743389816312,0.4333953941292472,
														0.6794095682990244,0.8650633666889845,0.9739065285171717 };
		static const Real w[] = { 0.2955242247147529,0.2692667193099963,
														0.2190863625159821,0.1494513491505806,0.0666713443086881 };
		Real xm = 0.5 * (b + a);
		Real xr = 0.5 * (b - a);
		Real s = 0;
		for (int j = 0; j < 5; j++) {
			Real dx = xr * x[j];
			s += w[j] * (func(xm + dx) + func(xm - dx));
		}
		return s *= xr;
	}

	struct SurfIntf2 : public IRealFunction
	{
		mutable Real xsav;
		IScalarFunction<2>& funcToInt;

		SurfIntf2(IScalarFunction<2>& func) : funcToInt(func) {}
		Real operator()(const Real y) const
		{
			VectorN<Real, 2> v{ xsav,y };
			return funcToInt(v);
		}
	};

	struct SurfIntf1 : public IRealFunction
	{
		mutable SurfIntf2 f2;
		IntegrationMethod method;

		IScalarFunction<2>& funcToInt;
		Real(*y1)(Real);
		Real(*y2)(Real);

		SurfIntf1(IScalarFunction<2>& func, IntegrationMethod inMethod, Real yy1(Real), Real yy2(Real)) : y1(yy1), y2(yy2), f2(func), funcToInt(func), method(inMethod)
		{}

		Real operator()(const Real x) const
		{
			f2.xsav = x;
			switch (method)
			{
				//case SIMPSON:
				//	return IntegrateSimpson(f2, y1(x), y2(x));
				//case ROMBERG:
				//	return IntegrateRomberg(f2, y1(x), y2(x));
			case GAUSS10:
				return IntegrateGauss10(f2, y1(x), y2(x));
			default:
				Real ret = IntegrateTrap(f2, y1(x), y2(x));
				//std::cout << "IntegrateTrap: " << ret << std::endl;
				return ret;
			}
		}
	};

	static Real IntegrateSurface(IScalarFunction<2>& func, IntegrationMethod method, const Real x1, const Real x2, Real y1(Real), Real y2(Real))
	{
		SurfIntf1 f1(func, method, y1, y2);

		switch (method)
		{
			//case SIMPSON:
			//	return IntegrateSimpson(f1, x1, x2);
			//case ROMBERG:
			//	return IntegrateRomberg(f1, x1, x2);
		case GAUSS10:
			return IntegrateGauss10(f1, x1, x2);
		default:
			return IntegrateTrap(f1, x1, x2);
		}
	}

	struct VolIntf3 : public IRealFunction
	{
		mutable Real xsav, ysav;
		IScalarFunction<3>& funcToInt;

		VolIntf3(IScalarFunction<3>& func) : funcToInt(func), xsav{ 0 }, ysav{ 0 } {}
		Real operator()(const Real z) const
		{
			VectorN<Real, 3> v{ xsav,ysav,z };
			return funcToInt(v);
		}
	};
	struct VolIntf2 : public IRealFunction
	{
		mutable VolIntf3 f3;

		IScalarFunction<3>& funcToInt;
		Real(*z1)(Real, Real);
		Real(*z2)(Real, Real);

		VolIntf2(IScalarFunction<3>& func, Real zz1(Real, Real), Real zz2(Real, Real)) : z1(zz1), z2(zz2), funcToInt(func), f3(func) {}

		Real operator()(const Real y) const
		{
			f3.ysav = y;
			return IntegrateGauss10(f3, z1(f3.xsav, y), z2(f3.xsav, y));
		}
	};
	struct VolIntf1 : public IRealFunction
	{
		mutable VolIntf2 f2;

		IScalarFunction<3>& funcToInt;
		Real(*y1)(Real);
		Real(*y2)(Real);

		VolIntf1(IScalarFunction<3>& func, Real yy1(Real), Real yy2(Real), Real z1(Real, Real),
			Real z2(Real, Real)) : y1(yy1), y2(yy2), f2(func, z1, z2), funcToInt(func)
		{}

		Real operator()(const Real x) const
		{
			f2.f3.xsav = x;
			return IntegrateGauss10(f2, y1(x), y2(x));
		}
	};

	static Real IntegrateVolume(IScalarFunction<3>& func, const Real x1, const Real x2, Real y1(Real), Real y2(Real),
		Real z1(Real, Real), Real z2(Real, Real))
	{
		VolIntf1 f1(func, y1, y2, z1, z2);

		return IntegrateGauss10(f1, x1, x2);
	}

	static inline Real(*Integrate)(const MML::IRealFunction& f, Real a, Real b, Real req_eps, Real* achieved_precision) = IntegrateTrap;

} // end namespace
#endif