#if !defined MML_FUNCTION_ANALYZER_H
#define MML_FUNCTION_ANALYZER_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/VectorN.h"

#include "core/Integration.h"
#include "core/FunctionHelpers.h"

#include "algorithms/RootFinding.h"


namespace MML
{
	// TODO - MED, function point analyzer at point and interval - is continuous, is derivation defined
	class RealFunctionAnalyzer
	{
		const IRealFunction& _f;
		std::string _funcName;
	public:
		RealFunctionAnalyzer(const IRealFunction& f) : _f(f) {}
		RealFunctionAnalyzer(const IRealFunction& f, std::string inName) : _f(f), _funcName(inName) {}

		void PrintPointAnalysis(Real x, Real eps = 1e-6)
		{
			std::cout << "Function analysis at point: " << x << ":" << std::endl;
			std::cout << "  Defined at point:    " << (isDefinedAtPoint(x) ? "yes" : "no") << std::endl;
			std::cout << "  Continuous at point: " << (isContinuousAtPoint(x, eps) ? "yes" : "no") << std::endl;
			std::cout << "  Is inflection point: " << (isInflectionPoint(x, eps) ? "yes" : "no") << std::endl;
		}

		void PrintIntervalAnalysis(Real x1, Real x2, int numPoints, Real eps = 1e-6)
		{
			if (_funcName != "")
				std::cout << std::fixed << "f(x) = " << _funcName << " - Function analysis in interval [" << x1 << ", " << x2 << "] with " << numPoints << " points:" << std::endl;
			else
				std::cout << std::fixed << "Function analysis in interval [" << x1 << ", " << x2 << "] with " << numPoints << " points:" << std::endl;

			bool isDef = true;
			bool isCont = true;
			std::vector<Real> _notDefinedPoints;
			std::vector<Real> _notContinuousPoints;

			for (int i = 0; i < numPoints; i++)
			{
				Real x = x1 + i * (x2 - x1) / numPoints;
				if (!isDefinedAtPoint(x))
				{
					isDef = false;
					_notDefinedPoints.push_back(x);
				}
				else if (!isContinuousAtPoint(x, eps))
				{
					isCont = false;
					_notContinuousPoints.push_back(x);
				}
			}

			std::cout << "  Defined    : " << (isDef ? "yes" : "no");
			if (!_notDefinedPoints.empty())
			{
				std::cout << "  Not defined at points: ";
				for (int i = 0; i < _notDefinedPoints.size(); i++)
					std::cout << _notDefinedPoints[i] << " ";
				std::cout << std::endl;
			}
			else
				std::cout << std::endl;
			std::cout << "  Continuous : " << (isCont ? "yes" : "no");
			if (!_notContinuousPoints.empty())
			{
				std::cout << "  Not continuous at points: ";
				for (int i = 0; i < _notContinuousPoints.size(); i++)
					std::cout << _notContinuousPoints[i] << " ";
				std::cout << std::endl;
			}
			else
				std::cout << std::endl;
			std::cout << "  Monotonic  : " << (isMonotonic(x1, x2, numPoints) ? "yes" : "no") << std::endl;
			std::cout << "  Min        : " << MinInNPoints(x1, x2, numPoints) << std::endl;
			std::cout << "  Max        : " << MaxInNPoints(x1, x2, numPoints) << std::endl;
		}
		void PrintDetailedIntervalAnalysis(Real x1, Real x2, int numPoints, Real eps = 1e-6)
		{
			std::cout << std::fixed << "Function analysis in interval [" << x1 << ", " << x2 << "] with " << numPoints << " points:" << std::endl;
			std::cout << " Point " << "           Value     " << "         First.der.     " << "     Sec.der.     " << "Defined " << " Continuous " << " Inflection " << std::endl;
			for (int i = 0; i < numPoints; i++)
			{
				Real x = x1 + i * (x2 - x1) / numPoints;

				std::cout << std::setw(6) << x << " :  ";
				std::cout << std::setw(16) << _f(x) << " :  ";
				std::cout << std::setw(16) << Derivation::NDer1(_f, x) << " :  ";
				std::cout << std::setw(16) << Derivation::NSecDer2(_f, x) << " :  ";
				std::cout << (isDefinedAtPoint(x) ? "   yes   " : "    no   ");
				std::cout << (isContinuousAtPoint(x, eps) ? "   yes   " : "    no   ");
				std::cout << (isInflectionPoint(x, eps) ? "     yes   " : "      no   ") << std::endl;
			}
		}

		std::vector<Real> GetRoots(Real x1, Real x2, Real eps)
		{
			std::vector<Real> roots;
			Real step = (x2 - x1) / 1000;
			Real prev = _f(x1);
			for (int i = 1; i < 1000; i++)
			{
				Real curr = _f(x1 + i * step);
				if (prev * curr < 0)
				{
					Real root = RootFinding::FindRootBisection(_f, x1 + (i - 1) * step, x1 + i * step, eps);
					roots.push_back(root);
				}
				prev = curr;
			}
			return roots;
		}
		Vector<Real> GetLocalOptimums(Real x1, Real x2)
		{
			Vector<Real> optimums;

			return optimums;
		}
		Vector<Real> GetInflectionPoints(Real x1, Real x2)
		{
			Vector<Real> inflection_points;

			return inflection_points;
		}
		
		bool isDefinedAtPoint(Real x)
		{
			Real y = _f(x);
			return !std::isnan(y) && !std::isinf(y);
		}
		bool isContinuousAtPoint(Real x, Real eps)
		{
			if (!isDefinedAtPoint(x))
				return false;

			Real h = eps;
			Real val = _f(x);
			Real left = _f(x - h);
			Real right = _f(x + h);

			// handling case of constant function
			if (val == left && val == right)
				return true;

			Real abs_dif = std::abs(left - right);

			// smanji tu abs razliku na pola, i nadji h za koji to vrijedi
			Real req_new_abs = abs_dif / 2;
			h /= 2.0;
			while (h > eps / 1000.0)
			{
				left = _f(x - h);
				right = _f(x + h);

				if (std::abs(left - right) < req_new_abs)
					return true;

				h /= 2.0;
			}
			return false;
		}
		bool isLocalOptimum(Real x, Real eps)
		{
			Real left_sec_der = Derivation::NSecDer4(_f, x - 4 * eps);
			Real right_sec_der = Derivation::NSecDer4(_f, x + 4 * eps);

			return left_sec_der * right_sec_der > 0;
		}
		bool isInflectionPoint(Real x, Real eps)
		{
			// TODO - FIX - mora biti first der == 0!
			Real left_sec_der = Derivation::NSecDer4(_f, x - 4 * eps);
			Real right_sec_der = Derivation::NSecDer4(_f, x + 4 * eps);

			return std::abs(left_sec_der - right_sec_der) < eps;
		}
		bool isContinuous(Real x1, Real x2, int numPoints)
		{
			// kroz sve tocke
			// vidjeti gdje se derivacije razlikuju u znaku (potencijalni min/max, ili discontinuity)
			// istraziti dalje oko te tocke
			return false;
		}
		bool isMonotonic(Real x1, Real x2, int numPoints)
		{
			Real step = (x2 - x1) / numPoints;
			Real prev = _f(x1 + step);
			if (_f(x1) < _f(x1 + step))
			{
				for (int i = 2; i < numPoints; i++)
				{
					Real curr = _f(x1 + i * step);
					if (curr < prev)
						return false;
					prev = curr;
				}
			}
			else
			{
				for (int i = 2; i < numPoints; i++)
				{
					Real curr = _f(x1 + i * step);
					if (curr > prev)
						return false;
					prev = curr;
				}
			}
			return true;
		}
		Real MinInNPoints(Real x1, Real x2, int numPoints)
		{
			Real step = (x2 - x1) / numPoints;
			Real min = _f(x1);

			for (int i = 1; i <= numPoints; i++)
			{
				Real curr = _f(x1 + i * step);
				if (curr < min)
					min = curr;
			}
			return min;
		}
		Real MaxInNPoints(Real x1, Real x2, int numPoints)
		{
			Real step = (x2 - x1) / numPoints;
			Real max = _f(x1);

			for (int i = 1; i < numPoints; i++)
			{
				Real curr = _f(x1 + i * step);
				if (curr > max)
					max = curr;
			}
			return max;
		}
	};

	class RealFunctionComparer
	{
		IRealFunction& _f1;
		IRealFunction& _f2;

	public:
		RealFunctionComparer(IRealFunction& f1, IRealFunction& f2) : _f1(f1), _f2(f2) {}

		Real getAbsDiffSum(Real a, Real b, int numPoints)
		{
			Real step = (b - a) / numPoints;
			Real sum = 0.0;

			for (int i = 0; i < numPoints; i++)
				sum += std::abs(_f1(a + i * step) - _f2(a + i * step));

			return sum;
		}
		Real getAbsDiffAvg(Real a, Real b, int numPoints)
		{
			return getAbsDiffSum(a, b, numPoints) / numPoints;
		}
		Real getAbsDiffMax(Real a, Real b, int numPoints)
		{
			Real step = (b - a) / numPoints;
			Real max = std::abs(_f1(a) - _f2(a));

			for (int i = 0; i < numPoints; i++)
			{
				Real diff = std::abs(_f1(a + i * step) - _f2(a + i * step));
				if (diff > max)
					max = diff;
			}
			return max;
		}
		Real getRelDiffSum(Real a, Real b, int numPoints)
		{
			Real step = (b - a) / numPoints;
			Real sum = 0.0;

			for (int i = 0; i < numPoints; i++)
				if (_f1(a + i * step) != 0.0)
					sum += std::abs(_f1(a + i * step) - _f2(a + i * step)) / std::abs(_f1(a + i * step));
				else
					--numPoints;

			return sum;
		}
		Real getRelDiffAvg(Real a, Real b, int numPoints)
		{
			return getRelDiffSum(a, b, numPoints) / numPoints;
		}
		Real getRelDiffMax(Real a, Real b, int numPoints)
		{
			Real step = (b - a) / numPoints;
			Real max = 0.0;

			for (int i = 0; i < numPoints; i++)
			{
				if (_f1(a + i * step) != 0.0)
				{
					Real diff = std::abs(_f1(a + i * step) - _f2(a + i * step)) / std::abs(_f1(a + i * step));
					if (diff > max)
						max = diff;
				}
			}
			return max;
		}

		///////////                  Integration measures                /////////
		Real getIntegratedDiff(Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			return getIntegratedDiff(_f1, _f2, a, b, method);
		}
		static Real getIntegratedDiff(IRealFunction& f1, IRealFunction& f2, Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			RealFuncDiffHelper helper(f1, f2);

			switch (method)
			{
			//case IntegrationMethod::SIMPSON:
			//	return IntegrateSimpson(helper, a, b);
			//case IntegrationMethod::ROMBERG:
			//	return IntegrateRomberg(helper, a, b);
			default:
				return IntegrateTrap(helper, a, b);
			}
		}

		Real getIntegratedAbsDiff(Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			return getIntegratedAbsDiff(_f1, _f2, a, b, method);
		}
		static Real getIntegratedAbsDiff(IRealFunction& f1, IRealFunction& f2, Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			RealFuncDiffAbsHelper helper(f1, f2);

			switch (method)
			{
			//case IntegrationMethod::SIMPSON:
			//	return IntegrateSimpson(helper, a, b);
			//case IntegrationMethod::ROMBERG:
			//	return IntegrateRomberg(helper, a, b);
			default:
				return IntegrateTrap(helper, a, b);
			}
		}

		Real getIntegratedSqrDiff(Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			return getIntegratedSqrDiff(_f1, _f2, a, b, method);
		}
		static Real getIntegratedSqrDiff(IRealFunction& f1, IRealFunction& f2, Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			RealFuncDiffSqrHelper helper(f1, f2);

			switch (method)
			{
			//case IntegrationMethod::SIMPSON:
			//	return IntegrateSimpson(helper, a, b);
			//case IntegrationMethod::ROMBERG:
			//	return IntegrateRomberg(helper, a, b);
			default:
				return IntegrateTrap(helper, a, b);
			}
		}
	};

	class ScalarFieldAnalyzer
	{
		IScalarFunction<3>& _f;
	public:
		ScalarFieldAnalyzer(IScalarFunction<3>& f) : _f(f) {}

		// kao out parametar vraca "mjeru" konzervativnosti
		bool IsConservative()
		{
			return false;
		}
	};

	class VectorFieldAnalyzer
	{
		IVectorFunction<3>& _f;
	public:
		VectorFieldAnalyzer(IVectorFunction<3>& f) : _f(f) {}

		bool IsSolenoidal()
		{
			return false;
		}
	};
}
#endif