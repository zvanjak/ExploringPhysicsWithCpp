#if !defined MML_BASE_H
#define MML_BASE_H

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <stdexcept>
#include <initializer_list>
#include <algorithm>
#include <memory>
#include <functional>

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <cmath>
#include <limits>
#include <complex>
#include <numbers>

#include "MMLExceptions.h"

#ifdef __APPLE__
#  include <TargetConditionals.h>
#  if (defined(TARGET_OS_OSX) && TARGET_OS_OSX == 1) || (defined(TARGET_OS_MAC) && TARGET_OS_MAC == 1)
#    define MML_PLATFORM_MAC
#  elif (defined(TARGET_OS_IPHONE) && TARGET_OS_IPHONE == 1)
#    define MML_PLATFORM_IPHONE
#  endif

#elif defined(linux) || defined(__linux) || defined(__linux__)
#  define MML_PLATFORM_LINUX

#elif defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) || defined(__MINGW32__)
#  define MML_PLATFORM_WINDOWS
#endif

// Complex must have the same underlaying type as Real
typedef double               Real;      // default real type
typedef std::complex<double> Complex;   // default complex type


namespace MML
{
	// Global paths for Visualizers
	static const std::string MML_GLOBAL_PATH = "E:\\Projects\\MinimalMathLibrary";

	static const std::string MML_PATH_ResultFiles = MML_GLOBAL_PATH + "\\results\\";
	static const std::string MML_PATH_RealFuncViz = MML_GLOBAL_PATH + "\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe";
	static const std::string MML_PATH_SurfaceViz  = MML_GLOBAL_PATH + "\\tools\\visualizers\\scalar_function_2d_visualizer\\MML_ScalarFunction2Visualizer.exe";
	static const std::string MML_PATH_ParametricCurveViz = MML_GLOBAL_PATH + "\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe";
	static const std::string MML_PATH_VectorFieldViz = MML_GLOBAL_PATH + "\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe";

	template<class Type>
	static Real Abs(const Type& a)
	{
		return std::abs(a);
	}
	template<class Type>
	static Real Abs(const std::complex<Type>& a)
	{
		return hypot(a.real(), a.imag());
	}
	template<class Type>
	static Real Sqrt(const std::complex<Type>& a)
	{
		return std::sqrt(a);
	}
	inline bool isWithinAbsPrec(Real a, Real b, Real eps)
	{
		return std::abs(a - b) < eps;
	}
	inline bool isWithinRelPrec(Real a, Real b, Real eps)
	{
		return std::abs(a - b) < eps * std::max(Abs(a), Abs(b));
	}

	template<class T> inline T POW2(const T a) { const T t=a; return t * t; }
	template<class T> inline T POW3(const T a) { const T t=a; return t * t * t; }
	template<class T> inline T POW4(const T a) { const T t=a; return t * t * t * t; }

	////////////                  Constants                ////////////////
	namespace Constants
	{
		static inline const Real PI = std::numbers::pi;
		static inline const Real E = std::numbers::e;

		static inline const Real Epsilon = std::numeric_limits<Real>::epsilon();
		static inline const Real PositiveInf = std::numeric_limits<Real>::max();
		static inline const Real NegativeInf = -std::numeric_limits<Real>::max();
	}

	// a * x^2 + b * x + c = 0
	static int  SolveQuadratic(Real a, Real b, Real c, 
														 Complex& x1, Complex& x2)
	{
		Real D = b * b - 4 * a * c;
		if (D >= 0)
		{
			Real sqrtD = sqrt(D);
			x1 = (-b + sqrtD) / (2 * a);
			x2 = (-b - sqrtD) / (2 * a);
			return 2;
		}
		else
		{
			Complex sqrtD = std::sqrt(Complex(D));
			x1 = (-b + sqrtD) / (2 * a);
			x2 = (-b - sqrtD) / (2 * a);
			return 0;
		}
	}
	static void SolveQuadratic(const Complex &a, const Complex &b, const Complex &c, 
														 Complex& x1, Complex& x2)
	{
		Complex D = b * b - 4.0 * a * c;
		Complex sqrtD = std::sqrt(D);
		x1 = (-b + sqrtD) / (2.0 * a);
		x2 = (-b - sqrtD) / (2.0 * a);
	}
	// CHECK!!!
	static void SolveCubic(Real a, Real b, Real c, Real d, 
												 Complex& x1, Complex& x2, Complex& x3)
	{
		// based on Cardano's method
		Real p = c / a - b * b / (3 * a * a);
		Real q = d / a + 2 * b * b * b / (27 * a * a * a) - b * c / (3 * a * a);
		Real D = q * q + 4 * p * p * p / 27;
		if (D >= 0)
		{
			Real sqrtD = sqrt(D);
			Real u = cbrt((-q + sqrtD) / 2);
			Real v = cbrt((-q - sqrtD) / 2);
			x1 = u + v - b / (3 * a);
			x2 = -0.5 * (u + v) - b / (3 * a) + 0.5 * sqrt(3) * (u - v) * Complex(0, 1);
			x3 = -0.5 * (u + v) - b / (3 * a) - 0.5 * sqrt(3) * (u - v) * Complex(0, 1);
		}
		else
		{
			Real sqrtD = std::sqrt(-D);
			Complex u = cbrt((-q + sqrtD) / 2.0);
			Complex v = cbrt((-q - sqrtD) / 2.0);
			x1 = u + v - b / (3 * a);
			x2 = -0.5 * (u + v) - b / (3 * a) + 0.5 * sqrt(3) * (u - v) * Complex(0, 1);
			x3 = -0.5 * (u + v) - b / (3 * a) - 0.5 * sqrt(3) * (u - v) * Complex(0, 1);
		}
	}
	// TODO!!!
	static void SolveQuartic(Real a, Real b, Real c, Real d, Real e, 
													 Complex& x1, Complex& x2, Complex& x3, Complex& x4)
	{
		// based on Ferrari's method
		Real p = -3 * b * b / (8 * a * a) + c / a;
		Real q = b * b * b / (8 * a * a * a) - b * c / (2 * a * a) + d / a;
		Real r = -3 * b * b * b * b / (256 * a * a * a * a) + b * b * c / (16 * a * a * a) - b * d / (4 * a * a) + e / a;
		Real D0 = q * q - 4 * p * p * p;
		Real D1 = 2 * q * q * q - 9 * p * p * q + 27 * p * p * p * r;
		Real C = cbrt((D1 + sqrt(D1 * D1 - 4 * D0 * D0 * D0)) / 2);
		Real u = 0;
		if (C == 0)
			u = 0;
		else
			u = cbrt((D1 - sqrt(D1 * D1 - 4 * D0 * D0 * D0)) / 2);
		Real y1 = -b / (4 * a) + (C + u + p / C) / 2;
		Real y2 = -b / (4 * a) + (-C - u + p / -C) / 2;
		Real y3 = -b / (4 * a) + (C - u + p / C) / 2;
		Real y4 = -b / (4 * a) + (-C + u + p / -C) / 2;
		//SolveQuadratic(1, y1, 0, 1, x1, x2);
		//SolveQuadratic(1, y2, 0, 1, x3, x4);
	}
	
	template<class T>
	inline T SIGN(const T& a, const T& b)
	{
		return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
	}

	namespace Defaults
	{
    static int VectorPrintWidth = 15;
    static int VectorPrintPrecision = 10;

		//////////               Default precisions             ///////////
		// TODO - how to make dependent on Real type 
		// //			(ie. different values for float, double and long double)
		static inline const double ComplexEqualityPrecision = 1e-15;
		static inline const double ComplexAbsEqualityPrecision = 1e-15;
		static inline const double MatrixEqualityPrecision = 1e-15;
		static inline const double VectorEqualityPrecision = 1e-15;

		static inline const double IsMatrixSymmetricPrecision = 1e-15;
		static inline const double IsMatrixDiagonalPrecision = 1e-15;
		static inline const double IsMatrixUnitPrecision = 1e-15;
		static inline const double IsMatrixOrthogonalPrecision = 1e-15;

		static inline const double DerivationDefaultStep = 1e-6;

		static inline const double TrapezoidIntegrationEPS = 1.0e-4;
		static inline const double SimpsonIntegrationEPS = 1.0e-5;
		static inline const double RombergIntegrationEPS = 1.0e-6;

		static inline const double WorkIntegralPrecision = 1e-05;
		static inline const double LineIntegralPrecision = 1e-05;

		//////////             Algorithm constants             ///////////	
		static inline const int    BisectionMaxSteps = 50;
		static inline const int    NewtonRaphsonMaxSteps = 20;

		static inline const int    TrapezoidIntegrationMaxSteps = 20;
		static inline const int    SimpsonIntegrationMaxSteps = 20;
		static inline const int    RombergIntegrationMaxSteps = 20;
		static inline const int    RombergIntegrationUsedPnts = 5;

		static inline const int    ODESolverMaxSteps = 50000;

	}
}
#endif