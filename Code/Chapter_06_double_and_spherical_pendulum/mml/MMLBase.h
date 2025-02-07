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

// Global paths for Visualizers
static const std::string GLOB_PATH_ResultFiles = "C:\\Projects\\MinimalMathLibrary\\results\\";
static const std::string GLOB_PATH_RealFuncViz = "C:\\Projects\\MinimalMathLibrary\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe";
static const std::string GLOB_PATH_SurfaceViz = "E:\\Projects\\MinimalMathLibrary\\tools\\visualizers\\scalar_function_2d_visualizer\\MML_ScalarFunction2Visualizer.exe";
static const std::string GLOB_PATH_ParametricCurveViz = "E:\\Projects\\MinimalMathLibrary\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe";
static const std::string GLOB_PATH_VectorFieldViz = "E:\\Projects\\MinimalMathLibrary\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe";

namespace MML
{
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

	static int  SolveQuadratic(Real a, Real b, Real c, Complex& x1, Complex& x2)
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
	static void SolveQuadratic(const Complex &a, const Complex &b, const Complex &c, Complex& x1, Complex& x2)
	{
		Complex D = b * b - 4.0 * a * c;
		Complex sqrtD = std::sqrt(D);
		x1 = (-b + sqrtD) / (2.0 * a);
		x2 = (-b - sqrtD) / (2.0 * a);
	}
	static void SolveCubic(Real a, Real b, Real c, Real d, Real& x1, Real& x2, Real& x3)
	{
		// based on Cardano's method
		Real p = b - a * a / 3;
		Real q = c + 2 * a * a * a / 27 - a * b / 3;
		Real D = q * q / 4 + p * p * p / 27;
		if (D > 0)
		{
			Real sqrtD = sqrt(D);
			Real u = cbrt(-q / 2 + sqrtD);
			Real v = cbrt(-q / 2 - sqrtD);
			x1 = u + v - a / 3;
			x2 = -0.5 * (u + v) - a / 3;
			x3 = -0.5 * (u + v) - a / 3;
		}
		else if (D == 0)
		{
			Real u = cbrt(-q / 2);
			x1 = 2 * u - a / 3;
			x2 = -u - a / 3;
			x3 = -u - a / 3;
		}
		else
		{
			Real phi = acos(-q / 2 / sqrt(-p * p * p / 27));
			Real u = 2 * sqrt(-p / 3) * cos(phi / 3);
			Real v = 2 * sqrt(-p / 3) * cos((phi + 2 * Constants::PI) / 3);
			Real w = 2 * sqrt(-p / 3) * cos((phi + 4 * Constants::PI) / 3);
			x1 = u - a / 3;
			x2 = v - a / 3;
			x3 = w - a / 3;
		}
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

		static inline const int    TrapezoidIntegrationMaxSteps = 20;
		static inline const double TrapezoidIntegrationEPS = 1.0e-4;
		static inline const int    SimpsonIntegrationMaxSteps = 20;
		static inline const double SimpsonIntegrationEPS = 1.0e-5;
		static inline const int    RombergIntegrationMaxSteps = 20;
		static inline const int    RombergIntegrationUsedPnts = 5;
		static inline const double RombergIntegrationEPS = 1.0e-6;

		static inline const double WorkIntegralPrecision = 1e-05;
		static inline const double LineIntegralPrecision = 1e-05;
	}
}
#endif