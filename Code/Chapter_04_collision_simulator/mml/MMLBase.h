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
static const std::string GLOB_PATH_ResultFiles = "E:\\Projects\\MinimalMathLibrary\\results\\";
static const std::string GLOB_PATH_RealFuncViz = "E:\\Projects\\MinimalMathLibrary\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe";
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
		return sqrt(a.real() * a.real() + a.imag() * a.imag());
	}

	template<class T> inline T POW2(const T &a) { return ((a) * (a)); }
	template<class T> inline T POW3(const T &a) { return ((a) * (a) * (a)); }
	template<class T> inline T POW4(const T &a) { return ((a) * (a) * (a) * (a)); }
	template<class T> inline T POW5(const T &a) { return ((a) * (a) * (a) * (a) * (a)); }
	template<class T> inline T POW6(const T &a) { return ((a) * (a) * (a) * (a) * (a) * (a)); }

	template<class T>
	inline T SIGN(const T& a, const T& b)
	{
		return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
	}

	////////////                  Constants                ////////////////
	namespace Constants
	{
		static inline const Real PI = std::numbers::pi;
		static inline const Real Epsilon = std::numeric_limits<Real>::epsilon();
		static inline const Real PositiveInf = std::numeric_limits<Real>::max();
		static inline const Real NegativeInf = -std::numeric_limits<Real>::max();
	}

	namespace Defaults
	{
    static int VectorPrintWidth = 15;
    static int VectorPrintPrecision = 10;

		//////////               Default precisions             ///////////
		// TODO - make dependent on Real type (ie. different values for float, double and long double)
		static inline const double ComplexEqualityPrecision = 1e-15;
		static inline const double ComplexAbsEqualityPrecision = 1e-15;
		static inline const double MatrixEqualityPrecision = 1e-15;
		static inline const double VectorEqualityPrecision = 1e-15;

		static inline const double IsMatrixSymmetricPrecision = 1e-15;
		static inline const double IsMatrixDiagonalPrecision = 1e-15;
		static inline const double IsMatrixUnitPrecision = 1e-15;
		static inline const double IsMatrixOrthogonalPrecision = 1e-15;

		static inline const double DerivationDefaultStep = 1e-6;

		static inline const int    IntegrateTrapMaxSteps = 20;
		static inline const double IntegrateTrapEPS = 1.0e-5;

		static inline const int    IntegrateSimpMaxSteps = 20;
		static inline const double IntegrateSimpEPS = 1.0e-5;

		static inline const int    IntegrateRombMaxSteps = 20;
		static inline const double IntegrateRombEPS = 1.0e-6;

		static inline const double WorkIntegralPrecision = 1e-05;
		static inline const double LineIntegralPrecision = 1e-05;
	}
}
#endif