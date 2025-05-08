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

#include <random>

#include "MMLExceptions.h"

// defines for different platforms
#ifdef __APPLE__
#  include <TargetConditionals.h>
#include <random>
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

// defines for different compilers
#if defined(_MSC_VER)
#  define MML_COMPILER_MSVC
#elif defined(__clang__)
#  define MML_COMPILER_CLANG
#elif defined(__GNUC__)
#  define MML_COMPILER_GCC
#elif defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__ECC)
#  define MML_COMPILER_INTEL
#elif defined(__IBMC__) || defined(__IBMCPP__)
#  define MML_COMPILER_IBM
#endif


// Complex must have the same underlaying type as Real
typedef double             Real;      // default real type (float, long double tested, __float128 TODO!!!)
typedef std::complex<Real> Complex;   // default complex type


namespace MML
{
	// Global paths for Visualizers
	static const std::string MML_GLOBAL_PATH = "E:/Projects/ExploringPhysicsWithCpp";
	//static const std::string MML_GLOBAL_PATH = "/usr/zvanjak/projects/ExploringPhysicsWithCpp";

	static const std::string MML_PATH_ResultFiles = MML_GLOBAL_PATH + "/results/";
	
	static const std::string MML_PATH_RealFuncViz = MML_GLOBAL_PATH + 
							"/tools/visualizers/real_function_visualizer/MML_RealFunctionVisualizer.exe";
	
	static const std::string MML_PATH_SurfaceViz = MML_GLOBAL_PATH + 
							"/tools/visualizers/scalar_function_2d_visualizer/MML_ScalarFunction2D_Visualizer.exe";
	
	static const std::string MML_PATH_ParametricCurve2DViz = MML_GLOBAL_PATH + 
							"/tools/visualizers/parametric_curve_2d_visualizer/MML_ParametricCurve2D_Visualizer.exe";
	static const std::string MML_PATH_ParametricCurve3DViz = MML_GLOBAL_PATH + 
							"/tools/visualizers/parametric_curve_3d_visualizer/MML_ParametricCurve3D_Visualizer.exe";
	
	static const std::string MML_PATH_VectorField2DViz = MML_GLOBAL_PATH +
							"/tools/visualizers/vector_field_2d_visualizer/MML_VectorField2D_Visualizer.exe";
	static const std::string MML_PATH_VectorField3DViz = MML_GLOBAL_PATH +
							"/tools/visualizers/vector_field_3d_visualizer/MML_VectorField3D_Visualizer.exe";

	static const std::string MML_PATH_Particle2DViz = MML_GLOBAL_PATH +
							"/tools/visualizers/particle_2d_visualizer/MML_ParticleVisualizer2D.exe";
	static const std::string MML_PATH_Particle3DViz = MML_GLOBAL_PATH +
							"/tools/visualizers/particle_3d_visualizer/MML_ParticleVisualizer3D.exe";

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

	template<class T> inline T POW2(const T& a) { const T& t = a; return t * t; }
	template<class T> inline T POW3(const T& a) { const T& t = a; return t * t * t; }
	template<class T> inline T POW4(const T& a) { const T& t = a; return t * t * t * t; }

	class Random
	{
	public:
		static Real UniformReal(Real min, Real max)
		{
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<Real> dis(min, max);
			return dis(gen);
		}
		static int UniformInt(int min, int max)
		{
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<int> dis(min, max);
			return dis(gen);
		}
	};
	////////////                  Constants                ////////////////
	namespace Constants
	{
		static inline const Real PI = std::numbers::pi;
		static inline const Real INV_PI = std::numbers::inv_pi;
		static inline const Real INV_SQRTPI = std::numbers::inv_sqrtpi;
		
		static inline const Real E = std::numbers::e;
		static inline const Real LN2 = std::numbers::ln2;
		static inline const Real LN10 = std::numbers::ln10;
		
		static inline const Real SQRT2 = std::numbers::sqrt2;
		static inline const Real SQRT3 = std::numbers::sqrt3;

		static inline const Real Eps = std::numeric_limits<Real>::epsilon();
		static inline const Real PosInf = std::numeric_limits<Real>::max();
		static inline const Real NegInf = -std::numeric_limits<Real>::max();
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
	static void SolveQuadratic(const Complex& a, const Complex& b, const Complex& c,
														 Complex& x1, Complex& x2)
	{
		Complex D = b * b - 4.0 * a * c;
		Complex sqrtD = std::sqrt(D);
		x1 = (-b + sqrtD) / (2.0 * a);
		x2 = (-b - sqrtD) / (2.0 * a);
	}
	// Solving cubic equation a * x^3 + b * x^2 + c * x + d = 0
	static int  SolveCubic(Real a, Real b, Real c, Real d,
												Complex& x1, Complex& x2, Complex& x3)
	{
		// Normalize the coefficients
		Real A = b / a;
		Real B = c / a;
		Real C = d / a;

		// Calculate the discriminant
		Real Q = (3.0 * B - POW2(A)) / 9.0;
		Real R = (9.0 * A * B - 27.0 * C - 2.0 * POW3(A)) / 54.0;
		Real D = POW3(Q) + POW2(R); // Discriminant

		if (D >= 0) // Complex or duplicate roots
		{
			Real S = std::cbrt(R + std::sqrt(D));
			Real T = std::cbrt(R - std::sqrt(D));

			x1 = -A / 3.0 + (S + T); // Real root
			x2 = -A / 3.0 - (S + T) / 2.0 + Complex(0, std::sqrt(3.0) * (S - T) / 2.0); // Complex root
			x3 = -A / 3.0 - (S + T) / 2.0 - Complex(0, std::sqrt(3.0) * (S - T) / 2.0); // Complex root

			return 1;
		}
		else // Three real roots
		{
			Real theta = std::acos(R / std::sqrt(-POW3(Q)));
			x1 = 2.0 * std::sqrt(-Q) * std::cos(theta / 3.0) - A / 3.0;
			x2 = 2.0 * std::sqrt(-Q) * std::cos((theta + 2.0 * Constants::PI) / 3.0) - A / 3.0;
			x3 = 2.0 * std::sqrt(-Q) * std::cos((theta + 4.0 * Constants::PI) / 3.0) - A / 3.0;

			return 3;
		}
	}
	// FIX - Solving quartic equation a * x^4 + b * x^3 + c * x^2 + d * x + e = 0
	static void SolveQuartic(Real a, Real b, Real c, Real d, Real e,
													 Complex& x1, Complex& x2, Complex& x3, Complex& x4)
	{
		// Normalize the coefficients
		Real A = b / a;
		Real B = c / a;
		Real C = d / a;
		Real D = e / a;

		// Depressed quartic substitution
		Real AA = A * A;
		Real P = -3.0 / 8.0 * AA + B;
		Real Q = 1.0 / 8.0 * AA * A - 1.0 / 2.0 * A * B + C;
		Real R = -3.0 / 256.0 * AA * AA + 1.0 / 16.0 * AA * B - 1.0 / 4.0 * A * C + D;

		if (std::abs(R) < Constants::Eps)
		{
			// Solve the biquadratic equation
			SolveQuadratic(1.0, P, Q, x1, x2);
			x3 = -x1;
			x4 = -x2;
		}
		else
		{
			// Solve the resolvent cubic
			Complex z1, z2, z3;
			SolveCubic(1.0, -P / 2.0, -R, R * P / 2.0 - Q * Q / 8.0, z1, z2, z3);

			Complex z = z1;
			if (z.imag() != 0.0)
			{
				z = z2;
				if (z.imag() != 0.0)
				{
					z = z3;
				}
			}

			Complex u = std::sqrt(z);
			Complex v = std::sqrt(z - P);

			// Solve the two quadratic equations
			SolveQuadratic(1.0, u, z - Q / (2.0 * u), x1, x2);
			SolveQuadratic(1.0, -u, z + Q / (2.0 * u), x3, x4);
		}

		// Adjust the roots
		x1 -= A / 4.0;
		x2 -= A / 4.0;
		x3 -= A / 4.0;
		x4 -= A / 4.0;
	}


	// Template struct for precision values
	template<typename T>
	struct PrecisionValues;

	// Specialization for float
	template<>
	struct PrecisionValues<float>
	{
		static constexpr float ComplexAreEqualTolerance = 1e-6f;
		static constexpr float ComplexAreEqualAbsTolerance = 1e-6f;

		static constexpr float MatrixIsEqualTolerance = 1e-6f;
		static constexpr float VectorIsEqualTolerance = 1e-6f;

		static constexpr float Pnt2CartIsEqualTolerance = 1e-6f;
		static constexpr float Pnt2PolarIsEqualTolerance = 1e-6f;
		static constexpr float Pnt3CartIsEqualTolerance = 1e-6f;
		static constexpr float Pnt3SphIsEqualTolerance = 1e-6f;
		static constexpr float Pnt3CylIsEqualTolerance = 1e-6f;

		static constexpr float Vec2CartIsEqualTolerance = 1e-6f;
		static constexpr float Vec3CartIsEqualTolerance = 1e-6f;
		static constexpr float Vec3CartIsParallelTolerance = 1e-6f;

		static constexpr float Line3DIsPerpendicularTolerance = 1e-6f;
		static constexpr float Line3DIsParallelTolerance = 1e-6f;
		static constexpr float Plane3DIsPointOnPlaneTolerance = 1e-6f;

		static constexpr float IsMatrixSymmetricTolerance = 1e-6f;
		static constexpr float IsMatrixDiagonalTolerance = 1e-6f;
		static constexpr float IsMatrixUnitTolerance = 1e-6f;
		static constexpr float IsMatrixOrthogonalTolerance = 1e-6f;
	};

	// Specialization for double
	template<>
	struct PrecisionValues<double>
	{
		static constexpr double ComplexAreEqualTolerance = 1e-10;
		static constexpr double ComplexAreEqualAbsTolerance = 1e-10;

		static constexpr double MatrixIsEqualTolerance = 1e-10;
		static constexpr double VectorIsEqualTolerance = 1e-10;

		static constexpr double Pnt2CartIsEqualTolerance = 1e-10;
		static constexpr double Pnt2PolarIsEqualTolerance = 1e-10;
		static constexpr double Pnt3CartIsEqualTolerance = 1e-10;
		static constexpr double Pnt3SphIsEqualTolerance = 1e-10;
		static constexpr double Pnt3CylIsEqualTolerance = 1e-10;

		static constexpr double Vec2CartIsEqualTolerance = 1e-10;
		static constexpr double Vec3CartIsEqualTolerance = 1e-10;
		static constexpr double Vec3CartIsParallelTolerance = 1e-10;

		static constexpr double Line3DIsPerpendicularTolerance = 1e-10;
		static constexpr double Line3DIsParallelTolerance = 1e-10;
		static constexpr double Plane3DIsPointOnPlaneTolerance = 1e-10;

		static constexpr double IsMatrixSymmetricTolerance = 1e-10;
		static constexpr double IsMatrixDiagonalTolerance = 1e-10;
		static constexpr double IsMatrixUnitTolerance = 1e-10;
		static constexpr double IsMatrixOrthogonalTolerance = 1e-10;
	};

	// Specialization for long double
	template<>
	struct PrecisionValues<long double>
	{
		static constexpr long double ComplexAreEqualTolerance = 1e-15L;
		static constexpr long double ComplexAreEqualAbsTolerance = 1e-15L;

		static constexpr long double MatrixIsEqualTolerance = 1e-15L;
		static constexpr long double VectorIsEqualTolerance = 1e-15L;

		static constexpr long double Pnt2CartIsEqualTolerance = 1e-15;
		static constexpr long double Pnt2PolarIsEqualTolerance = 1e-15;
		static constexpr long double Pnt3CartIsEqualTolerance = 1e-15;
		static constexpr long double Pnt3SphIsEqualTolerance = 1e-15;
		static constexpr long double Pnt3CylIsEqualTolerance = 1e-15;

		static constexpr long double Vec2CartIsEqualTolerance = 1e-15;
		static constexpr long double Vec3CartIsEqualTolerance = 1e-15;
		static constexpr long double Vec3CartIsParallelTolerance = 1e-15;

		static constexpr long double Line3DIsPerpendicularTolerance = 1e-15;
		static constexpr long double Line3DIsParallelTolerance = 1e-15;
		static constexpr long double Plane3DIsPointOnPlaneTolerance = 1e-15;

		static constexpr long double IsMatrixSymmetricTolerance = 1e-15;
		static constexpr long double IsMatrixDiagonalTolerance = 1e-15;
		static constexpr long double IsMatrixUnitTolerance = 1e-15;
		static constexpr long double IsMatrixOrthogonalTolerance = 1e-15;
	};

	namespace Defaults
	{
		// Output defaults
		static int VectorPrintWidth = 15;
		static int VectorPrintPrecision = 10;
		static int VectorNPrintWidth = 15;
		static int VectorNPrintPrecision = 10;

		//////////               Default precisions             ///////////
		// Use the precision values based on the Real type
		static inline const Real ComplexAreEqualTolerance    = PrecisionValues<Real>::ComplexAreEqualTolerance;
		static inline const Real ComplexAreEqualAbsTolerance = PrecisionValues<Real>::ComplexAreEqualAbsTolerance;
		static inline const Real VectorIsEqualTolerance			 = PrecisionValues<Real>::VectorIsEqualTolerance;
		static inline const Real MatrixIsEqualTolerance			 = PrecisionValues<Real>::MatrixIsEqualTolerance;

		static inline const Real Pnt2CartIsEqualTolerance  = PrecisionValues<Real>::Pnt2CartIsEqualTolerance;
		static inline const Real Pnt2PolarIsEqualTolerance = PrecisionValues<Real>::Pnt2PolarIsEqualTolerance;
		static inline const Real Pnt3CartIsEqualTolerance  = PrecisionValues<Real>::Pnt3CartIsEqualTolerance;
		static inline const Real Pnt3SphIsEqualTolerance	 = PrecisionValues<Real>::Pnt3SphIsEqualTolerance;
		static inline const Real Pnt3CylIsEqualTolerance   = PrecisionValues<Real>::Pnt3CylIsEqualTolerance;

		static inline const Real Vec2CartIsEqualTolerance		 = PrecisionValues<Real>::Vec2CartIsEqualTolerance;
		static inline const Real Vec3CartIsEqualTolerance		 = PrecisionValues<Real>::Vec3CartIsEqualTolerance;
		static inline const Real Vec3CartIsParallelTolerance = PrecisionValues<Real>::Vec3CartIsParallelTolerance;

		static inline const Real Line3DIsPerpendicularTolerance = PrecisionValues<Real>::Line3DIsPerpendicularTolerance;
		static inline const Real Line3DIsParallelTolerance			= PrecisionValues<Real>::Line3DIsParallelTolerance;
		static inline const Real Plane3DIsPointOnPlaneTolerance = PrecisionValues<Real>::Plane3DIsPointOnPlaneTolerance;

		static inline const Real IsMatrixSymmetricTolerance  = PrecisionValues<Real>::IsMatrixSymmetricTolerance;
		static inline const Real IsMatrixDiagonalTolerance	 = PrecisionValues<Real>::IsMatrixDiagonalTolerance;
		static inline const Real IsMatrixUnitTolerance			 = PrecisionValues<Real>::IsMatrixUnitTolerance;
		static inline const Real IsMatrixOrthogonalTolerance = PrecisionValues<Real>::IsMatrixOrthogonalTolerance;


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