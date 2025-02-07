#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#if !defined  MML_FUNCTIONS_H
#define MML_FUNCTIONS_H

#include "MMLBase.h"

namespace MML
{
    namespace Functions
    {
        // Functions of REAL domain
        // TODO - umjesto Real, naprosto staviti sve tri dostupne varijante?
        static inline Real Sin(Real x) { return sin(x); }
        static inline Real Cos(Real x) { return cos(x); }
        static inline Real Sec(Real x) { return 1.0 / cos(x); }
        static inline Real Csc(Real x) { return 1.0 / sin(x); }
        static inline Real Tan(Real x) { return tan(x); }
        static inline Real Ctg(Real x) { return 1.0 / tan(x); }
        
        static inline Real Exp(Real x) { return exp(x); }
        static inline Real Log(Real x) { return log(x); }
        static inline Real Log10(Real x){ return log10(x); }
        static inline Real Sqrt(Real x) { return sqrt(x); }
        static inline Real Pow(Real x, Real y) { return pow(x, y); }
        
        static inline Real Sinh(Real x) { return sinh(x); }
        static inline Real Cosh(Real x) { return cosh(x); }
        static inline Real Sech(Real x) { return 1.0 / cosh(x); }
        static inline Real Csch(Real x) { return 1.0 / sinh(x); }        
        static inline Real Tanh(Real x) { return tanh(x); }
        static inline Real Ctgh(Real x) { return 1.0 / tanh(x); }        
        
        static inline Real Asin(Real x) { return asin(x); }
        static inline Real Acos(Real x) { return acos(x); }
        static inline Real Atan(Real x) { return atan(x); }

        static inline Real Asinh(Real x) { return asinh(x); }
        static inline Real Acosh(Real x) { return acosh(x); }
        static inline Real Atanh(Real x) { return atanh(x); }

        static inline Real Erf(Real x)  { return std::erf(x); }
        static inline Real Erfc(Real x) { return std::erfc(x); }

        static inline Real TGamma(Real x) { return std::tgamma(x); }
        static inline Real LGamma(Real x) { return std::lgamma(x); }
        static inline Real RiemannZeta(Real x) { return std::riemann_zeta(x); }
        static inline Real Comp_ellint_1(Real x) { return std::comp_ellint_1(x); }
        static inline Real Comp_ellint_2(Real x) { return std::comp_ellint_2(x); }

        static inline Real Hermite(unsigned int n, Real x) { return std::hermite(n, x); }
        static inline Real Legendre(unsigned int n, Real x) { return std::legendre(n, x); }
        static inline Real Laguerre(unsigned int n, Real x) { return std::laguerre(n, x); }
        static inline Real SphBessel(unsigned int n, Real x) { return std::sph_bessel(n, x); }
        static inline Real SphLegendre(int n1, int n2, Real x) { return std::sph_legendre(n1, n2, x); }
        
        // Functions of COMPLEX domain
        static inline Complex Sin(Complex x) { return sin(x); }
        static inline Complex Cos(Complex x) { return cos(x); }
        static inline Complex Sec(Complex x) { return Real{1.0} / cos(x); }
        static inline Complex Csc(Complex x) { return Real{1.0} / sin(x); }
        static inline Complex Tan(Complex x) { return tan(x); }
        static inline Complex Ctg(Complex x) { return Real{1.0} / tan(x); }
        
        static inline Complex Exp(Complex x) { return exp(x); }
        static inline Complex Log(Complex x) { return log(x); }
        static inline Complex Sqrt(Complex x) { return sqrt(x); }
        static inline Complex Pow(Complex x, Complex y) { return pow(x, y); }
        
        static inline Complex Sinh(Complex x) { return sinh(x); }
        static inline Complex Cosh(Complex x) { return cosh(x); }
        static inline Complex Sech(Complex x) { return Complex(1.0) / cosh(x); }
        static inline Complex Csch(Complex x) { return Complex(1.0) / sinh(x); }        
        static inline Complex Tanh(Complex x) { return tanh(x); }
        static inline Complex Ctgh(Complex x) { return Complex(1.0) / tanh(x); } 

        static inline Complex Asin(Complex x) { return asin(x); }
        static inline Complex Acos(Complex x) { return acos(x); }
        static inline Complex Atan(Complex x) { return atan(x); }

        static inline Complex Asinh(Complex x) { return asinh(x); }
        static inline Complex Acosh(Complex x) { return acosh(x); }
        static inline Complex Atanh(Complex x) { return atanh(x); }    

        static inline Real Factorial(int n) {
            Real fact = 1.0;
            for( int i=2; i<=n; i++)
                fact *= i;
            return fact;
        }   
        static inline long long FactorialInt(int n) {
            long long fact = 1;
            for( int i=2; i<=n; i++)
                fact *= i;
            return fact;
        }           
    }
}

#endif