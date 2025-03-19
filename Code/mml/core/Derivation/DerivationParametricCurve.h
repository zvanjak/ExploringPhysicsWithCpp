#if !defined MML_DERIVATION_PARAMETRIC_CURVE_H
#define MML_DERIVATION_PARAMETRIC_CURVE_H

#include "MMLBase.h"

#include "DerivationBase.h"

#include "base/VectorN.h"
#include "base/MatrixNM.h"

namespace MML
{
	namespace Derivation
	{
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer1(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> y0 = f(t);
			VectorN<Real, N> diff = yh - y0;

			if (error)
			{
				VectorN<Real, N> ym = f(t - h);
				VectorN<Real, N> ypph_vec = yh - 2 * y0 + ym;

				Real ypph = ypph_vec.NormL2() / h;

				*error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * Constants::Eps / h;
			}
			return diff / h;
		}

		template <int N>
		static VectorN<Real, N> NDer1(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer1(f, t, NDer1_h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer2(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> diff = yh - ymh;

			if (error)
			{
				VectorN<Real, N> yth = f(t + 2 * h);
				VectorN<Real, N> ymth = f(t - 2 * h);

				*error = Constants::Eps * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
			}
			return diff / (2 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer2(f, t, NDer2_h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer4(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> y2h = f(t + 2 * h);
			VectorN<Real, N> ym2h = f(t - 2 * h);

			VectorN<Real, N> y2 = ym2h - y2h;
			VectorN<Real, N> y1 = yh - ymh;

			if (error)
			{
				VectorN<Real, N> y3h = f(t + 3 * h);
				VectorN<Real, N> ym3h = f(t - 3 * h);

				*error = std::abs((y3h - ym3h).NormL2() / 2 + 2 * (ym2h - y2h).NormL2() + 5 * (yh - ymh).NormL2() / 2) / (30 * h);
				*error += Constants::Eps * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer4(f, t, NDer4_h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SIXTH order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer6(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = f(t - 2 * h) - f(t + 2 * h);
			VectorN<Real, N> y3 = f(t + 3 * h) - f(t - 3 * h);

			if (error)
			{
				VectorN<Real, N> y7 = (f(t + 4 * h) - f(t - 4 * h) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2()) * Constants::Eps / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer6(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer6(f, t, NDer6_h, error);
		}	
		/********************************************************************************************************************/
		/********                               Numerical derivatives of EIGHTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer8(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = f(t - 2 * h) - f(t + 2 * h);
			VectorN<Real, N> y3 = f(t + 3 * h) - f(t - 3 * h);
			VectorN<Real, N> y4 = f(t - 4 * h) - f(t + 4 * h);

			VectorN<Real, N> tmp1 = 3 * y4 / 8 + 4 * y3;
			VectorN<Real, N> tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				VectorN<Real, N> f9 = (f(t + 5 * h) - f(t - 5 * h)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2()) * Constants::Eps / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer8(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer8(f, t, NDer8_h, error);
		}

		/********************************************************************************************************************/
		/********                                      SECOND DERIVATIVES                                            ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NSecDer1(const IParametricCurve<N>& f, Real x, Real h, Real* error = nullptr)
		{
			VectorN<Real, N>  yh = NDer2(f, x + h, h, error);
			VectorN<Real, N>  y0 = NDer2(f, x, h, error);
			VectorN<Real, N>  diff = yh - y0;
			if (error)
			{
				VectorN<Real, N> ym = NDer2(f, x - h, h, error);
				VectorN<Real, N> ypph_vec = (yh - 2 * y0 + ym) / h;

				Real ypph = ypph_vec.NormL2();

				*error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * Constants::Eps / h;
			}
			return diff / h;
		}
		template <int N>
		static VectorN<Real, N> NSecDer1(const IParametricCurve<N>& f, Real x, Real* error = nullptr)
		{
			return NSecDer1(f, x, NDer1_h, error);
		}

		template <int N>
		static VectorN<Real, N> NSecDer2(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NDer4(f, t + h, error);
			VectorN<Real, N> ymh = NDer4(f, t - h, error);
			VectorN<Real, N> diff = yh - ymh;

			if (error)
			{
				VectorN<Real, N> yth = NDer4(f, t + 2 * h, error);
				VectorN<Real, N> ymth = NDer4(f, t - 2 * h, error);

				*error = Constants::Eps * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
			}
			return diff / (2 * h);
		}
		template <int N>
		static VectorN<Real, N> NSecDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer2(f, t, NDer2_h, error);
		}

		template <int N>
		static VectorN<Real, N> NSecDer4(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NDer6(f, t + h, error);
			VectorN<Real, N> ymh = NDer6(f, t - h, error);
			VectorN<Real, N> y2h = NDer6(f, t + 2 * h, error);
			VectorN<Real, N> ym2h = NDer6(f, t - 2 * h, error);

			VectorN<Real, N> y2 = ym2h - y2h;
			VectorN<Real, N> y1 = yh - ymh;

			if (error)
			{
				VectorN<Real, N> y3h = NDer6(f, t + 3 * h, error);
				VectorN<Real, N> ym3h = NDer6(f, t - 3 * h, error);

				*error = std::abs((y3h - ym3h).NormL2() / 2 + 2 * (ym2h - y2h).NormL2() + 5 * (yh - ymh).NormL2() / 2) / (30 * h);
				*error += Constants::Eps * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}
		template <int N>
		static VectorN<Real, N> NSecDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer4(f, t, NDer4_h, error);
		}

		template <int N>
		static VectorN<Real, N> NSecDer6(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NDer8(f, t + h, error);
			VectorN<Real, N> ymh = NDer8(f, t - h, error);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = NDer8(f, t - 2 * h, error) - NDer8(f, t + 2 * h, error);
			VectorN<Real, N> y3 = NDer8(f, t + 3 * h, error) - NDer8(f, t - 3 * h, error);

			if (error)
			{
				VectorN<Real, N> y7 = (NDer8(f, t + 4 * h, error) - NDer8(f, t - 4 * h, error) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2()) * Constants::Eps / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}
		template <int N>
		static VectorN<Real, N> NSecDer6(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer6(f, t, NDer6_h, error);
		}

		template <int N>
		static VectorN<Real, N> NSecDer8(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NDer8(f, t + h, error);
			VectorN<Real, N> ymh = NDer8(f, t - h, error);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = NDer8(f, t - 2 * h, error) - NDer8(f, t + 2 * h, error);
			VectorN<Real, N> y3 = NDer8(f, t + 3 * h, error) - NDer8(f, t - 3 * h, error);
			VectorN<Real, N> y4 = NDer8(f, t - 4 * h, error) - NDer8(f, t + 4 * h, error);

			VectorN<Real, N> tmp1 = 3 * y4 / 8 + 4 * y3;
			VectorN<Real, N> tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				VectorN<Real, N> f9 = (NDer8(f, t + 5 * h, error) - NDer8(f, t - 5 * h, error)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2()) * Constants::Eps / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}
		template <int N>
		static VectorN<Real, N> NSecDer8(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer8(f, t, NDer8_h, error);
		}

		/********************************************************************************************************************/
		/********                                       THIRD DERIVATIVES                                            ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NThirdDer1(const IParametricCurve<N>& f, Real x, Real h, Real* error = nullptr)
		{
			VectorN<Real, N>  yh = NSecDer2(f, x + h, h, error);
			VectorN<Real, N>  y0 = NSecDer2(f, x, h, error);
			VectorN<Real, N>  diff = yh - y0;
			if (error)
			{
				VectorN<Real, N> ym = NSecDer2(f, x - h, h, error);
				VectorN<Real, N> ypph_vec = (yh - 2 * y0 + ym) / h;

				Real ypph = ypph_vec.NormL2();

				*error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * Constants::Eps / h;
			}
			return diff / h;
		}
		template <int N>
		static VectorN<Real, N> NThirdDer1(const IParametricCurve<N>& f, Real x, Real* error = nullptr)
		{
			return NThirdDer1(f, x, NDer1_h, error);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer2(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NSecDer4(f, t + h, error);
			VectorN<Real, N> ymh = NSecDer4(f, t - h, error);
			VectorN<Real, N> diff = yh - ymh;

			if (error)
			{
				VectorN<Real, N> yth = NSecDer4(f, t + 2 * h, error);
				VectorN<Real, N> ymth = NSecDer4(f, t - 2 * h, error);

				*error = Constants::Eps * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
			}
			return diff / (2 * h);
		}
		template <int N>
		static VectorN<Real, N> NThirdDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer2(f, t, NDer2_h, error);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer4(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NSecDer6(f, t + h, error);
			VectorN<Real, N> ymh = NSecDer6(f, t - h, error);
			VectorN<Real, N> y2h = NSecDer6(f, t + 2 * h, error);
			VectorN<Real, N> ym2h = NSecDer6(f, t - 2 * h, error);

			VectorN<Real, N> y2 = ym2h - y2h;
			VectorN<Real, N> y1 = yh - ymh;

			if (error)
			{
				VectorN<Real, N> y3h = NSecDer6(f, t + 3 * h, error);
				VectorN<Real, N> ym3h = NSecDer6(f, t - 3 * h, error);

				*error = std::abs((y3h - ym3h).NormL2() / 2 + 2 * (ym2h - y2h).NormL2() + 5 * (yh - ymh).NormL2() / 2) / (30 * h);
				*error += Constants::Eps * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}
		template <int N>
		static VectorN<Real, N> NThirdDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer4(f, t, NDer4_h, error);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer6(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NSecDer8(f, t + h, error);
			VectorN<Real, N> ymh = NSecDer8(f, t - h, error);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = NSecDer8(f, t - 2 * h, error) - NSecDer8(f, t + 2 * h, error);
			VectorN<Real, N> y3 = NSecDer8(f, t + 3 * h, error) - NSecDer8(f, t - 3 * h, error);

			if (error)
			{
				VectorN<Real, N> y7 = (NSecDer8(f, t + 4 * h, error) - NSecDer8(f, t - 4 * h, error) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2()) * Constants::Eps / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}
		template <int N>
		static VectorN<Real, N> NThirdDer6(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer6(f, t, NDer6_h, error);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer8(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NSecDer8(f, t + h, error);
			VectorN<Real, N> ymh = NSecDer8(f, t - h, error);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = NSecDer8(f, t - 2 * h, error) - NSecDer8(f, t + 2 * h, error);
			VectorN<Real, N> y3 = NSecDer8(f, t + 3 * h, error) - NSecDer8(f, t - 3 * h, error);
			VectorN<Real, N> y4 = NSecDer8(f, t - 4 * h, error) - NSecDer8(f, t + 4 * h, error);

			VectorN<Real, N> tmp1 = 3 * y4 / 8 + 4 * y3;
			VectorN<Real, N> tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				VectorN<Real, N> f9 = (NSecDer8(f, t + 5 * h, error) - NSecDer8(f, t - 5 * h, error)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2()) * Constants::Eps / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}
		template <int N>
		static VectorN<Real, N> NThirdDer8(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer8(f, t, NDer8_h, error);
		}
		
		/********************************************************************************************************************/
		/********                            Definitions of default derivation functions                             ********/
		/********************************************************************************************************************/
		template<int N>
		static inline VectorN<Real, N>(*DeriveCurve)(const IParametricCurve<N>& f, 
																								 Real x, Real* error) = Derivation::NDer4;
		template<int N>
		static inline VectorN<Real, N>(*DeriveCurveSec)(const IParametricCurve<N>& f, 
																										Real x, Real* error) = Derivation::NSecDer4;
		template<int N>
		static inline VectorN<Real, N>(*DeriveCurveThird)(const IParametricCurve<N>& f, 
																											Real x, Real* error) = Derivation::NThirdDer4;

	}
}

#endif // MML_DERIVATION_PARAMETRIC_CURVE_H