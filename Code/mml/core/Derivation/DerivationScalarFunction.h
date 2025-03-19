#if !defined MML_DERIVATION_SCALAR_FUNCTION_H
#define MML_DERIVATION_SCALAR_FUNCTION_H

#include "MMLBase.h"

#include "DerivationBase.h"

#include "base/VectorN.h"

namespace MML
{
	namespace Derivation
	{
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer1Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, 
														 Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x = point;
			Real y0 = f(x);

			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = orig_x - h;
				Real ym = f(x);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Eps / h;
			}
			return diff / h;
		}
		template <int N>
		static Real NDer1Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, 
														 Real* error = nullptr)
		{
			return NDer1Partial(f, deriv_index, point, NDer1_h, error);
		}
		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, 
																							Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer1Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer1Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, 
																							VectorN<Real, N>* error = nullptr)
		{
			return NDer1PartialByAll(f, point, NDer1_h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer2Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			auto    x = point;
			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			x[deriv_index] = orig_x - h;
			Real ymh = f(x);

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f(x);

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f(x);

				*error = Constants::Eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}
		template <int N>
		static Real NDer2Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, deriv_index, point, NDer2_h, error);
		}
		
		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer2Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer2Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer2PartialByAll(f, point, NDer2_h, error);
		}
		
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer4Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, 
														 Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			x[deriv_index] = orig_x - h;
			Real ymh = f(x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f(x);

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f(x);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Eps * (std::abs(y2h) + std::abs(ym2h) + 
																				8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}
		template <int N>
		static Real NDer4Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer4Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer4Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer4PartialByAll(f, point, NDer4_h, error);
		}
		
		/********************************************************************************************************************/
		/********                               Numerical derivatives of SIXTH order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer6Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };

			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			x[deriv_index] = orig_x - h;
			Real ymh = f(x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x);

			x[deriv_index] = orig_x + 3 * h;
			Real y3h = f(x);

			x[deriv_index] = orig_x - 3 * h;
			Real ym3h = f(x);

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;

			if (error)
			{
				x[deriv_index] = orig_x + 4 * h;
				Real y4h = f(x);

				x[deriv_index] = orig_x - 4 * h;
				Real ym4h = f(x);

				Real y7 = (y4h - ym4h - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Eps / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}
		template <int N>
		static Real NDer6Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer6Partial(f, deriv_index, point, NDer6_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer6Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer6Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer6PartialByAll(f, point, NDer6_h, error);
		}
		
		/********************************************************************************************************************/
		/********                               Numerical derivatives of EIGHTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer8Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };

			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			x[deriv_index] = orig_x - h;
			Real ymh = f(x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x);

			x[deriv_index] = orig_x + 3 * h;
			Real y3h = f(x);

			x[deriv_index] = orig_x - 3 * h;
			Real ym3h = f(x);

			x[deriv_index] = orig_x + 4 * h;
			Real y4h = f(x);

			x[deriv_index] = orig_x - 4 * h;
			Real ym4h = f(x);

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;
			Real y4 = ym4h - y4h;

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				x[deriv_index] = orig_x + 5 * h;
				Real y5h = f(x);

				x[deriv_index] = orig_x - 5 * h;
				Real ym5h = f(x);

				Real f9 = (y5h - ym5h) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Eps / h;
			}

			return (tmp1 + tmp2) / (105 * h);
		}
		template <int N>
		static Real NDer8Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer8Partial(f, deriv_index, point, NDer8_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer8Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer8Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer8PartialByAll(f, point, NDer8_h, error);
		}
		
		/********************************************************************************************************************/
		/********                                      SECOND DERIVATIVES                                            ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NSecDer1Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
																const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real x_orig_val = point[der_ind2];

			auto x_eval_pos = point;
			Real y0 = NDer2Partial(f, der_ind1, x_eval_pos, error);
			x_eval_pos[der_ind2] = x_orig_val + h;
			Real yh = NDer2Partial(f, der_ind1, x_eval_pos, error);

			Real diff = yh - y0;
			if (error)
			{
				x_eval_pos[der_ind2] = x_orig_val - h;

				Real ym = NDer2Partial(f, der_ind1, x_eval_pos, error);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;

				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Eps / h;
			}
			return diff / h;
		}
		template <int N>
		static Real NSecDer1Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
																const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer1Partial(f, der_ind1, der_ind2, point, NDer1_h, error);
		}

		template <int N>
		static Real NSecDer2Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer2Partial(f, der_ind1, der_ind2, point, NDer2_h, error);
		}
		template <int N>
		static Real NSecDer2Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real orig_x = point[der_ind2];
			auto x_eval_pos = point;

			x_eval_pos[der_ind2] = orig_x + h;
			Real yh = NDer4Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - h;
			Real ymh = NDer4Partial(f, der_ind1, x_eval_pos, error);

			Real diff = yh - ymh;

			if (error)
			{
				x_eval_pos[der_ind2] = orig_x + 2 * h;
				Real y2h = NDer4Partial(f, der_ind1, x_eval_pos, error);

				x_eval_pos[der_ind2] = orig_x - 2 * h;
				Real ym2h = NDer4Partial(f, der_ind1, x_eval_pos, error);

				*error = Constants::Eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}
		template <int N>

		static Real NSecDer4Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[der_ind2];
			auto x_eval_pos = point;

			x_eval_pos[der_ind2] = orig_x + h;
			Real yh = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - h;
			Real ymh = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x + 2 * h;
			Real y2h = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - 2 * h;
			Real ym2h = NDer6Partial(f, der_ind1, x_eval_pos, error);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x_eval_pos[der_ind2] = orig_x + 3 * h;
				Real y3h = NDer6Partial(f, der_ind1, x_eval_pos, error);

				x_eval_pos[der_ind2] = orig_x - 3 * h;
				Real ym3h = NDer6Partial(f, der_ind1, x_eval_pos, error);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}
		template <int N>
		static Real NSecDer4Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer4Partial(f, der_ind1, der_ind2, point, NDer4_h, error);
		}

		template <int N>
		static Real NSecDer6Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer6Partial(f, der_ind1, der_ind2, point, NDer6_h, error);
		}
		template <int N>
		static Real NSecDer6Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[der_ind2];
			auto x_eval_pos = point;

			x_eval_pos[der_ind2] = orig_x + h;
			Real yh = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - h;
			Real ymh = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x + 2 * h;
			Real y2h = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - 2 * h;
			Real ym2h = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x + 3 * h;
			Real y3h = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - 3 * h;
			Real ym3h = NDer6Partial(f, der_ind1, x_eval_pos, error);

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;

			if (error)
			{
				x_eval_pos[der_ind2] = orig_x + 4 * h;
				Real y4h = NDer6Partial(f, der_ind1, x_eval_pos, error);

				x_eval_pos[der_ind2] = orig_x - 4 * h;
				Real ym4h = NDer6Partial(f, der_ind1, x_eval_pos, error);

				Real y7 = (y4h - ym4h - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Eps / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		template <int N>
		static Real NSecDer8Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer8Partial(f, der_ind1, der_ind2, point, NDer8_h, error);
		}
		template <int N>
		static Real NSecDer8Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[der_ind2];
			auto x_eval_pos = point;

			x_eval_pos[der_ind2] = orig_x + h;
			Real yh = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - h;
			Real ymh = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x + 2 * h;
			Real y2h = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - 2 * h;
			Real ym2h = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x + 3 * h;
			Real y3h = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - 3 * h;
			Real ym3h = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x + 4 * h;
			Real y4h = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - 4 * h;
			Real ym4h = NDer8Partial(f, der_ind1, x_eval_pos, error);

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;
			Real y4 = ym4h - y4h;

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				x_eval_pos[der_ind2] = orig_x + 5 * h;
				Real y5h = NDer8Partial(f, der_ind1, x_eval_pos, error);

				x_eval_pos[der_ind2] = orig_x - 5 * h;
				Real ym5h = NDer8Partial(f, der_ind1, x_eval_pos, error);

				Real f9 = (y5h - ym5h) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Eps / h;
			}

			return (tmp1 + tmp2) / (105 * h);
		}

		/********************************************************************************************************************/
		/********                            Definitions of default derivation functions                             ********/
		/********************************************************************************************************************/
		template<int N>
		static inline Real(*DerivePartial)(const IScalarFunction<N>& f, int deriv_index, 
																			 const VectorN<Real, N>& point, Real* error) = Derivation::NDer4Partial;
		template<int N>
		static inline Real(*DeriveSecPartial)(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
																					const VectorN<Real, N>& point, Real* error) = Derivation::NSecDer4Partial;
		template<int N>
		static inline VectorN<Real, N>(*DerivePartialAll)(const IScalarFunction<N>& f, const VectorN<Real, N>& point, 
																											VectorN<Real, N>* error) = Derivation::NDer4PartialByAll;

	}
}

#endif // MML_DERIVATION_SCALAR_FUNCTION_H
