#if !defined MML_DERIVATION_H
#define MML_DERIVATION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"
#include "interfaces/ITensorField.h"

#include "base/VectorN.h"
#include "base/MatrixNM.h"

namespace MML
{
	// TODO - make this a NAMESPACE, instead of class
	// - all second and third derivation calculations should be moved to the end of the file
	class Derivation
	{
	public:
		static inline const Real NDer1_h = 2 * std::sqrt(Constants::Epsilon);
		static inline const Real NDer2_h = std::pow(3 * Constants::Epsilon, 1.0 / 3.0);
		static inline const Real NDer4_h = std::pow(11.25 * Constants::Epsilon, 1.0 / 5.0);     // 0.0012009323661373839 for double!
		static inline const Real NDer6_h = std::pow(Constants::Epsilon / 168.0, 1.0 / 7.0);
		static inline const Real NDer8_h = std::pow(551.25 * Constants::Epsilon, 1.0 / 9.0);

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/
		
		//////////////////////////              RealFunction            //////////////////////////
		static Real NDer1(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real y0 = f(x);
			Real diff = yh - y0;
			if (error)
			{
				Real ym = f(x - h);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;

				// h*|f''(x)|*0.5 + (|f(x+h)+|f(x)|) * eps/h
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}
		static Real NDer1(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^1/2
			// Note that this estimate of h differs from the best estimate by a factor of sqrt((|f(x)| + |f(x+h)|)/|f''(x)|).
			// Since this factor is invariant under the scaling f -> kf, then we are somewhat justified in approximating it by 1.
			// This approximation will get better as we move to higher orders of accuracy.
			return NDer1(f, x, NDer1_h, error);
		}
		static Real NDer1(const IRealFunction& f, Real x)
		{
			return NDer1(f, x, NDer1_h, nullptr);
		}
		
		static Real NDer1Left(const IRealFunction& f, Real x, Real* error = nullptr) 
		{ return NDer1(f, x - 2 * NDer1_h, NDer1_h, error); }
		static Real NDer1Right(const IRealFunction& f, Real x, Real* error = nullptr) 
		{ return NDer1(f, x + 2 * NDer1_h, NDer1_h, error); }
		static Real NDer1Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) 
		{ return NDer1(f, x - 2 * h, h, error); }
		static Real NDer1Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) 
		{ return NDer1(f, x + 2 * h, h, error); }

		static Real NSecDer1(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NDer2(f, x + h, h, error);
			Real y0 = NDer2(f, x, h, error);
			Real diff = yh - y0;
			if (error)
			{
				Real ym = NDer2(f, x - h, h, error);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;

				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}
		static Real NSecDer1(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer1(f, x, NDer1_h, error);
		}

		static Real NThirdDer1(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NSecDer2(f, x + h, h, error);
			Real y0 = NSecDer2(f, x, h, error);
			Real diff = yh - y0;
			if (error)
			{
				Real ym = NSecDer2(f, x - h, h, error);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;

				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}
		static Real NThirdDer1(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer1(f, x, NDer1_h, error);
		}

		//////////////////////////             ScalarFunction           //////////////////////////
		template <int N>
		static Real NDer1Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, 
														 Real* error = nullptr)
		{
			return NDer1Partial(f, deriv_index, point, NDer1_h, error);
		}

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
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
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

				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, 
																							VectorN<Real, N>* error = nullptr)
		{
			return NDer1PartialByAll(f, point, NDer1_h, error);
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

		//////////////////////////             VectorFunction           //////////////////////////
		template <int N>
		static Real NDer1Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, func_index, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static Real NDer1Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto x = point;

			Real x_orig = x[deriv_index];
			Real y0 = f(x)[func_index];

			x[deriv_index] = x_orig + h;
			Real yh = f(x)[func_index];

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = x_orig - h;
				Real ym = f(x)[func_index];
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer1PartialByAll(f, func_index, point, NDer1_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer1Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer1Partial(f, func_index, i, point, h);
			}

			return ret;
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer1PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer1PartialAllByAll(f, point, NDer1_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer1PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer1Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer1Partial(f, i, j, point, h);
				}

			return ret;
		}

		//////////////////////////             TensorField           /////////////////////////////
		template <int N>
		static Real NDer1Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, i, j, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto x = point;

			Real x_orig = x[deriv_index];
			Real y0 = f.Component(i, j, x);

			x[deriv_index] = x_orig + h;
			Real yh = f.Component(i, j, x);

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = x_orig - h;
				Real ym = f.Component(i, j, x);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static Real NDer1Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, i, j, k, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto x = point;

			Real x_orig = x[deriv_index];
			Real y0 = f.Component(i, j, k, x);

			x[deriv_index] = x_orig + h;
			Real yh = f.Component(i, j, k, x);

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = x_orig - h;
				Real ym = f.Component(i, j, k, x);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static Real NDer1Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, i, j, k, l, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto x = point;

			Real x_orig = x[deriv_index];
			Real y0 = f.Component(i, j, k, l, x);

			x[deriv_index] = x_orig + h;
			Real yh = f.Component(i, j, k, l, x);

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = x_orig - h;
				Real ym = f.Component(i, j, k, l, x);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		/////////////////////////             ParametricCurve           //////////////////////////
		template <int N>
		static VectorN<Real, N> NDer1(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer1(f, t, NDer1_h, error);
		}

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

				*error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static VectorN<Real, N> NSecDer1(const IParametricCurve<N>& f, Real x, Real* error = nullptr)
		{
			return NSecDer1(f, x, NDer1_h, error);
		}

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

				*error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static VectorN<Real, N> NThirdDer1(const IParametricCurve<N>& f, Real x, Real* error = nullptr)
		{
			return NThirdDer1(f, x, NDer1_h, error);
		}

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

				*error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/

		//////////////////////////              RealFunction            //////////////////////////
		static Real NDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real diff = yh - ymh;
			if (error)
			{
				Real y2h = f(x + 2 * h);
				Real ym2h = f(x - 2 * h);
				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}
		static Real NDer2(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^2/3
			// See the previous discussion to understand determination of h and the error bound.
			// Series[(f[x+h] - f[x-h])/(2*h), {h, 0, 4}]

			return NDer2(f, x, NDer2_h, error);
		}
		static Real NDer2(const IRealFunction& f, Real x)
		{
			return NDer2(f, x, NDer2_h, nullptr);
		}
		
		static Real NDer2Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer2(f, x - 2 * NDer2_h, NDer2_h, error); }
		static Real NDer2Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer2(f, x + 2 * NDer2_h, NDer2_h, error); }
		static Real NDer2Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer2(f, x - 3 * h, h, error); }
		static Real NDer2Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer2(f, x + 3 * h, h, error); }

		static Real NSecDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NDer4(f, x + h, error);
			Real ymh = NDer4(f, x - h, error);
			Real diff = yh - ymh;
			if (error)
			{
				Real y2h = NDer4(f, x + 2 * h, error);
				Real ym2h = NDer4(f, x - 2 * h, error);

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}
		static Real NSecDer2(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer2(f, x, NDer2_h, error);
		}

		static Real NThirdDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NSecDer4(f, x + h, error);
			Real ymh = NSecDer4(f, x - h, error);
			Real diff = yh - ymh;
			if (error)
			{
				Real y2h = NSecDer4(f, x + 2 * h, error);
				Real ym2h = NSecDer4(f, x - 2 * h, error);

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}
		static Real NThirdDer2(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer2(f, x, NDer2_h, error);
		}

		//////////////////////////             ScalarFunction           //////////////////////////
		template <int N>
		static Real NDer2Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, deriv_index, point, NDer2_h, error);
		}

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

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
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

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer2PartialByAll(f, point, NDer2_h, error);
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

		//////////////////////////             VectorFunction           //////////////////////////
		template <int N>
		static Real NDer2Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, func_index, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static Real NDer2Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f(x)[func_index];

			x[deriv_index] = orig_x - h;
			Real ymh = f(x)[func_index];

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f(x)[func_index];

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f(x)[func_index];

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer2PartialByAll(f, func_index, point, NDer2_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer2Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer2Partial(f, func_index, i, point, h);
			}

			return ret;
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer2PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer2PartialAllByAll(f, point, NDer2_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer2PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer2Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer2Partial(f, i, j, point, h);
				}

			return ret;
		}

		//////////////////////////             TensorField           //////////////////////////
		template <int N>
		static Real NDer2Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, i, j, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, x);

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f.Component(i, j, x);

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f.Component(i, j, x);

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, i, j, k, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, k, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, k, x);

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f.Component(i, j, k, x);

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f.Component(i, j, k, x);

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, i, j, k, l, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, k, l, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, k, l, x);

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f.Component(i, j, k, l, x);

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f.Component(i, j, k, l, x);

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}
		/////////////////////////             ParametricCurve           /////////////////////////
		template <int N>
		static VectorN<Real, N> NDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer2(f, t, NDer2_h, error);
		}

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

				*error = Constants::Epsilon * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
			}
			return diff / (2 * h);
		}

		template <int N>
		static VectorN<Real, N> NSecDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer2(f, t, NDer2_h, error);
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

				*error = Constants::Epsilon * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
			}
			return diff / (2 * h);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer2(f, t, NDer2_h, error);
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

				*error = Constants::Epsilon * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
			}
			return diff / (2 * h);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/

		//////////////////////////              RealFunction            //////////////////////////
		static Real NDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real y2h = f(x + 2 * h);
			Real ym2h = f(x - 2 * h);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				// Mathematica code to extract the remainder:
				// Series[(f[x-2*h]+ 8*f[x+h] - 8*f[x-h] - f[x+2*h])/(12*h), {h, 0, 7}]
				Real y3h = f(x + 3 * h);
				Real ym3h = f(x - 3 * h);

				// Error from fifth derivative:
				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				// Error from function evaluation:
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}
		static Real NDer4(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^4/5
			return NDer4(f, x, NDer4_h, error);
		}
		static Real NDer4(const IRealFunction& f, Real x)
		{
			return NDer4(f, x, NDer4_h, nullptr);
		}

		static Real NDer4Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer4(f, x - 4 * NDer4_h, NDer4_h, error); }
		static Real NDer4Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer4(f, x + 4 * NDer4_h, NDer4_h, error); }
		static Real NDer4Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer4(f, x - 4 * h, h, error); }
		static Real NDer4Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer4(f, x + 4 * h, h, error); }

		static Real NSecDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NDer6(f, x + h, error);
			Real ymh = NDer6(f, x - h, error);
			Real y2h = NDer6(f, x + 2 * h, error);
			Real ym2h = NDer6(f, x - 2 * h, error);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				Real y3h = NDer6(f, x + 3 * h, error);
				Real ym3h = NDer6(f, x - 3 * h, error);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}
		static Real NSecDer4(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer4(f, x, NDer4_h, error);
		}

		static Real NThirdDer4(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer4(f, x, NDer4_h, error);
		}
		static Real NThirdDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NSecDer6(f, x + h, error);
			Real ymh = NSecDer6(f, x - h, error);
			Real y2h = NSecDer6(f, x + 2 * h, error);
			Real ym2h = NSecDer6(f, x - 2 * h, error);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				Real y3h = NSecDer6(f, x + 3 * h, error);
				Real ym3h = NSecDer6(f, x - 3 * h, error);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		//////////////////////////             ScalarFunction           //////////////////////////
		template <int N>
		static Real NDer4Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static Real NDer4Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
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
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		template <int N>
		static Real NSecDer4Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer4Partial(f, der_ind1, der_ind2, point, NDer4_h, error);
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
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer4PartialByAll(f, point, NDer4_h, error);
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

		//////////////////////////             VectorFunction           //////////////////////////
		template <int N>
		static Real NDer4Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, 
														 Real* error = nullptr)
		{
			return NDer4Partial(f, func_index, deriv_index, point, NDer4_h, error);
		}
		template <int N>
		static Real NDer4Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, 
														 Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f(x)[func_index];

			x[deriv_index] = orig_x - h;
			Real ymh = f(x)[func_index];

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x)[func_index];

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x)[func_index];

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f(x)[func_index];

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f(x)[func_index];

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}
		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, 
																							VectorN<Real, N>* error = nullptr)
		{
			return NDer4PartialByAll(f, func_index, point, NDer4_h, error);
		}
		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, 
																							VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer4Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer4Partial(f, func_index, i, point, h);
			}

			return ret;
		}
		template <int N>
		static MatrixNM<Real, N, N> NDer4PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, 
																										 MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer4PartialAllByAll(f, point, NDer4_h, error);
		}
		template <int N>
		static MatrixNM<Real, N, N> NDer4PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, 
																										 MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer4Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer4Partial(f, i, j, point, h);
				}

			return ret;
		}

		//////////////////////////             TensorField           //////////////////////////
		template <int N>
		static Real NDer4Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, i, j, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f.Component(i, j, x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f.Component(i, j, x);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f.Component(i, j, x);

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f.Component(i, j, x);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, i, j, k, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, k, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, k, x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f.Component(i, j, k, x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f.Component(i, j, k, x);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f.Component(i, j, k, x);

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f.Component(i, j, k, x);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}


		template <int N>
		static Real NDer4Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, i, j, k, l, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, k, l, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, k, l, x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f.Component(i, j, k, l, x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f.Component(i, j, k, l, x);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f.Component(i, j, k, l, x);

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f.Component(i, j, k, l, x);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		/////////////////////////             ParametricCurve           /////////////////////////
		template <int N>
		static VectorN<Real, N> NDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer4(f, t, NDer4_h, error);
		}

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
				*error += Constants::Epsilon * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		template <int N>
		static VectorN<Real, N> NSecDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer4(f, t, NDer4_h, error);
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
				*error += Constants::Epsilon * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}


		template <int N>
		static VectorN<Real, N> NThirdDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer4(f, t, NDer4_h, error);
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
				*error += Constants::Epsilon * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SIXTH order                                 ********/
		/********************************************************************************************************************/

		//////////////////////////              RealFunction            //////////////////////////
		static Real NDer6(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			const Real eps = (std::numeric_limits<Real>::epsilon)();

			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real y1 = yh - ymh;
			Real y2 = f(x - 2 * h) - f(x + 2 * h);
			Real y3 = f(x + 3 * h) - f(x - 3 * h);

			if (error)
			{
				// Mathematica code to generate fd scheme for 7th derivative:
				// Sum[(-1)^i*Binomial[7, i]*(f[x+(3-i)*h] + f[x+(4-i)*h])/2, {i, 0, 7}]
				// Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
				// Series[(f[x+4*h]-f[x-4*h] + 6*(f[x-3*h] - f[x+3*h]) + 14*(f[x-h] - f[x+h] + f[x+2*h] - f[x-2*h]))/2, {h, 0, 15}]
				Real y7 = (f(x + 4 * h) - f(x - 4 * h) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}
		static Real NDer6(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^6/7
			// Error: h^6f^(7)(x)/140 + 5|f(x)|eps/h
			return NDer6(f, x, NDer6_h, error);
		}
		static Real NDer6(const IRealFunction& f, Real x)
		{
			// Error bound ~eps^6/7
			// Error: h^6f^(7)(x)/140 + 5|f(x)|eps/h
			return NDer6(f, x, NDer6_h, nullptr);
		}

		static Real NDer6Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer6(f, x - 5 * NDer6_h, NDer6_h, error); }
		static Real NDer6Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer6(f, x + 5 * NDer6_h, NDer6_h, error); }
		static Real NDer6Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer6(f, x - 5 * h, h, error); }
		static Real NDer6Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer6(f, x + 5 * h, h, error); }

		static Real NSecDer6(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NDer8(f, x + h, error);
			Real ymh = NDer8(f, x - h, error);
			Real y1 = yh - ymh;
			Real y2 = NDer8(f, x - 2 * h, error) - NDer8(f, x + 2 * h, error);
			Real y3 = NDer8(f, x + 3 * h, error) - NDer8(f, x - 3 * h, error);

			if (error)
			{
				Real y7 = (NDer8(f, x + 4 * h, error) - NDer8(f, x - 4 * h, error) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}
		static Real NSecDer6(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer6(f, x, NDer6_h, error);
		}

		static Real NThirdDer6(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer6(f, x, NDer6_h, error);
		}
		static Real NThirdDer6(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NSecDer8(f, x + h, error);
			Real ymh = NSecDer8(f, x - h, error);
			Real y1 = yh - ymh;
			Real y2 = NSecDer8(f, x - 2 * h, error) - NSecDer8(f, x + 2 * h, error);
			Real y3 = NSecDer8(f, x + 3 * h, error) - NSecDer8(f, x - 3 * h, error);

			if (error)
			{
				Real y7 = (NSecDer8(f, x + 4 * h, error) - NSecDer8(f, x - 4 * h, error) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		//////////////////////////             ScalarFunction           //////////////////////////
		template <int N>
		static Real NDer6Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer6Partial(f, deriv_index, point, NDer6_h, error);
		}

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

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
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

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer6PartialByAll(f, point, NDer6_h, error);
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

		//////////////////////////             VectorFunction           //////////////////////////
		template <int N>
		static Real NDer6Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer6Partial(f, func_index, deriv_index, point, NDer6_h, error);
		}

		template <int N>
		static Real NDer6Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };

			x[deriv_index] = orig_x + h;
			Real yh = f(x)[func_index];

			x[deriv_index] = orig_x - h;
			Real ymh = f(x)[func_index];

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x)[func_index];

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x)[func_index];

			x[deriv_index] = orig_x + 3 * h;
			Real y3h = f(x)[func_index];

			x[deriv_index] = orig_x - 3 * h;
			Real ym3h = f(x)[func_index];

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;

			if (error)
			{
				x[deriv_index] = orig_x + 4 * h;
				Real y4h = f(x)[func_index];

				x[deriv_index] = orig_x - 4 * h;
				Real ym4h = f(x)[func_index];

				Real y7 = (y4h - ym4h - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer6PartialByAll(f, func_index, point, NDer6_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer6Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer6Partial(f, func_index, i, point, h);
			}

			return ret;
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer6PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer6PartialAllByAll(f, point, NDer6_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer6PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer6Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer6Partial(f, i, j, point, h);
				}

			return ret;
		}

		/////////////////////////             ParametricCurve           /////////////////////////
		template <int N>
		static VectorN<Real, N> NDer6(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer6(f, t, NDer6_h, error);
		}

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

				*error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2()) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		template <int N>
		static VectorN<Real, N> NSecDer6(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer6(f, t, NDer6_h, error);
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

				*error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2()) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer6(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer6(f, t, NDer6_h, error);
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

				*error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2()) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of EIGHTH order                                ********/
		/********************************************************************************************************************/

		//////////////////////////              RealFunction            //////////////////////////
		static Real NDer8(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real y1 = yh - ymh;
			Real y2 = f(x - 2 * h) - f(x + 2 * h);
			Real y3 = f(x + 3 * h) - f(x - 3 * h);
			Real y4 = f(x - 4 * h) - f(x + 4 * h);

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				// Mathematica code to generate fd scheme for 7th derivative:
				// Sum[(-1)^i*Binomial[9, i]*(f[x+(4-i)*h] + f[x+(5-i)*h])/2, {i, 0, 9}]
				// Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
				// Series[(f[x+5*h]-f[x- 5*h])/2 + 4*(f[x-4*h] - f[x+4*h]) + 27*(f[x+3*h] - f[x-3*h])/2 + 24*(f[x-2*h]  - f[x+2*h]) + 21*(f[x+h] - f[x-h]), {h, 0, 15}]
				Real f9 = (f(x + 5 * h) - f(x - 5 * h)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}
		static Real NDer8(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^8/9.
			// In Real precision, we only expect to lose two digits of precision while using this formula, at the cost of 8 function evaluations.
			// Error: h^8|f^(9)(x)|/630 + 7|f(x)|eps/h assuming 7 unstabilized additions.
			// Mathematica code to get the error:
			// Series[(f[x+h]-f[x-h])*(4/5) + (1/5)*(f[x-2*h] - f[x+2*h]) + (4/105)*(f[x+3*h] - f[x-3*h]) + (1/280)*(f[x-4*h] - f[x+4*h]), {h, 0, 9}]
			// If we used Kahan summation, we could get the max error down to h^8|f^(9)(x)|/630 + |f(x)|eps/h.

			return NDer8(f, x, NDer8_h, error);
		}
		static Real NDer8(const IRealFunction& f, Real x)
		{
			return NDer8(f, x, NDer8_h, nullptr);
		}

		static Real NDer8Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer8(f, x - 6 * NDer8_h, NDer8_h, error); }
		static Real NDer8Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer8(f, x + 6 * NDer8_h, NDer8_h, error); }
		static Real NDer8Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer8(f, x - 6 * h, h, error); }
		static Real NDer8Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer8(f, x + 6 * h, h, error); }

		static Real NSecDer8(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer8(f, x, NDer8_h, error);
		}
		static Real NSecDer8(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NDer8(f, x + h, error);
			Real ymh = NDer8(f, x - h, error);
			Real y1 = yh - ymh;
			Real y2 = NDer8(f, x - 2 * h, error) - NDer8(f, x + 2 * h, error);
			Real y3 = NDer8(f, x + 3 * h, error) - NDer8(f, x - 3 * h, error);
			Real y4 = NDer8(f, x - 4 * h, error) - NDer8(f, x + 4 * h, error);

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				Real f9 = (NDer8(f, x + 5 * h, error) - NDer8(f, x - 5 * h, error)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}

		static Real NThirdDer8(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer8(f, x, NDer8_h, error);
		}
		static Real NThirdDer8(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NSecDer8(f, x + h, error);
			Real ymh = NSecDer8(f, x - h, error);
			Real y1 = yh - ymh;
			Real y2 = NSecDer8(f, x - 2 * h, error) - NSecDer8(f, x + 2 * h, error);
			Real y3 = NSecDer8(f, x + 3 * h, error) - NSecDer8(f, x - 3 * h, error);
			Real y4 = NSecDer8(f, x - 4 * h, error) - NSecDer8(f, x + 4 * h, error);

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				Real f9 = (NSecDer8(f, x + 5 * h, error) - NSecDer8(f, x - 5 * h, error)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}

		//////////////////////////             ScalarFunction           //////////////////////////
		template <int N>
		static Real NDer8Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer8Partial(f, deriv_index, point, NDer8_h, error);
		}

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

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}

			return (tmp1 + tmp2) / (105 * h);
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

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}

			return (tmp1 + tmp2) / (105 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer8PartialByAll(f, point, NDer8_h, error);
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

		//////////////////////////             VectorFunction           //////////////////////////
		template <int N>
		static Real NDer8Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer8Partial(f, func_index, deriv_index, point, NDer8_h, error);
		}

		template <int N>
		static Real NDer8Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };

			x[deriv_index] = orig_x + h;
			Real yh = f(x)[func_index];

			x[deriv_index] = orig_x - h;
			Real ymh = f(x)[func_index];

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x)[func_index];

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x)[func_index];

			x[deriv_index] = orig_x + 3 * h;
			Real y3h = f(x)[func_index];

			x[deriv_index] = orig_x - 3 * h;
			Real ym3h = f(x)[func_index];

			x[deriv_index] = orig_x + 4 * h;
			Real y4h = f(x)[func_index];

			x[deriv_index] = orig_x - 4 * h;
			Real ym4h = f(x)[func_index];

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;
			Real y4 = ym4h - y4h;

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				x[deriv_index] = orig_x + 5 * h;
				Real y5h = f(x)[func_index];

				x[deriv_index] = orig_x - 5 * h;
				Real ym5h = f(x)[func_index];

				Real f9 = (y5h - ym5h) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}

			return (tmp1 + tmp2) / (105 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer8PartialByAll(f, func_index, point, NDer8_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer8Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer8Partial(f, func_index, i, point, h);
			}

			return ret;
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer8PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer8PartialAllByAll(f, point, NDer8_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer8PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer8Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer8Partial(f, i, j, point, h);
				}

			return ret;
		}

		/////////////////////////             ParametricCurve           /////////////////////////
		template <int N>
		static VectorN<Real, N> NDer8(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer8(f, t, NDer8_h, error);
		}

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

				*error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2()) * Constants::Epsilon / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}

		template <int N>
		static VectorN<Real, N> NSecDer8(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer8(f, t, NDer8_h, error);
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

				*error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2()) * Constants::Epsilon / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer8(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer8(f, t, NDer8_h, error);
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

				*error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2()) * Constants::Epsilon / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}

		/********************************************************************************************************************/
		/********************************************************************************************************************/
		/********                                      SECOND DERIVATIVES                                            ********/
		/********************************************************************************************************************/
		/********************************************************************************************************************/

		/********************************************************************************************************************/
		/********************************************************************************************************************/
		/********                                       THIRD DERIVATIVES                                            ********/
		/********************************************************************************************************************/
		/********************************************************************************************************************/



		/********************************************************************************************************************/
		/********                            Definitions of default derivation functions                             ********/
		/********************************************************************************************************************/
		static inline Real(*Derive)(const IRealFunction& f, Real x) = Derivation::NDer4;
		static inline Real(*DeriveErr)(const IRealFunction& f, Real x, Real* error) = Derivation::NDer4;

		static inline Real(*DeriveSec)(const IRealFunction& f, Real x, Real* error) = Derivation::NSecDer4;

		static inline Real(*DeriveThird)(const IRealFunction& f, Real x, Real* error) = Derivation::NThirdDer2;


		template<int N>
		static inline Real(*DerivePartial)(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error) = Derivation::NDer4Partial;

		template<int N>
		static inline Real(*DeriveSecPartial)(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real* error) = Derivation::NSecDer4Partial;

		template<int N>
		static inline VectorN<Real, N>(*DerivePartialAll)(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error) = Derivation::NDer4PartialByAll;


		template<int N>
		static inline Real(*DeriveVecPartial)(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error) = Derivation::NDer4Partial;

		template<int N>
		static inline VectorN<Real, N>(*DeriveVecPartialAll)(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error) = Derivation::NDer4PartialByAll;

		template<int N>
		static inline MatrixNM<Real, N, N>(*DeriveVecPartialAllByAll)(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error) = Derivation::NDer4PartialAllByAll;


		template<int N>
		static inline VectorN<Real, N>(*DeriveCurve)(const IParametricCurve<N>& f, Real x, Real* error) = Derivation::NDer4;

		template<int N>
		static inline VectorN<Real, N>(*DeriveCurveSec)(const IParametricCurve<N>& f, Real x, Real* error) = Derivation::NSecDer4;

		template<int N>
		static inline VectorN<Real, N>(*DeriveCurveThird)(const IParametricCurve<N>& f, Real x, Real* error) = Derivation::NThirdDer4;
	

		/********************************************************************************************************************/
		/********************************************************************************************************************/
		/********                                      CALCULATING JACOBIANS                                         ********/
		/********************************************************************************************************************/
		/********************************************************************************************************************/

  	template<int N>
		static MatrixNM<Real, N, N> calcJacobian(const IVectorFunction<N>& func, const VectorN<Real, N>& pos)
		{
			MatrixNM<Real, N, N> jac;

			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
				{
					jac(i, j) = NDer4Partial(func, i, j, pos);
				}

			return jac;
		}

    template<int N, int M>
		static MatrixNM<Real, M, N> calcJacobian(const IVectorFunctionNM<N, M>& func, const VectorN<Real, N>& pos)
		{
			MatrixNM<Real, M, N> jac;

			for (int i = 0; i < M; ++i)
				for (int j = 0; j < N; ++j)
				{
					jac(i, j) = Derivation::NDer4Partial(func, i, j, pos);
				}

			return jac;
		}    
  };
}

#endif