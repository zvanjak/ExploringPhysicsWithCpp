#if !defined  MML_METRIC_TENSOR_H
#define MML_METRIC_TENSOR_H

#include "interfaces/IFunction.h"
#include "interfaces/ICoordTransf.h"

#include "base/VectorN.h"
#include "base/MatrixNM.h"

#include "core/Derivation.h"

namespace MML
{
	template<int N>
	class MetricTensorField : public ITensorField2<N>
	{
	public:
		MetricTensorField() : ITensorField2<N>(2, 0) { }
		MetricTensorField(int numContra, int numCo) : ITensorField2<N>(numContra, numCo) { }

		// implementing operator() required by IFunction interface
		Tensor2<N>   operator()(const VectorN<Real, N>& pos) const
		{
			Tensor2<N> ret(this->getNumContravar(), this->getNumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					ret(i, j) = this->Component(i, j, pos);

			return ret;
		}

		Real GetChristoffelSymbolFirstKind(int i, int j, int k, const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N>& g = *this;

			Real gamma_ijk = 0.0;
			for (int m = 0; m < N; m++)
			{
				gamma_ijk += g.Component(m, k, pos) * GetChristoffelSymbolSecondKind(m, i, j, pos);
			}
			return gamma_ijk;
		}
		Real GetChristoffelSymbolSecondKind(int i, int j, int k, const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N>& g = *this;

			// TODO - depends on covar-contravar of metric tensor!
			Real gamma_ijk = 0.0;
			for (int l = 0; l < N; l++)
			{
				Real coef1 = Derivation::DerivePartial<N>(g, i, l, j, pos, nullptr);
				Real coef2 = Derivation::DerivePartial<N>(g, j, l, i, pos, nullptr);
				Real coef3 = Derivation::DerivePartial<N>(g, i, j, l, pos, nullptr);

				gamma_ijk += 0.5 * g.Component(i, l, pos) * (coef1 + coef2 - coef3);
			}
			return gamma_ijk;
		}

		VectorN<Real, N> CovariantDerivativeContravar(const IVectorFunction<N>& func, int j, 
																									const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N> &g = *this;

			VectorN<Real, N> ret;
			VectorN<Real, N> vec_val = func(pos);

			for (int i = 0; i < N; i++) {
				Real comp_val = Derivation::DeriveVecPartial<N>(func, i, j, pos, nullptr);

				for (int k = 0; k < N; k++)
					comp_val += GetChristoffelSymbolSecondKind(i, k, j, pos) * vec_val[k];

				ret[i] = comp_val;
			}
			return ret;
		}
		Real CovariantDerivativeContravarComp(const IVectorFunction<N>& func, int i, int j, 
																					const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N>& g = *this;

			Real ret = Derivation::DeriveVecPartial<N>(func, i, j, pos, nullptr);

			for (int k = 0; k < N; k++)
				ret += GetChristoffelSymbolSecondKind(i, k, j, pos) * func(pos)[k];

			return ret;
		}

		VectorN<Real, N> CovariantDerivativeCovar(const IVectorFunction<N>& func, int j, 
																							const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N>& g = *this;

			VectorN<Real, N> ret;
			VectorN<Real, N> vec_val = func(pos);

			for (int i = 0; i < N; i++) {
				Real comp_val = Derivation::DeriveVecPartial<N>(func, i, j, pos, nullptr);

				for (int k = 0; k < N; k++)
					comp_val -= GetChristoffelSymbolSecondKind(k, i, j, pos) * vec_val[k];

				ret[i] = comp_val;
			}
			return ret;
		}
		Real CovariantDerivativeCovarComp(const IVectorFunction<N>& func, int i, int j, 
																			const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N>& g = *this;

			Real comp_val = Derivation::DeriveVecPartial<N>(func, i, j, pos, nullptr);

			for (int k = 0; k < N; k++)
				comp_val -= GetChristoffelSymbolSecondKind(k, i, j, pos) * func(pos)[k];

			return comp_val;
		}
	};

	class MetricTensorCartesian3D : public MetricTensorField<3>
	{
	public:
		MetricTensorCartesian3D() : MetricTensorField<3>(2, 0) { }

		Real Component(int i, int j, const VectorN<Real, 3>& pos) const
		{
			if (i == j)
				return 1.0;
			else
				return 0.0;
		}
	};

	class MetricTensorSpherical : public MetricTensorField<3>
	{
	public:
		MetricTensorSpherical() : MetricTensorField<3>(0, 2) { }

		virtual  Real Component(int i, int j, const VectorN<Real, 3>& pos) const override
		{
			if (i == 0 && j == 0)
				return 1.0;
			else if (i == 1 && j == 1)
				return POW2(pos[0]);
			else if (i == 2 && j == 2)
				return pos[0] * pos[0] * sin(pos[1]) * sin(pos[1]);
			else
				return 0.0;
		}
	};
	class MetricTensorSphericalContravar : public MetricTensorField<3>
	{
	public:
		MetricTensorSphericalContravar() : MetricTensorField<3>(2, 0) { }

		virtual Real Component(int i, int j, const VectorN<Real, 3>& pos) const
		{
			if (i == 0 && j == 0)
				return 1.0;
			else if (i == 1 && j == 1)
				return 1 / (pos[0] * pos[0]);
			else if (i == 2 && j == 2)
				return 1 / (pos[0] * pos[0] * sin(pos[1]) * sin(pos[1]));
			else
				return 0.0;
		}
	};
	class MetricTensorCylindrical : public MetricTensorField<3>
	{
	public:
		MetricTensorCylindrical() : MetricTensorField<3>(2, 0) { }

		virtual Real Component(int i, int j, const VectorN<Real, 3>& pos) const
		{
			if (i == 0 && j == 0)
				return 1.0;
			else if (i == 1 && j == 1)
				return pos[0] * pos[0];
			else if (i == 2 && j == 2)
				return 1.0;
			else
				return 0.0;
		}
	};

	template<typename VectorFrom, typename VectorTo, int N>
	class MetricTensorFromCoordTransf : public MetricTensorField<N>
	{
		const ICoordTransfWithInverse<VectorFrom, VectorTo, N>& _coordTransf;

	public:
		MetricTensorFromCoordTransf(ICoordTransfWithInverse<VectorFrom, VectorTo, N>& inTransf) : _coordTransf(inTransf)
		{ }

		virtual Real Component(int i, int j, const VectorN<Real, N>& pos) const
		{
			Real g_ij = 0.0;
			for (int k = 0; k < N; k++)
			{
				auto der_k_by_i = Derivation::DerivePartial<N>(_coordTransf.coordTransfFunc(k), i, pos, nullptr);
				auto der_k_by_j = Derivation::DerivePartial<N>(_coordTransf.coordTransfFunc(k), j, pos, nullptr);

				g_ij += der_k_by_i * der_k_by_j;
			}
			return g_ij;
		}
	};

	class MetricTensorMinkowski : public MetricTensorField<4>
	{
	public:
		MetricTensorMinkowski() : MetricTensorField<4>(2, 0) {}

		virtual Real Component(int i, int j, const VectorN<Real, 4>& pos) const override
		{
			if (i == 0 && j == 0)
				return -1.0;
			else if (i == 1 && j == 1)
				return 1.0;
			else if (i == 2 && j == 2)
				return 1.0;
			else if (i == 3 && j == 3)
				return 1.0;
			else
				return 0.0;
		}
	};
}
#endif