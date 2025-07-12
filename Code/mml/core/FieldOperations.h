// grad
// - cart
// - spher
// - cyl

#if !defined MML_FIELD_OPERATIONS_H
#define MML_FIELD_OPERATIONS_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Geometry.h"
#include "base/VectorN.h"
#include "base/VectorTypes.h"
#include "base/MatrixNM.h"
#include "base/Tensor.h"

#include "core/Derivation.h"
#include "core/MetricTensor.h"

namespace MML
{
	namespace ScalarFieldOperations
	{
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////                 GENERAL COORD. FIELD OPERATIONS                    ///////////////////
		template<int N>
		static VectorN<Real, N> Gradient(IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos, 
																		 const MetricTensorField<N>& metricTensorField)
		{
			VectorN<Real, N> derivsAtPoint = Derivation::DerivePartialAll<N>(scalarField, pos, nullptr);

			// TODO - depends on covar-contravar of metric tensor!!!
			Tensor2<N> metricAtPoint(2, 0);
			metricTensorField.ValueAtPoint(pos, metricAtPoint);

			VectorN<Real, N> ret = metricAtPoint * derivsAtPoint;

			return ret;
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////                   GRADIENT                     /////////////////////////////
		template<int N>
		static VectorN<Real, N> GradientCart(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos)
		{
			return Derivation::DerivePartialAll<N>(scalarField, pos, nullptr);
		}
		template<int N>
		static VectorN<Real, N> GradientCart(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos, 
																				 int der_order)
		{
			switch (der_order)
			{
			case 1: return Derivation::NDer1PartialByAll<N>(scalarField, pos, nullptr);
			case 2: return Derivation::NDer2PartialByAll<N>(scalarField, pos, nullptr);
			case 4: return Derivation::NDer4PartialByAll<N>(scalarField, pos, nullptr);
			case 6: return Derivation::NDer6PartialByAll<N>(scalarField, pos, nullptr);
			case 8: return Derivation::NDer8PartialByAll<N>(scalarField, pos, nullptr);
			default:
				throw std::invalid_argument("GradientCart: der_order must be in 1, 2, 4, 6 or 8");
			}
		}

		static Vec3Sph GradientSpher(const IScalarFunction<3>& scalarField, const Vec3Sph& pos)
		{
			Vector3Spherical ret = Derivation::DerivePartialAll<3>(scalarField, pos, nullptr);

			ret[1] = ret[1] / pos[0];
			ret[2] = ret[2] / (pos[0] * sin(pos[1]));

			return ret;
		}
		static Vec3Sph GradientSpher(const IScalarFunction<3>& scalarField, const Vec3Sph& pos, 
																 int der_order)
		{
			Vector3Spherical ret;

			switch (der_order)
			{
			case 1: ret = Derivation::NDer1PartialByAll<3>(scalarField, pos, nullptr); break;
			case 2: ret = Derivation::NDer2PartialByAll<3>(scalarField, pos, nullptr); break;
			case 4: ret = Derivation::NDer4PartialByAll<3>(scalarField, pos, nullptr); break;
			case 6: ret = Derivation::NDer6PartialByAll<3>(scalarField, pos, nullptr); break;
			case 8: ret = Derivation::NDer8PartialByAll<3>(scalarField, pos, nullptr); break;
			default:
				throw std::invalid_argument("GradientSpher: der_order must be in 1, 2, 4, 6 or 8");
			}

			ret[1] = ret[1] / pos[0];
			ret[2] = ret[2] / (pos[0] * sin(pos[1]));

			return ret;
		}

		static Vec3Cyl GradientCyl(const IScalarFunction<3>& scalarField, const Vec3Cyl& pos)
		{
			Vector3Cylindrical ret = Derivation::DerivePartialAll<3>(scalarField, pos, nullptr);

			ret[1] = ret[1] / pos[0];

			return ret;
		}
		static Vec3Cyl GradientCyl(const IScalarFunction<3>& scalarField, const Vec3Cyl& pos, 
															 int der_order)
		{
			Vector3Cylindrical ret;

			switch (der_order)
			{
			case 1: ret = Derivation::NDer1PartialByAll<3>(scalarField, pos, nullptr); break;
			case 2: ret = Derivation::NDer2PartialByAll<3>(scalarField, pos, nullptr); break;
			case 4: ret = Derivation::NDer4PartialByAll<3>(scalarField, pos, nullptr); break;
			case 6: ret = Derivation::NDer6PartialByAll<3>(scalarField, pos, nullptr); break;
			case 8: ret = Derivation::NDer8PartialByAll<3>(scalarField, pos, nullptr); break;
			default:
				throw std::invalid_argument("GradientCyl: der_order must be in 1, 2, 4, 6 or 8");
			}
			ret[1] = ret[1] / pos[0];

			return ret;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////                  LAPLACIAN                     /////////////////////////////
		template<int N>
		static Real LaplacianCart(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos)
		{
			Real lapl = 0.0;
			for (int i = 0; i < N; i++)
				lapl += Derivation::DeriveSecPartial<N>(scalarField, i, i, pos, nullptr);

			return lapl;
		}
		static Real LaplacianSpher(const IScalarFunction<3>& scalarField, const Vec3Sph& pos)
		{
			const Real r = pos.R();
			const Real phi = pos.Phi();
			const Real theta = pos.Theta();

			Real first = Derivation::DeriveSecPartial<3>(scalarField, 0, 0, pos, nullptr);
			Real second = 2 / pos.R() * Derivation::DerivePartial<3>(scalarField, 0, pos, nullptr);
			Real third = 1 / (r * r * sin(theta)) * (cos(theta) * Derivation::DerivePartial<3>(scalarField, 1, pos, nullptr) + sin(theta) * Derivation::DeriveSecPartial<3>(scalarField, 1, 1, pos, nullptr));
			Real fourth = 1 / (r * r * sin(theta) * sin(theta));

			return first + second + third;
		}
		static Real LaplacianCyl(const IScalarFunction<3>& scalarField, const Vec3Cyl& pos)
		{
			const Real r = pos[0];

			Real first = 1 / r * (Derivation::DerivePartial<3>(scalarField, 0, pos, nullptr) + r * Derivation::DeriveSecPartial<3>(scalarField, 0, 0, pos, nullptr));
			Real second = 1 / (r * r) * Derivation::DeriveSecPartial<3>(scalarField, 1, 1, pos, nullptr);
			Real third = Derivation::DeriveSecPartial<3>(scalarField, 2, 2, pos, nullptr);

			return first + second + third;
		}
	};

	namespace VectorFieldOperations
	{
		template<int N>
		static Real Divergence(const IVectorFunction<N>& vectorField, const VectorN<Real, N>& pos,
													const MetricTensorField<N>& metricTensorField)
		{
			Real div = 0.0;
			VectorN<Real, N> vec_val = vectorField(pos);

			for (int i = 0; i < N; i++)
			{
				div += Derivation::DeriveVecPartial<N>(vectorField, i, i, pos, nullptr);

				// correction for general coordinates
				for (int k = 0; k < N; k++)
				{
					div += vec_val[k] * metricTensorField.GetChristoffelSymbolSecondKind(i, i, k, pos);
				}
			}
			return div;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////                  DIVERGENCE                    /////////////////////////////
		template<int N>
		static Real DivCart(const IVectorFunction<N>& vectorField, const VectorN<Real, N>& pos)
		{
			Real div = 0.0;
			for (int i = 0; i < N; i++)
				div += Derivation::DeriveVecPartial<N>(vectorField, i, i, pos, nullptr);

			return div;
		}
		static Real DivSpher(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& x)
		{
			VectorN<Real, 3> vals = vectorField(x);

			VectorN<Real, 3> derivs;
			for (int i = 0; i < 3; i++)
				derivs[i] = Derivation::DeriveVecPartial<3>(vectorField, i, i, x, nullptr);

			Real div = 0.0;
			div += 1 / (x[0] * x[0]) * (2 * x[0] * vals[0] + x[0] * x[0] * derivs[0]);
			div += 1 / (x[0] * sin(x[1])) * (cos(x[1]) * vals[1] + sin(x[1]) * derivs[1]);
			div += 1 / (x[0] * sin(x[1])) * derivs[2];

			return div;
		}
		static Real DivCyl(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& x)
		{
			VectorN<Real, 3> vals = vectorField(x);

			VectorN<Real, 3> derivs;
			for (int i = 0; i < 3; i++)
				derivs[i] = Derivation::DeriveVecPartial<3>(vectorField, i, i, x, nullptr);

			Real div = 0.0;
			div += 1 / x[0] * (vals[0] + x[0] * derivs[0]);
			div += 1 / x[0] * derivs[1];
			div += derivs[2];

			return div;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////                     CURL                       /////////////////////////////
		static Vec3Cart CurlCart(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos)
		{
			Real dzdy = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);
			Real dydz = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);

			Real dxdz = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);
			Real dzdx = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);

			Real dydx = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);
			Real dxdy = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);

			Vector3Cartesian curl{ dzdy - dydz, dxdz - dzdx, dydx - dxdy };

			return curl;
		}
		static Vec3Sph  CurlSpher(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos)
		{
			VectorN<Real, 3> vals = vectorField(pos);

			Real dphidtheta = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);
			Real dthetadphi = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);

			Real drdphi = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);
			Real dphidr = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);

			Real dthetadr = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);
			Real drdtheta = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);

			Vector3Spherical ret;
			const Real& r = pos[0];
			const Real& theta = pos[1];
			const Real& phi = pos[2];

			ret[0] = 1 / (r * sin(theta)) * (cos(theta) * vals[2] + sin(theta) * dphidtheta - dthetadphi);
			ret[1] = 1 / r * (1 / sin(theta) * drdphi - vals[2] - r * dphidr);
			ret[2] = 1 / r * (vals[1] + r * dthetadr - drdtheta);

			return ret;
		}
		static Vec3Cyl  CurlCyl(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos)
		{
			VectorN<Real, 3> vals = vectorField(pos);

			Real dzdphi = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);
			Real dphidz = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);

			Real drdz = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);
			Real dzdr = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);

			Real dphidr = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);
			Real drdphi = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);

			Vector3Cylindrical ret{(1 / pos[0] * dzdphi - dphidz), drdz - dzdr, 1 / pos[0] * (vals[1] + pos[0] * dphidr - drdphi)};

			return ret;
		}
	};
}
#endif
