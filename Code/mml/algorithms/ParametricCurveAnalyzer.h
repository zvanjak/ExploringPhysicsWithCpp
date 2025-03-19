 #if !defined MML_DIFF_GEOMETRY_ALGORITHMS_H
#define MML_DIFF_GEOMETRY_ALGORITHMS_H

#include "MMLBase.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/Geometry3D.h"

#include "core/Derivation.h"

#include "core/Integration/PathIntegration.h"

namespace MML
{
	// TODO - LOW, Curve analyzer - cuspoid i slicne tocke

	class ParametricCurveAnalyzer
	{
	public:
		// helper class that returns unit tangent vector
		template<int N>
		class CurveTangentUnit : public IParametricCurve<N>
		{
			const IParametricCurve<N>& _curve;
		public:
			CurveTangentUnit(const IParametricCurve<N>& curve) : _curve(curve) {}

			VectorN<Real, N> operator()(Real t) const
			{
				auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
				return tangent_vec / tangent_vec.NormL2();
			}
		};

		// Functions for calculating curve parameters
		template<int N>
		static VectorN<Real, N> getTangent(const IParametricCurve<N>& curve, Real t)
		{
			return Derivation::DeriveCurve<N>(curve, t, nullptr);
		}
		template<int N>
		static VectorN<Real, N> getTangentUnit(const IParametricCurve<N>& curve, Real t)
		{
			auto tangent = getTangent(curve, t);
			return tangent / tangent.NormL2();
		}

		template<int N>
		static VectorN<Real, N> getNormal(const IParametricCurve<N>& curve, Real t)
		{
			return Derivation::DeriveCurveSec<N>(curve, t, nullptr);
		}

		template<int N>
		static VectorN<Real, N> getNormalUnit(const IParametricCurve<N>& curve, Real t)
		{
			auto normal = getNormal(curve, t);
			return normal / normal.NormL2();
		}
		template<int N>
		static VectorN<Real, N> getPrincipalNormal(const IParametricCurve<N>& curve, Real t)
		{
			auto y_der_1 = Vector3Cartesian(getTangent(curve, t));
			auto y_der_2 = Vector3Cartesian(Derivation::DeriveCurveSec<3>(curve, t, nullptr));

			Vector3Cartesian vec_prod1 = VectorProduct(y_der_2, y_der_1);
			Vector3Cartesian res_vec = VectorProduct(y_der_1, vec_prod1);

			return res_vec / (y_der_1.NormL2() * vec_prod1.NormL2());
		}

		template<int N>
		static VectorN<Real, N> getBinormal(const IParametricCurve<N>& curve, Real t)
		{
			auto y_der_1 = Vector3Cartesian(getTangent(curve, t));
			auto y_der_2 = Vector3Cartesian(Derivation::DeriveCurveSec<3>(curve, t, nullptr));

			Vector3Cartesian vec_prod1 = VectorProduct(y_der_2, y_der_1);
			return vec_prod1 / vec_prod1.NormL2();
		}

		template<int N>
		static VectorN<Real, N> getCurvatureVector(const IParametricCurve<N>& curve, Real t)
		{
			auto y_der_1 = getTangent(curve, t);
			auto y_der_2 = Derivation::DeriveCurveSec<3>(curve, t, nullptr);

			Real res1 = pow(y_der_1.NormL2(), -2.0);
			auto   vec2 = y_der_2 - res1 * Utils::ScalarProduct<N>(y_der_1, y_der_2) * y_der_1;

			return vec2 / res1;
		}

		template<int N>
		static Real getCurvature(const IParametricCurve<N>& curve, Real t)
		{
			auto y_der_1 = getTangent(curve, t);
			auto y_der_2 = Derivation::DeriveCurveSec<3>(curve, t, nullptr);

			Real res1 = pow(y_der_1.NormL2(), -2.0);
			auto vec2 = y_der_2 - res1 * Utils::ScalarProduct(y_der_1, y_der_2) * y_der_1;
			Real res2 = vec2.NormL2();

			return res1 * res2;
		}

		static Real getCurvature3(const IParametricCurve<3>& curve, Real t)
		{
			auto curve_first_der = Vector3Cartesian(getTangent(curve, t));
			auto curve_sec_der = Vector3Cartesian(Derivation::DeriveCurveSec<3>(curve, t, nullptr));

			auto prod = VectorProduct(curve_first_der, curve_sec_der);

			return prod.NormL2() / pow(curve_first_der.NormL2(), 3);
		}

		static Real getTorsion3(const IParametricCurve<3>& curve, Real t)
		{
			auto curve_first_der = Vector3Cartesian(getTangent(curve, t));
			auto curve_sec_der = Vector3Cartesian(Derivation::DeriveCurveSec<3>(curve, t, nullptr));
			auto curve_third_der = Vector3Cartesian(Derivation::DeriveCurveThird<3>(curve, t, nullptr));

			auto prod = VectorProduct(curve_first_der, curve_sec_der);

			Real temp = ScalarProduct(prod, curve_third_der);

			return -temp / pow(prod.NormL2(), 2);
		}

		static Plane3D getOsculationPlane(const IParametricCurve<3>& curve, Real t)
		{
			return Plane3D(Vector3Cartesian(curve(t)).getAsPoint(), Vector3Cartesian(getNormal(curve, t)));
		}
		static Plane3D getNormalPlane(const IParametricCurve<3>& curve, Real t)
		{
			return Plane3D(Vector3Cartesian(curve(t)).getAsPoint(), Vector3Cartesian(getTangentUnit(curve, t)));
		}
		static Plane3D getRectifyingPlane(const IParametricCurve<3>& curve, Real t)
		{
			return Plane3D(Vector3Cartesian(curve(t)).getAsPoint(), Vector3Cartesian(getBinormal(curve, t)));
		}

		static void getMovingTrihedron(const IParametricCurve<3>& curve, Real t, Vector3Cartesian& tangent, Vector3Cartesian& normal, Vector3Cartesian& binormal)
		{
			tangent = Vector3Cartesian(getTangentUnit(curve, t));
			normal = Vector3Cartesian(getPrincipalNormal(curve, t));
			binormal = Vector3Cartesian(getBinormal(curve, t));
		}

		static bool isArcLengthParametrized(const IParametricCurve<3>& curve, Real t1, Real t2)
		{
			int numPnt = 100;
			Real delta = (t2 - t1) / numPnt;
			for (Real t = t1 + delta; t < t2; t += delta)
			{
				Real len = PathIntegration::ParametricCurveLength(curve, t1, t);
				if (fabs(len - (t - t1)) > 1e-03)
					return false;
			}

			return true;
		}
	};
}

#endif
