#if !defined MML_COORD_TRANSF_SPHERICAL_H
#define MML_COORD_TRANSF_SPHERICAL_H

#include "MMLBase.h"

#include "core/CoordTransf.h"


namespace MML
{
	// TODO 0.9 - VIDJETI STO S MATH vs PHY konvencijama o redoslijedu koordinata
	class CoordTransfSphericalToCartesian : public CoordTransfWithInverse<Vector3Spherical, Vector3Cartesian, 3>
	{
	private:
		// q[0] = r     - radial distance
		// q[1] = theta - inclination
		// q[2] = phi   - azimuthal angle
		static Real x(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * cos(q[2]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * sin(q[2]); }
		static Real z(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }

		// q[0] = x
		// q[1] = y
		// q[2] = z
		static Real r(const VectorN<Real, 3>& q)		 { return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]); }
		static Real theta(const VectorN<Real, 3>& q) { return acos(q[2] / sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2])); }
		static Real phi(const VectorN<Real, 3>& q)	 { return atan2(q[1], q[0]); }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{x},
																								 ScalarFunction<3>{y},
																								 ScalarFunction<3>{z}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{r},
																												ScalarFunction<3>{theta},
																												ScalarFunction<3>{phi}
		};
	public:
		Vector3Cartesian     transf(const Vector3Spherical& q)				const { return Vector3Cartesian{ x(q), y(q), z(q) }; }
		Vector3Spherical     transfInverse(const Vector3Cartesian& q) const { return Vector3Spherical{ r(q), theta(q), phi(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i)				const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }

		Vector3Cartesian getBasisVectorExplicit(int ind, const Vector3Spherical& pos)
		{
			const Real r = pos[0];
			const Real theta = pos[1];
			const Real phi = pos[2];
			switch (ind)
			{
			case 0: return Vector3Cartesian{      sin(theta) * cos(phi),     sin(theta) * sin(phi),      cos(theta) };
			case 1: return Vector3Cartesian{  r * cos(theta) * cos(phi), r * cos(theta) * sin(phi), -r * sin(theta) };
			case 2: return Vector3Cartesian{ -r * sin(theta) * sin(phi), r * sin(theta) * cos(phi),						  0.0 };
			default:
				return Vector3Cartesian{ 0.0, 0.0, 0.0 };
			}
		}
		Vector3Cartesian getUnitBasisVectorExplicit(int ind, const Vector3Spherical& pos)
		{
			const Real r = pos[0];
			const Real theta = pos[1];
			const Real phi = pos[2];
			switch (ind)
			{
			case 0: return Vector3Cartesian{ sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) };
			case 1: return Vector3Cartesian{ cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta) };
			case 2: return Vector3Cartesian{ -sin(phi), cos(phi), 0.0 };
			default:
				return Vector3Cartesian{ 0.0, 0.0, 0.0 };
			}
		}

		Vector3Spherical getInverseBasisVectorExplicit(int ind, const Vector3Spherical& pos)
		{
			const Real r = pos[0];
			const Real theta = pos[1];
			const Real phi = pos[2];
			switch(ind)
			{
			case 0: return Vector3Spherical{ sin(theta) * cos(phi), r * cos(theta) * cos(phi), -r * sin(theta) * sin(phi) };
			case 1: return Vector3Spherical{ sin(theta) * sin(phi), r * cos(theta) * sin(phi),  r * sin(theta) * cos(phi) };
			case 2: return Vector3Spherical{ cos(theta)           ,-r * sin(theta)           ,                        0.0 };
			default: 
				return Vector3Spherical{ 0.0, 0.0, 0.0 };
			}
		}
		Vector3Spherical getInverseUnitBasisVectorExplicit(int ind, const Vector3Spherical& pos)
		{
			const Real r = pos[0];
			const Real theta = pos[1];
			const Real phi = pos[2];
			switch(ind)
			{
			case 0: return Vector3Spherical{ sin(theta) * cos(phi),  cos(theta) * cos(phi), -sin(phi) };
			case 1: return Vector3Spherical{ sin(theta) * sin(phi),  cos(theta) * sin(phi),  cos(phi) };
			case 2: return Vector3Spherical{ cos(theta)           , -sin(theta)           ,  0.0 };
			default: 
				return Vector3Spherical{ 0.0, 0.0, 0.0 };
			}
		}
	};

	class CoordTransfCartesianToSpherical : public CoordTransfWithInverse<Vector3Cartesian, Vector3Spherical, 3>
	{
	private:
		// q[0] = x
		// q[1] = y
		// q[2] = z
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]); }
		static Real theta(const VectorN<Real, 3>& q) { return acos(q[2] / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2])); }
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }

		// q[0] = r     - radial distance
		// q[1] = theta - inclination
		// q[2] = phi   - azimuthal angle
		static Real x(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * cos(q[2]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * sin(q[2]); }
		static Real z(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{r},
																								 ScalarFunction<3>{theta},
																								 ScalarFunction<3>{phi}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{x},
																												ScalarFunction<3>{y},
																												ScalarFunction<3>{z}
		};
	public:
		Vector3Spherical     transf(const Vector3Cartesian& q) const { return Vector3Spherical{ r(q), theta(q), phi(q) }; }
		Vector3Cartesian     transfInverse(const Vector3Spherical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};


	static CoordTransfSphericalToCartesian      CoordTransfSpherToCart;
	static CoordTransfCartesianToSpherical      CoordTransfCartToSpher;
}

#endif
