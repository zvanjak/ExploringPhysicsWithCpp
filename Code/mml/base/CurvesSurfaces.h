#if !defined  MML_CURVES_SURFACES_H
#define MML_CURVES_SURFACES_H

#include "MMLBase.h"

#include "base/VectorN.h"
#include "base/Vector.h"

#include "base/Function.h"
#include "base/Geometry3D.h"

namespace MML
{
	namespace Curves3D
	{
		///////////////////////////////             SPHERICAL SPACE CURVES                  ///////////////////////////////
		class Circle3DXYSpherical : public IParametricCurve<3> {
			Real _radius;
		public:
			Circle3DXYSpherical() : _radius(1) {}
			Circle3DXYSpherical(Real radius) : _radius(radius) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{_radius, Constants::PI/2, t}; }
		};

	} // end namespace Curves3D

	namespace Surfaces
	{
		// TODO - PlaneSurface3D, given by point and normal

		// MonkeySaddle
		class MonkeySaddle : public IParametricSurfaceRect<3>
		{
		public:
			Real getMinU() const { return -10; }
			Real getMaxU() const { return 10; }
			Real getMinW() const { return -10; }
			Real getMaxW() const { return 10; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{u, w, u* (u* u - 3 * w * w)}; }
		};
		// MobiusStrip
		class MobiusStrip : public IParametricSurfaceRect<3>
		{
		public:
			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return -1; }
			Real getMaxW() const { return 1; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{(1 + w * cos(u / 2))* cos(u), (1 + w * cos(u / 2))* sin(u), w* sin(u / 2)}; }
		};
		// Torus
		class Torus : public IParametricSurfaceRect<3>
		{
			Real _R, _r;
		public:
			Torus() : _R(1), _r(0.5) {}
			Torus(Real R, Real r) : _R(R), _r(r) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{(_R + _r * cos(w))* cos(u), (_R + _r * cos(w))* sin(u), _r* sin(w)}; }
		};
		// Sphere
		class Sphere : public IParametricSurfaceRect<3>
		{
			Real _R;
		public:
			Sphere() : _R(1) {}
			Sphere(Real R) : _R(R) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{_R* sin(u)* cos(w), _R* sin(u)* sin(w), _R* cos(u)}; }
		};
		// Ellipsoid
		class Ellipsoid : public IParametricSurfaceRect<3>
		{
			Real _a, _b, _c;
		public:
			Ellipsoid() : _a(1), _b(1), _c(1) {}
			Ellipsoid(Real a, Real b, Real c) : _a(a), _b(b), _c(c) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{_a* sin(u)* cos(w), _b* sin(u)* sin(w), _c* cos(u)}; }
		};
		// Cylinder
		class Cylinder : public IParametricSurfaceRect<3>
		{
			Real _R, _H;
		public:
			Cylinder() : _R(1), _H(1) {}
			Cylinder(Real R, Real H) : _R(R), _H(H) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return _H; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{_R* cos(u), _R* sin(u), w}; }
		};
	} // end namespace Surfaces
}

#endif