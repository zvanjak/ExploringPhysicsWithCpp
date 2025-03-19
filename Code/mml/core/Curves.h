#if !defined  MML_CURVES_H
#define MML_CURVES_H

#include "MMLBase.h"

#include "base/BaseUtils.h"
#include "base/VectorN.h"
#include "base/Vector.h"

#include "base/Function.h"
#include "base/Geometry3D.h"

#include "core/Derivation.h"

namespace MML
{
	namespace Curves
	{
		/////////////////////////////              CARTESIAN PLANAR CURVES                  ///////////////////////////////
		// abstract class, providing basic curve formulas in 2D
		class ICurveCartesian2D : public IParametricCurve<2>
		{
		public:
			Vec2Cart getTangent(Real t)
			{
				return Derivation::DeriveCurve<2>(*this, t, nullptr);
			}
			Vec2Cart getTangentUnit(Real t)
			{
				return getTangent(t).GetAsUnitVector();
			}
			Vec2Cart getNormal(Real t)
			{
				return Derivation::DeriveCurveSec<2>(*this, t, nullptr);
			}
			Vec2Cart getNormalUnit(Real t)
			{
				return getNormal(t).GetAsUnitVector();
			}
		};

		class Circle2DCurve : public ICurveCartesian2D
		{
			Real _radius;
			Pnt2Cart _center;
		public:
			Circle2DCurve() : _radius(1), _center(0, 0) {}
			Circle2DCurve(Real radius) : _radius(radius), _center(0,0) {}
			Circle2DCurve(Real radius, const Pnt2Cart& center) : _radius(radius), _center(center) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 2> operator()(Real t) const { 
				return MML::VectorN<Real, 2>{_center.X()+_radius* cos(t), _center.Y()+_radius* sin(t)}; 
			}
		};

		class LogSpiralCurve : public ICurveCartesian2D
		{
			Real _lambda, _c;
		public:
			LogSpiralCurve() : _lambda(-1), _c(1) {}
			LogSpiralCurve(Real lambda, Real c) : _lambda(lambda), _c(c) {
				if (lambda >= 0) throw std::invalid_argument("LogSpiralCurve: lambda must be negative.");
				if (c == 0) throw std::invalid_argument("LogSpiralCurve: c must not be zero.");
			}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{exp(_lambda* t)* cos(t), exp(_lambda* t)* sin(t)}; }
		};

		class LemniscateCurve : public ICurveCartesian2D
		{
		public:
			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{cos(t) / (1 + sin(t) * sin(t)), sin(t)* cos(t) / (1 + sin(t) * sin(t))}; }
		};

		class DeltoidCurve : public ICurveCartesian2D
		{
			int _n;
		public:
			DeltoidCurve() : _n(1) {}
			DeltoidCurve(int n) : _n(n) {}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{2 * _n * cos(t) * (1 + cos(t)), 2 * _n * sin(t) * (1 - cos(t))}; }
		};

		class AstroidCurve : public ICurveCartesian2D
		{
			Real _c;
		public:
			AstroidCurve() : _c(1) {}
			AstroidCurve(Real c) : _c(c) {
				if (c <= 0) throw std::invalid_argument("AstroidCurve: c must be positive.");
			}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{_c* cos(t)* cos(t)* cos(t), _c* sin(t)* sin(t)* sin(t)}; }
		};

		class EpitrochoidCurve : public ICurveCartesian2D
		{
			Real _radius, _c;
			int _n;
		public:
			EpitrochoidCurve() : _radius(1), _c(1), _n(1) {}
			EpitrochoidCurve(Real radius, Real c, int n) : _radius(radius), _c(c), _n(n) {}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{cos(t) - _c * cos(_n * t), sin(t) - _c * sin(_n * t) }; }
		};

		class ArchimedeanSpiralCurve : public ICurveCartesian2D
		{
			Real _a;
		public:
			ArchimedeanSpiralCurve() : _a(1) {}
			ArchimedeanSpiralCurve(Real a) : _a(a) {}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{_a* t* cos(t), _a* t* sin(t)}; }
		};

		///////////////////////////////             CARTESIAN SPACE CURVES                  ///////////////////////////////
		// abstract class, providing basic curve formulas in 3D
		class ICurveCartesian3D : public IParametricCurve<3>
		{
		public:
			Vec3Cart getTangent(Real t) const
			{
				return Derivation::DeriveCurve<3>(*this, t, nullptr);
			}
			Vec3Cart getTangentUnit(Real t) const
			{
				return getTangent(t).GetAsUnitVector();
			}
			Vec3Cart getNormal(Real t) const
			{
				return Derivation::DeriveCurveSec<3>(*this, t, nullptr);
			}
			Vec3Cart getNormalUnit(Real t) const
			{
				return getNormal(t).GetAsUnitVector();
			}
			Vec3Cart getBinormal(Real t) const
			{
				Vec3Cart tangent = getTangentUnit(t);
				Vec3Cart normal  = getNormalUnit(t);

				return VectorProduct(tangent, normal);
			}

			Vec3Cart getCurvatureVector(Real t) const
			{
				Vec3Cart y_der_1 = getTangent(t);
				Vec3Cart y_der_2 = getNormal(t);

				Real		 res1 = pow(y_der_1.NormL2(), -2.0);
				Vec3Cart vec2 = y_der_2 - res1 * Utils::ScalarProduct(y_der_1, y_der_2) * y_der_1;

				return vec2 / res1;
			}
			virtual Real getCurvature(Real t) const
			{
				Vec3Cart y_der_1 = getTangent(t);
				Real				res1 = pow(y_der_1.NormL2(), -2.0);
				
				Vec3Cart y_der_2 = getNormal(t);
				Vec3Cart    vec2 = y_der_2 - res1 * Utils::ScalarProduct(y_der_1, y_der_2) * y_der_1;

				Real res2 = vec2.NormL2();

				return res1 * res2;
			}
			virtual Real getTorsion( Real t) const
			{
				Vec3Cart r_prime = getTangent(t); // First derivative
				Vec3Cart r_double_prime = getNormal(t); // Second derivative
				Vec3Cart r_triple_prime = Derivation::DeriveCurveThird<3>(*this, t, nullptr); // Third derivative

				Vec3Cart cross_product = VectorProduct(r_prime, r_double_prime);
				Real numerator = Utils::ScalarProduct(cross_product, r_triple_prime);
				Real denominator = pow(cross_product.NormL2(), 2);

				return numerator / denominator;
			}
		};

			//Real getCurvature3(Real t) const
			//{
			//	Vec3Cart curve_first_der = getTangent(t);
			//	Vec3Cart curve_sec_der	 = getNormal(t);

			//	auto prod = VectorProduct(curve_first_der, curve_sec_der);

			//	return prod.NormL2() / pow(curve_first_der.NormL2(), 3);
			//}


		// concrete class, that can be initialized with function pointers or lambdas
		class CurveCartesian3D : public ICurveCartesian3D
		{
			Real _minT;
			Real _maxT;
			VectorN<Real, 3>(*_func)(Real);
		public:
			CurveCartesian3D(VectorN<Real, 3>(*inFunc)(Real)) 
				: _func(inFunc), _minT(Constants::NegInf), _maxT(Constants::PosInf) {}
			CurveCartesian3D(Real minT, Real maxT, VectorN<Real, 3>(*inFunc)(Real)) 
				: _func(inFunc), _minT(minT), _maxT(maxT) {}

			Real getMinT() const { return _minT; }
			Real getMaxT() const { return _maxT; }

			virtual VectorN<Real, 3> operator()(Real x) const { return _func(x); }
		};

		class LineCurve : public ICurveCartesian3D
		{
			Line3D  _line;
			Real _minT;
			Real _maxT;
		public:
			LineCurve(Real minT, Real maxT, const Point3Cartesian& pnt, const Vector3Cartesian& dir) : _line(pnt, dir), _minT(minT), _maxT(maxT) {}
			LineCurve(Real t1, const Point3Cartesian& pnt1, Real t2, const Point3Cartesian& pnt2)
			{
				// tocno samo ako je t1 = 0.0!!!
				_line.StartPoint() = pnt1;
				Vec3Cart dir = Vec3Cart(pnt1, pnt2);
				_line.Direction() = dir.NormL2() / (t2 - t1) * dir.GetAsUnitVector();
				_minT = t1;
				_maxT = t2;
			}

			Real getMinT() const { return _minT; }
			Real getMaxT() const { return _maxT; }

			VectorN<Real, 3> operator()(Real t) const
			{
				//if (t < _minT || t > _maxT)
				//	throw std::invalid_argument("LineCurve: t is out of range.");

				auto pnt = _line(t);
				return VectorN<Real, 3>{pnt.X(), pnt.Y(), pnt.Z()};
			}
		};

		class Circle3DXY : public ICurveCartesian3D {
			Real _radius;
		public:
			Circle3DXY() : _radius(1) {}
			Circle3DXY(Real radius) : _radius(radius) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{_radius* cos(t), _radius* sin(t), 0}; }
		};

		class Circle3DXZ : public ICurveCartesian3D {
			Real _radius;
		public:
			Circle3DXZ() : _radius(1) {}
			Circle3DXZ(Real radius) : _radius(radius) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{_radius* cos(t), 0, _radius* sin(t)}; }
		};

		class Circle3DYZ : public ICurveCartesian3D {
			Real _radius;
		public:
			Circle3DYZ() : _radius(1) {}
			Circle3DYZ(Real radius) : _radius(radius) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{0, _radius* cos(t), _radius* sin(t)}; }
		};

		class HelixCurve : public ICurveCartesian3D
		{
			Real _radius, _b;
		public:
			HelixCurve() : _radius(1.0), _b(1.0) {}
			HelixCurve(Real radius, Real b) : _radius(radius), _b(b) {}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 3> operator()(Real t) const { 
				return MML::VectorN<Real, 3>{_radius * cos(t), _radius * sin(t), _b * t}; 
			}

			Real getCurvature(Real t) const { return _radius / (POW2(_radius) + POW2(_b)); }
			Real getTorsion(Real t) const { return _b / (POW2(_radius) + POW2(_b)); }
		};

		class TwistedCubicCurve : public ICurveCartesian3D
		{
		public:
			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{t, t* t, t* t* t}; }
		};

		class ToroidalSpiralCurve : public ICurveCartesian3D
		{
			int _n;
			Real _scale = 1.0;

		public:
			ToroidalSpiralCurve() : _n(1) {}
			ToroidalSpiralCurve(int n) : _n(n) {}
			ToroidalSpiralCurve(Real scale) : _n(1), _scale(scale) {}
			ToroidalSpiralCurve(int n, Real scale) : _n(n), _scale(scale) {}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{(_scale* (4 + sin(_n * t))* cos(t)), _scale* (4 + sin(_n * t))* sin(t), _scale* cos(_n* t)}; }
		};

		// TODO - add general circle, lying in plane with given normal, with given center position


		//////////////////////////////////////////////////////////////////////////////
		// abstract class, providing basic curve formulas in 4D
		class CurveLorentz4D : public IParametricCurve<4>
		{
		};
	}
}

#endif // MML_CURVES_H