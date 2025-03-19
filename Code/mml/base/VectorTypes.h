#if !defined MML_VECTOR_TYPES_H__
#define MML_VECTOR_TYPES_H__

#include "MMLBase.h"

#include "base/VectorN.h"
#include "base/Geometry.h"

namespace MML
{
	class Vector2Cartesian : public VectorN<Real, 2>
	{
	public:
		Vector2Cartesian() {}
		Vector2Cartesian(Real x, Real y)
		{
			_val[0] = x;
			_val[1] = y;
		}
		Vector2Cartesian(const VectorN<Real, 2>& b) : VectorN<Real, 2>{ b[0], b[1] } {}
		Vector2Cartesian(const Point2Cartesian& a, const Point2Cartesian& b)
		{
			_val[0] = b.X() - a.X();
			_val[1] = b.Y() - a.Y();
		}

		Real  X() const { return _val[0]; }
		Real& X() { return _val[0]; }
		Real  Y() const { return _val[1]; }
		Real& Y() { return _val[1]; }
			
		// For Cartesian vector, we will enable operator* to represent standard scalar product
		Real operator*(const Vector2Cartesian& b) const
		{
			return X()*b.X() + Y()*b.Y();
		}

		Vector2Cartesian operator*(Real b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = _val[i] * b;
			return ret;
		}
		Vector2Cartesian operator/(Real b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = _val[i] / b;
			return ret;
		}
		friend Vector2Cartesian operator*(Real a, const Vector2Cartesian& b)
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = a * b[i];
			return ret;
		}

		// equality operators
		bool operator==(const Vector2Cartesian& b) const
		{
			return (X() == b.X()) && (Y() == b.Y());
		}
		bool operator!=(const Vector2Cartesian& b) const
		{
			return (X() != b.X()) || (Y() != b.Y());
		}
		bool IsEqual(const Vector2Cartesian& b, Real absEps = Defaults::Vec2CartIsEqualTolerance) const
		{
			return (std::abs(X() - b.X()) < absEps) && (std::abs(Y() - b.Y()) < absEps);
		}

		Vector2Cartesian GetAsUnitVector() const
		{
			VectorN<Real, 2> res = (*this) / NormL2();

			return Vector2Cartesian(res[0], res[1]);
		}
		Vector2Cartesian GetAsUnitVectorAtPos(const Vector2Cartesian& pos) const
		{
			return Vector2Cartesian{ (*this) / NormL2() };
		}

		friend Real ScalarProduct(const Vector2Cartesian& a, const Vector2Cartesian& b)
		{
			return a * b;
		}

		friend Point2Cartesian operator+(const Point2Cartesian& a, const Vector2Cartesian& b) { return Point2Cartesian(a.X() + b[0], a.Y() + b[1]); }
		friend Point2Cartesian operator-(const Point2Cartesian& a, const Vector2Cartesian& b) { return Point2Cartesian(a.X() - b[0], a.Y() - b[1]); }
	};

	class Vector2Polar : public VectorN<Real, 2>
	{
	public:
		Real  R() const { return _val[0]; }
		Real& R() { return _val[0]; }
		Real  Phi() const { return _val[1]; }
		Real& Phi() { return _val[1]; }

		Vector2Polar() {}
		Vector2Polar(Real r, Real phi)
		{
			_val[0] = r;
			_val[1] = phi;
		}
		Vector2Polar(const VectorN<Real, 2>& b) : VectorN<Real, 2>{ b[0], b[1] } {}

		Vector2Polar GetAsUnitVectorAtPos(const Vector2Polar& pos) const
		{
			// TODO 1.1 - BIG!!!
			return Vector2Polar{ (*this) / NormL2() };
		}
	};

	class Vector3Cartesian : public VectorN<Real, 3>
	{
	public:
		Real  X() const { return _val[0]; }
		Real& X()				{ return _val[0]; }
		Real  Y() const { return _val[1]; }
		Real& Y()				{ return _val[1]; }
		Real  Z() const { return _val[2]; }
		Real& Z()				{ return _val[2]; }

		Vector3Cartesian() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Cartesian(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b } {}
		Vector3Cartesian(Real x, Real y, Real z) : VectorN<Real, 3>{ x, y, z } {}
		Vector3Cartesian(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }
		Vector3Cartesian(const Point3Cartesian& a, const Point3Cartesian& b)
		{
			_val[0] = b.X() - a.X();
			_val[1] = b.Y() - a.Y();
			_val[2] = b.Z() - a.Z();
		}
		
		Point3Cartesian getAsPoint()
		{
			return Point3Cartesian(_val[0], _val[1], _val[2]);
		}
		
		Vector3Cartesian GetAsUnitVector() const
		{
			return Vector3Cartesian{ (*this) / NormL2() };
		}
		Vector3Cartesian GetAsUnitVectorAtPos(const Vector3Cartesian& pos) const
		{
			return Vector3Cartesian{ (*this) / NormL2() };
		}

		// For Cartesian vector, we will enable operator* to represent standard scalar product
		Real operator*(const Vector3Cartesian& b) const
		{
			return X()*b.X() + Y()*b.Y() + Z()*b.Z();
		}

		Vector3Cartesian operator*(Real b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = _val[i] * b;
			return ret;
		}
		friend Vector3Cartesian operator*(Real a, const Vector3Cartesian& b)
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = a * b[i];
			return ret;
		}
		Vector3Cartesian operator/(Real b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = _val[i] / b;
			return ret;
		}

		// equality operators
		bool operator==(const Vector3Cartesian& b) const
		{
			return (X() == b.X()) && (Y() == b.Y()) && (Z() == b.Z());
		}
		bool operator!=(const Vector3Cartesian& b) const
		{
			return (X() != b.X()) || (Y() != b.Y()) || (Z() != b.Z());
		}
		bool IsEqual(const Vector3Cartesian& b, Real absEps = Defaults::Vec3CartIsEqualTolerance) const
		{
			return (std::abs(X() - b.X()) < absEps) && (std::abs(Y() - b.Y()) < absEps) && (std::abs(Z() - b.Z()) < absEps);
		}

		friend Point3Cartesian operator+(const Point3Cartesian& a, const Vector3Cartesian& b) { return Point3Cartesian(a.X() + b[0], a.Y() + b[1], a.Z() + b[2]); }
		friend Point3Cartesian operator-(const Point3Cartesian& a, const Vector3Cartesian& b) { return Point3Cartesian(a.X() - b[0], a.Y() - b[1], a.Z() - b[2]); }

		bool IsParallelTo(const Vector3Cartesian& b, Real eps = Defaults::Vec3CartIsParallelTolerance) const
		{
			Real norm1 = NormL2();
			Real norm2 = b.NormL2();

			return std::abs(X() / norm1 - b.X() / norm2) < eps &&
				std::abs(Y() / norm1 - b.Y() / norm2) < eps &&
				std::abs(Z() / norm1 - b.Z() / norm2) < eps;
		}
		bool IsPerpendicularTo(const Vector3Cartesian& b, Real eps = 1e-15) const
		{
			if (std::abs(ScalarProduct(*this, b)) < eps)
				return true;
			else
				return false;
		}
		Real AngleToVector(const Vector3Cartesian& b)
		{
			Real cos_phi = ScalarProduct(*this, b) / (NormL2() * b.NormL2());

			return acos(cos_phi);
		}

		friend Real ScalarProduct(const Vector3Cartesian& a, const Vector3Cartesian& b)
		{
			return a * b;
		}
		friend Vector3Cartesian VectorProduct(const Vector3Cartesian& a, const Vector3Cartesian& b)
		{
			Vector3Cartesian ret;

			ret.X() = a.Y() * b.Z() - a.Z() * b.Y();
			ret.Y() = a.Z() * b.X() - a.X() * b.Z();
			ret.Z() = a.X() * b.Y() - a.Y() * b.X();

			return ret;
		}
	};

	class Vector3Spherical : public VectorN<Real, 3>
	{
	public:
		Real  R()     const { return _val[0]; }
		Real& R()						{ return _val[0]; }
		Real  Theta() const { return _val[1]; }
		Real& Theta()				{ return _val[1]; }
		Real  Phi()   const { return _val[2]; }
		Real& Phi()					{ return _val[2]; }

		Vector3Spherical() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Spherical(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b[0], b[1], b[2] } {}
		Vector3Spherical(Real r, Real theta, Real phi) : VectorN<Real, 3>{ r, theta, phi } {}
		Vector3Spherical(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }

		Vector3Spherical GetAsUnitVectorAtPos(const Vector3Spherical& pos) const
		{
			// TODO 1.1 - VERIFY this!!!
			return Vector3Spherical{ R(), Theta() / pos.R(), Phi() / (pos.R() * sin(pos.Theta())) };
		}

		std::ostream& PrintDeg(std::ostream& stream, int width, int precision) const
		{
			stream << "[ ";
			stream << std::fixed << std::setw(width) << std::setprecision(precision);
			stream << R();
			stream << ", " << Theta() * 180.0 / Constants::PI;
			stream << ", " << Phi() * 180.0 / Constants::PI << " ]" << std::endl;

			return stream;
		}
	};

	class Vector3Cylindrical : public VectorN<Real, 3>
	{
	public:
		Real    R()   const { return _val[0]; }
		Real& R()						{ return _val[0]; }
		Real    Phi() const { return _val[1]; }
		Real& Phi()					{ return _val[1]; }
		Real    Z()   const { return _val[2]; }
		Real& Z()						{ return _val[2]; }

		Vector3Cylindrical() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Cylindrical(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b[0], b[1], b[2] } {}
		Vector3Cylindrical(Real r, Real phi, Real z) : VectorN<Real, 3>{ r, phi, z } {}
		Vector3Cylindrical(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }

		Vector3Cylindrical GetAsUnitVectorAtPos(const Vector3Cylindrical& pos) const
		{
			return Vector3Cylindrical{ R(), Phi() / pos.R(), Z() };
		}
	};

	class Vector4Lorentz : public VectorN<Real, 4>
	{
	public:
		Real  T() const { return _val[0]; }
		Real& T()				{ return _val[0]; }
		Real  X() const { return _val[1]; }
		Real& X()				{ return _val[1]; }
		Real  Y() const { return _val[2]; }
		Real& Y()				{ return _val[2]; }
		Real  Z() const { return _val[3]; }
		Real& Z()				{ return _val[3]; }

		Vector4Lorentz() : VectorN<Real, 4>{ 0.0, 0.0, 0.0, 0.0 } {}
		Vector4Lorentz(std::initializer_list<Real> list) : VectorN<Real, 4>(list) { }
	};

	typedef Vector2Cartesian    Vec2Cart;
	typedef Vector3Cartesian    Vec3Cart;
	typedef Vector3Spherical    Vec3Sph;
	typedef Vector3Cylindrical  Vec3Cyl;
	typedef Vector4Lorentz      Vec4Lor;
}
#endif
