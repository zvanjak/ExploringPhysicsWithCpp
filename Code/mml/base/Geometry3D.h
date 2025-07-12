#if !defined MML_GEOMETRY_3D_H
#define MML_GEOMETRY_3D_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/VectorN.h"
#include "base/VectorTypes.h"
#include "base/Geometry.h"

namespace MML
{
	class Line3D
	{
	private:
		Pnt3Cart  _point;
		Vec3Cart _direction;

	public:
		Line3D() {}
		// by default, direction vector is normalized to unit vector (but, it need not be such!)
		Line3D(const Pnt3Cart& pnt, const Vec3Cart dir)
		{
			// check for null vector as direction
			if (dir.X() == 0.0 && dir.Y() == 0.0 && dir.Z() == 0.0)
				throw std::runtime_error("Line3D ctor - null vector as direction");

			_point = pnt;
			_direction = dir.GetAsUnitVector();
		}
		Line3D(const Pnt3Cart& a, const Pnt3Cart& b)
		{
			// check for same points
			if (a == b)
				throw std::runtime_error("Line3D ctor - same points");

			Vec3Cart dir(a, b);
			_point = a;
			_direction = dir.GetAsUnitVector();
		}

		Pnt3Cart   StartPoint() const { return _point; }
		Pnt3Cart& StartPoint() { return _point; }

		Vec3Cart  Direction() const { return _direction; }
		Vec3Cart& Direction() { return _direction; }

		Pnt3Cart operator()(Real t) const { return _point + t * _direction; }

		bool AreEqual(const Line3D& b, Real eps = Defaults::Line3DAreEqualTolerance) const
		{
			return IsPointOnLine(b.StartPoint()) && b.IsPointOnLine(StartPoint()) && Direction().IsEqual(b.Direction(), eps);
		}
		bool operator==(const Line3D& b) const
		{
			return AreEqual(b, Defaults::Line3DAreEqualTolerance);
		}

		bool IsPerpendicular(const Line3D& b, Real eps = Defaults::Line3DIsPerpendicularTolerance) const
		{
			return ScalarProduct(Direction(), b.Direction()) < eps;
		}
		bool IsParallel(const Line3D& b, Real eps = Defaults::Line3DIsParallelTolerance) const
		{
			return Direction().IsEqual(b.Direction(), eps);
		}

		bool IsPointOnLine(const Pnt3Cart& pnt, Real eps = Defaults::Line3DIsPointOnLineTolerance) const
		{
			if (pnt == StartPoint())
				return true;	// point is start point of line

			// check if point is on line by checking if the vector from start point to point is parallel to direction vector
			Vec3Cart vecFromStartToPnt(_point, pnt);
			return vecFromStartToPnt.IsParallelTo(_direction, eps);
		}

		// distance between point and line
		Real Dist(const Pnt3Cart& pnt) const
		{
			// Bronshtein 3.394
			const Real a = pnt.X();
			const Real b = pnt.Y();
			const Real c = pnt.Z();

			const Real x1 = StartPoint().X();
			const Real y1 = StartPoint().Y();
			const Real z1 = StartPoint().Z();

			const Real l = Direction().X();
			const Real m = Direction().Y();
			const Real n = Direction().Z();

			Real numer = POW2((a - x1) * m - (b - y1) * l) + POW2((b - y1) * n - (c - z1) * m) + POW2((c - z1) * l - (a - x1) * n);
			Real denom = l * l + m * m + n * n;

			return sqrt(numer / denom);
		}
		// nearest point on line to given point
		Pnt3Cart NearestPointOnLine(const Pnt3Cart& pnt) const
		{
			// https://math.stackexchange.com/questions/1521128/given-a-line-and-a-point-in-3d-how-to-find-the-closest-point-on-the-line         
			Vec3Cart line_dir = this->Direction();
			Vec3Cart rad_vec_AP(StartPoint(), pnt);

			double t = rad_vec_AP * line_dir / POW2(line_dir.NormL2());

			return StartPoint() + t * line_dir;
		}

		// distance between two lines
		Real Dist(const Line3D& line) const
		{
			// https://math.stackexchange.com/questions/2213165/distance-between-two-lines-in-3d-space
			// https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
			Pnt3Cart  p1 = StartPoint();
			Vec3Cart d1 = Direction();
			Pnt3Cart  p2 = line.StartPoint();
			Vec3Cart d2 = line.Direction();

			Vec3Cart n = VectorProduct(d1, d2);

			if (n.IsNullVec())
			{
				// parallel lines
				return Vec3Cart(p1, p2) * d1 / d2.NormL2();
			}
			else
			{
				// skew lines
				return Abs(n * Vec3Cart(p1, p2)) / n.NormL2();
			}
		}
		// distance between two lines, while also returning nearest points on both lines
		// if lines are parallel, returns false
		bool Dist(const Line3D& line, Real& out_dist, Pnt3Cart& out_line1_pnt, Pnt3Cart& out_line2_pnt) const
		{
			// https://math.stackexchange.com/questions/2213165/distance-between-two-lines-in-3d-space
			// https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
			Pnt3Cart  p1 = StartPoint();
			Vec3Cart d1 = Direction();
			Pnt3Cart  p2 = line.StartPoint();
			Vec3Cart d2 = line.Direction();
			out_dist = 0.0;

			Vec3Cart n = VectorProduct(d1, d2);

			if (n.IsNullVec())				// check for parallel lines
			{
				out_dist = Vec3Cart(p1, p2) * d1 / d2.NormL2();
				return false;
			}
			else
			{
				// skew lines
				out_dist = Abs(n * Vec3Cart(p1, p2)) / n.NormL2();

				Real t1 = (VectorProduct(d2, n) * Vec3Cart(p1, p2)) / POW2(n.NormL2());
				Real t2 = (VectorProduct(d1, n) * Vec3Cart(p1, p2)) / POW2(n.NormL2());

				out_line1_pnt = p1 + t1 * d1;
				out_line2_pnt = p2 + t2 * d2;
			}

			return true;
		}
		// intersection of two lines
		bool Intersection(const Line3D& line, Pnt3Cart& out_inter_pnt) const
		{
			// https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
			// https://math.stackexchange.com/questions/2213165/distance-between-two-lines-in-3d-space
			Pnt3Cart  p1 = StartPoint();
			Vec3Cart d1 = Direction();
			Pnt3Cart  p2 = line.StartPoint();
			Vec3Cart d2 = line.Direction();
			Real dist = 0.0;

			Vec3Cart n = VectorProduct(d1, d2);

			if (n.IsNullVec())				// check for parallel lines
			{
				if (*this == line)
				{
					// lines are equal, return any point on the line
					out_inter_pnt = p1;
					return true;
				}
				else
					return false;
			}
			else				// skew lines
			{
				dist = Abs(n * Vec3Cart(p1, p2)) / n.NormL2();

				if (dist > Defaults::Line3DIntersectionTolerance)
					return false;	// lines are skew, no intersection

				Real t1 = (VectorProduct(d2, n) * Vec3Cart(p1, p2)) / POW2(n.NormL2());

				out_inter_pnt = p1 + t1 * d1;
			}

			return true;
			return true;
		}

		// perpendicular line that goes through givenpoint
		// direction is from given pnt to point on line that is closest to given point
		Line3D PerpendicularLineThroughPoint(const Pnt3Cart& pnt)
		{
			if (IsPointOnLine(pnt))
				throw std::runtime_error("Line3D::PerpendicularLineThroughPoint - point is on the line");

			return Line3D(pnt, Vec3Cart(pnt, NearestPointOnLine(pnt)).GetAsUnitVector());
		}
	};

	class SegmentLine3D
	{
	private:
		Pnt3Cart _point1;
		Pnt3Cart _point2;

	public:
		SegmentLine3D(Pnt3Cart pnt1, Pnt3Cart pnt2) : _point1(pnt1), _point2(pnt2)
		{
		}
		SegmentLine3D(Pnt3Cart pnt1, Vec3Cart direction, Real t)
		{
			_point1 = pnt1;
			_point2 = pnt1 + t * direction;
		}

		Pnt3Cart   StartPoint() const { return _point1; }
		Pnt3Cart& StartPoint() { return _point1; }

		Pnt3Cart   EndPoint() const { return _point2; }
		Pnt3Cart& EndPoint() { return _point2; }

		Pnt3Cart		operator()(Real t)
		{
			if (t < 0.0 || t > 1.0)
				throw std::runtime_error("SegmentLine3D::PointOnSegment - t not in [0,1]");

			return _point1 + t * Direction();
		}

		Real              Length()    const { return _point1.Dist(_point2); }
		Vec3Cart  Direction() const { return Vec3Cart(_point1, _point2); }
	};

	class Plane3D
	{
	private:
		Real _A, _B, _C, _D;

	public:
		Plane3D(const Pnt3Cart& a, const Vec3Cart& normal)
		{
			if (normal.IsNullVec())
				throw std::runtime_error("Plane3D ctor - normal is null vector");

			Vec3Cart unitNormal = normal.GetAsUnitVector();

			_A = unitNormal.X();
			_B = unitNormal.Y();
			_C = unitNormal.Z();
			_D = -(a.X() * unitNormal.X() + a.Y() * unitNormal.Y() + a.Z() * unitNormal.Z());
		}
		Plane3D(const Pnt3Cart& a, const Pnt3Cart& b, const Pnt3Cart& c)
			: Plane3D(a, VectorProduct(Vec3Cart(a, b), Vec3Cart(a, c)))
		{
		}
		// Hesse normal form
		Plane3D(Real alpha, Real beta, Real gamma, Real d)
		{
			_A = cos(alpha);
			_B = cos(beta);
			_C = cos(gamma);
			_D = -d;
		}
		// segments on coordinate axes
		Plane3D(Real seg_x, Real seg_y, Real seg_z)
		{
			if (seg_x == 0 || seg_y == 0 || seg_z == 0)
				throw std::runtime_error("Plane3D ctor - zero segment");

			_A = 1 / seg_x;
			_B = 1 / seg_y;
			_C = 1 / seg_z;
			_D = -1;
		}

		static Plane3D GetXYPlane() { return Plane3D(Pnt3Cart(0, 0, 0), Vec3Cart(0, 0, 1)); }
		static Plane3D GetXZPlane() { return Plane3D(Pnt3Cart(0, 0, 0), Vec3Cart(0, 1, 0)); }
		static Plane3D GetYZPlane() { return Plane3D(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0)); }

		Real  A() const { return _A; }
		Real& A() { return _A; }
		Real  B() const { return _B; }
		Real& B() { return _B; }
		Real  C() const { return _C; }
		Real& C() { return _C; }
		Real  D() const { return _D; }
		Real& D() { return _D; }

		Vec3Cart	Normal() const { return Vec3Cart(_A, _B, _C); }
		Pnt3Cart	GetPointOnPlane() const {
			if (_A != 0.0)
				return Pnt3Cart(-_D / _A, 0, 0);
			else  if (_B != 0.0)
				return Pnt3Cart(0, -_D / _B, 0);
			else
				return Pnt3Cart(0, 0, -_D / _C);
		}

		void GetCoordAxisSegments(Real& outseg_x, Real& outseg_y, Real& outseg_z)
		{
			if (_A != 0.0)
				outseg_x = -_D / _A;
			else
				outseg_x = Constants::PosInf;

			if (_B != 0.0)
				outseg_y = -_D / _B;
			else
				outseg_y = Constants::PosInf;

			if (_C != 0.0)
				outseg_z = -_D / _C;
			else
				outseg_z = Constants::PosInf;
		}

		// point to plane operations
		bool IsPointOnPlane(const Pnt3Cart& pnt, Real defEps = Defaults::Plane3DIsPointOnPlaneTolerance) const
		{
			return std::abs(_A * pnt.X() + _B * pnt.Y() + _C * pnt.Z() + _D) < defEps;
		}
		Real DistToPoint(const Pnt3Cart& pnt) const
		{
			Real a = _A * pnt.X() + _B * pnt.Y() + _C * pnt.Z() + _D;
			Real b = sqrt(_A * _A + _B * _B + _C * _C);

			return std::abs(a / b);
		}

		Pnt3Cart ProjectionToPlane(const Pnt3Cart& pnt) const
		{
			if (IsPointOnPlane(pnt))
				return pnt;

			Line3D line(pnt, Normal());

			// find intersection of this line with plane
			Pnt3Cart ret;
			IntersectionWithLine(line, ret);

			return ret;
		}

		// line to plane operations
		bool IsLineOnPlane(const Line3D& line) const
		{
			// get two points
			const Pnt3Cart pnt1 = line.StartPoint();
			const Pnt3Cart pnt2 = line(1.0);

			if (IsPointOnPlane(pnt1) && IsPointOnPlane(pnt2))
				return true;
			else
				return false;
		}
		Real AngleToLine(const Line3D& line) const
		{
			// angle between line and normal to plane
			return Constants::PI / 2.0 - line.Direction().AngleToVector(this->Normal());
		}
		// TODO - check
		bool IntersectionWithLine(const Line3D& line, Pnt3Cart& out_inter_pnt) const
		{
			Vec3Cart line_dir = line.Direction();
			Vec3Cart plane_normal = Normal();

			// Check if the line is parallel to the plane
			if (ScalarProduct(line_dir, plane_normal) == 0) {
				return false;
			}

			// Calculate the distance between the point on the line and the plane
			double dist = ScalarProduct(Vec3Cart(GetPointOnPlane(), line.StartPoint()), plane_normal) / ScalarProduct(line_dir, plane_normal);

			// Calculate the point of intersection
			Pnt3Cart inter_pnt = line.StartPoint() + line_dir * dist;

			// Set the output parameter to the point of intersection
			out_inter_pnt = inter_pnt;

			return true;
		}

		// plane to plane operations
		bool IsParallelToPlane(const Plane3D& plane) const
		{
			Vec3Cart norm1(_A, _B, _C);
			Vec3Cart norm2(plane._A, plane._B, plane._C);

			return norm1.IsParallelTo(norm2);
		}
		bool IsPerpendicularToPlane(const Plane3D& plane) const
		{
			Vec3Cart norm1(_A, _B, _C);
			Vec3Cart norm2(plane._A, plane._B, plane._C);

			return norm1.IsPerpendicularTo(norm2);
		}
		Real AngleToPlane(const Plane3D& plane) const
		{
			// angle between normals of two planes
			return this->Normal().AngleToVector(plane.Normal());
		}
		Real DistToPlane(const Plane3D& plane) const
		{
			if (IsParallelToPlane(plane))
			{
				// distance between two parallel planes is the distance between any point on one plane and the other plane
				return DistToPoint(plane.GetPointOnPlane());
			}
			else
				return 0.0;
		}
		// TODO - check
		bool IntersectionWithPlane(const Plane3D& plane, Line3D& out_inter_line) const
		{
			Vec3Cart inter_dir = VectorProduct(Normal(), plane.Normal());

			// Check if the planes are parallel
			if (inter_dir.NormL2() == 0.0) {
				return false;
			}

			// Calculate a point on the intersection line
			Pnt3Cart inter_pnt = GetPointOnPlane();

			// Calculate the distance between the intersection line and the two planes
			double dist1 = ScalarProduct(Vec3Cart(inter_pnt, GetPointOnPlane()), plane.Normal());
			double dist2 = ScalarProduct(Vec3Cart(inter_pnt, plane.GetPointOnPlane()), Normal());

			// Calculate the point of intersection
			inter_pnt = inter_pnt - inter_dir * (dist1 / ScalarProduct(inter_dir, plane.Normal()));

			// Set the output parameter to the intersection line
			out_inter_line = Line3D(inter_pnt, inter_dir);

			return true;
			return false;
		}
	};

	class Triangle3D
	{
	protected:
		Pnt3Cart _pnt1, _pnt2, _pnt3;
	public:
		Triangle3D() {}
		Triangle3D(Pnt3Cart pnt1, Pnt3Cart pnt2, Pnt3Cart pnt3)
			: _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3)
		{
		}

		Real A() const { return _pnt1.Dist(_pnt2); }
		Real B() const { return _pnt2.Dist(_pnt3); }
		Real C() const { return _pnt3.Dist(_pnt1); }

		Pnt3Cart	Pnt1() const { return _pnt1; }
		Pnt3Cart& Pnt1()			 { return _pnt1; }
		Pnt3Cart	Pnt2() const { return _pnt2; }
		Pnt3Cart& Pnt2()			 { return _pnt2; }
		Pnt3Cart	Pnt3() const { return _pnt3; }
		Pnt3Cart& Pnt3()			 { return _pnt3; }

		Real Area() const
		{
			Real s = (A() + B() + C()) / 2.0;
			return sqrt(s * (s - A()) * (s - B()) * (s - C()));
		}

		bool IsRight() const
		{
			return (hypot(A(), B()) == C() || hypot(A(), C()) == B() || hypot(B(), C()) == A());
		}
		bool IsIsosceles() const		// two sides are the same length
		{
			return (A() == B()) || (A() == C()) || (B() == C());
		}
		bool IsEquilateral() const	// all sides are the same length
		{
			return (A() == B()) && (A() == C());
		}

		Plane3D getDefinedPlane() const
		{
			return Plane3D(_pnt1, _pnt2, _pnt3);
		}
	};

	// Represents triangular surface in 3D and defines all needed mappings 
	// to represent Triangle3D as IParametricSurface
	class TriangleSurface3D : public Triangle3D, IParametricSurface<3>
	{
	public:
		Real _minX, _maxX, _minY, _maxY;
		Pnt3Cart _origin;
		Pnt3Cart _center;
		Vec3Cart _localX, _localY;
		Vec3Cart _normal;
		Real _pnt3XCoord;

		// pnt1-pnt2 should be hypothenuse of the triangle!!!
		// but we will handle it, if it is not
		TriangleSurface3D(Pnt3Cart pnt1, Pnt3Cart pnt2, Pnt3Cart pnt3)
		{
			// CHECK THAT 1-2 side is THE LONGEST ONE!!!
			Real a = pnt1.Dist(pnt2);
			Real b = pnt2.Dist(pnt3);
			Real c = pnt3.Dist(pnt1);

			if (a >= b && a >= c)
			{
				; // nice, we are all good
			}
			else if (b >= a && b >= c)
			{
				// rotate points one place, so that 'b' side (pnt2-pnt3) is at the beginning
				Pnt3Cart tmp = pnt1;
				pnt1 = pnt2;
				pnt2 = pnt3;
				pnt3 = tmp;
			}
			else
			{
				// rotate points two places
				Pnt3Cart tmp = pnt1;
				pnt1 = pnt3;
				pnt2 = tmp;
				pnt3 = pnt2;
			}

			_pnt1 = pnt1;		// initializing points in base Triangle3D class
			_pnt2 = pnt2;
			_pnt3 = pnt3;

			// calculate min and max
			// calculate center
			_center = (pnt1 + pnt2 + pnt3) / 3.0;

			// calculate local coordinate system
			// we will take x-axis to be from pnt1 to pnt2
			_localX = Vec3Cart(pnt1, pnt2).GetAsUnitVector();

			// we will calculate y-axis as following
			// calculate perpendicular vector to x-axis, that goes through pnt3
			Line3D	 lineX(pnt1, pnt2);
			Pnt3Cart pntOnLine = lineX.NearestPointOnLine(pnt3);

			_localY = Vec3Cart(pntOnLine, pnt3).GetAsUnitVector();

			_minX = 0.0;
			_maxX = pnt1.Dist(pnt2);
			_minY = 0.0;
			_maxY = pntOnLine.Dist(pnt3);
			_pnt3XCoord = pntOnLine.Dist(pnt1);

			_normal = VectorProduct(_localX, _localY).GetAsUnitVector();
		}

		virtual Real getMinU() const { return _minX; }
		virtual Real getMaxU() const { return _maxX; }
		virtual Real getMinW(Real u) const { return _minY; }
		virtual Real getMaxW(Real u) const
		{
			// this depends on value of u (which is localX)
			if (u < _pnt3XCoord)
				return _minY + (u / _pnt3XCoord) * (_maxY - _minY);
			else
				return _minY + ((_maxX - u) / (_maxX - _pnt3XCoord)) * (_maxY - _minY);
		}

		// for a given (u,w), which is basically (x,y) in local coordinate system
		// return point in global coordinate system
		VectorN<Real, 3> operator()(Real u, Real w) const
		{
			Pnt3Cart ret = _origin + u * _localX + w * _localY;
			return VectorN<Real, 3>({ ret.X(), ret.Y(), ret.Z() });
		}
	};

	// represents rectangular surface in 3D
	class RectSurface3D : public IParametricSurfaceRect<3>
	{
	public:
		// all constructors are expected to set properly these values
		Pnt3Cart _pnt1, _pnt2, _pnt3, _pnt4;
		Real _minX, _maxX, _minY, _maxY;
		Pnt3Cart _center;
		Vec3Cart _localX, _localY;

		RectSurface3D() {}
		// MUST SET NORMAL, FOR ORIENTATION! (uzeti right-hand ... defined by order of points?)
		// zadati i centralnom tockom, vektorom normale, uz dodatni a i b!
		RectSurface3D(Pnt3Cart pnt1, Pnt3Cart pnt2, Pnt3Cart pnt3, Pnt3Cart pnt4)
			: _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3), _pnt4(pnt4)
		{
			// KLJUCNO - provjeriti da li su sve u ravnini!

			// podesiti min i max
			Real lenX = _pnt1.Dist(_pnt2);
			Real lenY = _pnt1.Dist(_pnt4);
			_minX = -lenX / 2;
			_maxX = lenX / 2;
			_minY = -lenY / 2;
			_maxY = lenY / 2;

			// calculate center
			_center = (_pnt1 + _pnt2 + _pnt3 + _pnt4) / 4.0;

			// calculate local coordinate system
			_localX = Vec3Cart(_pnt1, _pnt2).GetAsUnitVector();
			_localY = Vec3Cart(_pnt1, _pnt4).GetAsUnitVector();
		}

		virtual Real getMinU() const { return _minX; }
		virtual Real getMaxU() const { return _maxX; }
		virtual Real getMinW() const { return _minY; }
		virtual Real getMaxW() const { return _maxY; }

		Vec3Cart getNormal() const {
			return VectorProduct(Vec3Cart(_pnt1, _pnt2), Vec3Cart(_pnt1, _pnt4)).GetAsUnitVector();
		}
		Pnt3Cart getCenter() const { return _center; }

		Real getArea() const {
			return VectorProduct(Vec3Cart(_pnt1, _pnt2), Vec3Cart(_pnt1, _pnt4)).NormL2();
		}

		// vraca dva Triangle3D - i orijentacija je parametar! (kako ce odabrati tocke)

		VectorN<Real, 3> operator()(Real u, Real w) const {
			Pnt3Cart ret = _center + u * _localX + w * _localY;
			return VectorN<Real, 3>({ ret.X(), ret.Y(), ret.Z() });
		}
	};

	class IBody
	{
	public:
		virtual bool isInside(const Pnt3Cart& pnt) const = 0;
	};

	class ISolidBodyWithBoundary : public IBody
	{
	public:
		Real _x1, _x2;

		Real(*_y1)(Real);
		Real(*_y2)(Real);

		Real(*_z1)(Real, Real);
		Real(*_z2)(Real, Real);

	public:
		ISolidBodyWithBoundary(Real x1, Real x2,
			Real(*y1)(Real), Real(*y2)(Real),
			Real(*z1)(Real, Real), Real(*z2)(Real, Real))
			: _x1(x1), _x2(x2), _y1(y1), _y2(y2), _z1(z1), _z2(z2)
		{	}

		virtual Real getDensity(const VectorN<Real, 3>& x) const = 0;

		virtual bool isInside(const Pnt3Cart& pnt) const
		{
			const Real x = pnt.X();
			const Real y = pnt.Y();
			const Real z = pnt.Z();

			if (_x1 < x && x < _x2)
			{
				// check y bounds
				if (_y1(x) < y && y < _y2(x))
				{
					// check z bounds
					if (_z1(x, y) < z && z < _z2(x, y))
					{
						return true;	// point is inside the solid body
					}
				}
			}
			return false;
		}

	};

	class SolidBodyWithBoundary : public ISolidBodyWithBoundary
	{
		Real(*_density)(const VectorN<Real, 3>& x);
	public:
		SolidBodyWithBoundary(Real x1, Real x2, Real(*y1)(Real), Real(*y2)(Real),
			Real(*z1)(Real, Real), Real(*z2)(Real, Real),
			Real(*density)(const VectorN<Real, 3>& x))
			: ISolidBodyWithBoundary(x1, x2, y1, y2, z1, z2), _density(density)
		{	}

		virtual Real getDensity(const VectorN<Real, 3>& x) const
		{
			return _density(x);
		}
	};

	class SolidBodyWithBoundaryConstDensity : public ISolidBodyWithBoundary
	{
		Real _density;
	public:
		SolidBodyWithBoundaryConstDensity(Real x1, Real x2, Real(*y1)(Real), Real(*y2)(Real),
			Real(*z1)(Real, Real), Real(*z2)(Real, Real),
			Real density)
			: ISolidBodyWithBoundary(x1, x2, y1, y2, z1, z2), _density(density)
		{	}

		virtual Real getDensity(const VectorN<Real, 3>& x) const
		{
			return _density;
		}
	};

	// solid body in 3D defined by triangular surfaces
	class BodyWithTriangleSurfaces : public IBody
	{
	protected:
		std::vector<TriangleSurface3D> _surfaces;
	public:
		int getSurfaceCount() const
		{
			return _surfaces.size();
		}
		const TriangleSurface3D& getSurface(int index) const
		{
			if (index < 0 || index >= getSurfaceCount())
				throw std::out_of_range("BodyWithTriangleSurfaces::getSurface - index out of range");
			return _surfaces[index];
		}
	};

	// solid body in 3D defined by rectangular surfaces
	class BodyWithRectSurfaces : public IBody
	{
	protected:
		std::vector<RectSurface3D> _surfaces;

	public:
		int getSurfaceCount() const
		{
			return static_cast<int>(_surfaces.size());
		}
		const RectSurface3D& getSurface(int index) const
		{
			if (index < 0 || index >= getSurfaceCount())
				throw std::out_of_range("BodyWithRectSurfaces::getSurface - index out of range");
			return _surfaces[index];
		}

		// isClosed() - iz centra mase (?) odasilje zrake u svim smje
		// rovima i gleda da li je pogodio iti jednu povrsinu
		// vraca listu povrsina, koje bi TREBALE omedjivati tijelo!

		// i po tome se moze obaviti surface integracija!!!
	};

	// represents solid that is composed of other solids
	// dobar nacin za provjeriti je li composed solid "tight" je izracunati fluks kroz sve povrsine
	class ComposedSolidSurfaces3D
	{
		bool IsInside(const Pnt3Cart& pnt) const
		{
			// TODO - implement
			return false;
		}
	};

	class Cube3D : public BodyWithRectSurfaces
	{
		// kocka
		Real _a;
		Pnt3Cart _center;
	public:
		Cube3D(Real a) : _a(a)
		{
			Pnt3Cart pnt1( a / 2, -a / 2, -a / 2);
			Pnt3Cart pnt2( a / 2,  a / 2, -a / 2);
			Pnt3Cart pnt3(-a / 2,  a / 2, -a / 2);
			Pnt3Cart pnt4(-a / 2, -a / 2, -a / 2);
			Pnt3Cart pnt5( a / 2, -a / 2,  a / 2);
			Pnt3Cart pnt6( a / 2,  a / 2,  a / 2);
			Pnt3Cart pnt7(-a / 2,  a / 2,  a / 2);
			Pnt3Cart pnt8(-a / 2, -a / 2,  a / 2);

			// dodati svih 6 stranica u popis povrsina
			_surfaces.push_back(RectSurface3D(pnt1, pnt4, pnt3, pnt2));     // lower side in xy plane
			_surfaces.push_back(RectSurface3D(pnt5, pnt6, pnt7, pnt8));     // upper side in xy plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt2, pnt6, pnt5));     // front side in yz plane
			_surfaces.push_back(RectSurface3D(pnt4, pnt8, pnt7, pnt3));     // back side in yz plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt5, pnt8, pnt4));     // left side in xz plane
			_surfaces.push_back(RectSurface3D(pnt2, pnt3, pnt7, pnt6));     // right side in xz plane
		}
		Cube3D(Real a, const Pnt3Cart& center) : Cube3D(a)
		{
			_center = center;

			// create 8 points for corners of the cube, centered at the given center point
			Pnt3Cart pnt1(center.X() + a / 2, center.Y() - a / 2, center.Z() - a / 2);
			Pnt3Cart pnt2(center.X() + a / 2, center.Y() + a / 2, center.Z() - a / 2);
			Pnt3Cart pnt3(center.X() - a / 2, center.Y() + a / 2, center.Z() - a / 2);
			Pnt3Cart pnt4(center.X() - a / 2, center.Y() - a / 2, center.Z() - a / 2);
			Pnt3Cart pnt5(center.X() + a / 2, center.Y() - a / 2, center.Z() + a / 2);
			Pnt3Cart pnt6(center.X() + a / 2, center.Y() + a / 2, center.Z() + a / 2);
			Pnt3Cart pnt7(center.X() - a / 2, center.Y() + a / 2, center.Z() + a / 2);
			Pnt3Cart pnt8(center.X() - a / 2, center.Y() - a / 2, center.Z() + a / 2);

			// dodati svih 6 stranica u popis povrsina
			_surfaces.push_back(RectSurface3D(pnt1, pnt4, pnt3, pnt2));     // lower side in xy plane
			_surfaces.push_back(RectSurface3D(pnt5, pnt6, pnt7, pnt8));     // upper side in xy plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt2, pnt6, pnt5));     // front side in yz plane
			_surfaces.push_back(RectSurface3D(pnt4, pnt8, pnt7, pnt3));     // back side in yz plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt5, pnt8, pnt4));     // left side in xz plane
			_surfaces.push_back(RectSurface3D(pnt2, pnt3, pnt7, pnt6));     // right side in xz plane
		}

		bool isInside(const Pnt3Cart& pnt) const
		{
			// check if point is inside the cube
			// assuming that cube is centered at origin
			Real half_a = _a / 2.0;
			return (pnt.X() >= -half_a && pnt.X() <= half_a &&
							pnt.Y() >= -half_a && pnt.Y() <= half_a &&
							pnt.Z() >= -half_a && pnt.Z() <= half_a);
		}
	};

	class Pyramid3D : public BodyWithTriangleSurfaces
	{
		Real _a;			// base side length
		Real _h;			// pyramid height
		Vec3Cart _center;
	public:
		Pyramid3D(Real a, Real h) : _a(a), _h(h)
		{
			Pnt3Cart pnt1( a / 2, -a / 2, -a / 2);
			Pnt3Cart pnt2( a / 2,  a / 2, -a / 2);
			Pnt3Cart pnt3(-a / 2,  a / 2, -a / 2);
			Pnt3Cart pnt4(-a / 2, -a / 2, -a / 2);
			Pnt3Cart pnt5(     0,      0,      h);	// apex of the pyramid
			
			_surfaces.push_back(TriangleSurface3D(pnt1, pnt4, pnt3));     // lower side in xy plane
			_surfaces.push_back(TriangleSurface3D(pnt5, pnt1, pnt2));     // front side in yz plane
			_surfaces.push_back(TriangleSurface3D(pnt5, pnt2, pnt3));     // right side in yz plane
			_surfaces.push_back(TriangleSurface3D(pnt5, pnt3, pnt4));     // back side in yz plane
			_surfaces.push_back(TriangleSurface3D(pnt5, pnt4, pnt1));     // left side in xz plane
		}
		Pyramid3D(Real a, Real h, const Vec3Cart& center) : Pyramid3D(a, h)
		{
			_center = center;
		}

		bool isInside(const Pnt3Cart& pnt) const
		{
			// check if point is inside the pyramid
			// assuming that pyramid is centered at origin
			Real half_a = _a / 2.0;
			if (pnt.Z() < 0 || pnt.Z() > _a / sqrt(3))
				return false;	// point is outside the height of the pyramid
			
			// TODO - finish
			return false;
		}
	};

	class PyramidEquilateral3D : public Pyramid3D
	{
	public:
		PyramidEquilateral3D(Real a) : Pyramid3D(a, a / sqrt(3))
		{ }
	};
}


#endif
