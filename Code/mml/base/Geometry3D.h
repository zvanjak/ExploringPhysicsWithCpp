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
		Point3Cartesian  _point;
		Vector3Cartesian _direction;

	public:
		Line3D() {}
		// by default, direction vector is normalized to unit vector (but, it need not be such!)
		Line3D(const Point3Cartesian& pnt, const Vector3Cartesian dir)
		{
			// check for null vector as direction
			if (dir.X() == 0.0 && dir.Y() == 0.0 && dir.Z() == 0.0)
				throw std::runtime_error("Line3D ctor - null vector as direction");

			_point = pnt;
			_direction = dir.GetAsUnitVector();
		}
		Line3D(const Point3Cartesian& a, const Point3Cartesian& b)
		{
			// check for same points
			if (a == b)
				throw std::runtime_error("Line3D ctor - same points");

			Vector3Cartesian dir(a, b);
			_point = a;
			_direction = dir.GetAsUnitVector();
		}

		Point3Cartesian   StartPoint() const { return _point; }
		Point3Cartesian&  StartPoint() { return _point; }

		Vector3Cartesian  Direction() const { return _direction; }
		Vector3Cartesian& Direction() { return _direction; }

		Point3Cartesian operator()(Real t) const { return _point + t * _direction; }

		bool IsPerpendicular(const Line3D& b) const
		{
			return ScalarProduct(Direction(), b.Direction()) == 0.0f;
		}
		bool IsParallel(const Line3D& b) const
		{
			return Direction() == b.Direction();
		}

		// distance between point and line
		Real Dist(const Point3Cartesian& pnt) const
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
		// distance between two lines
		Real Dist(const Line3D& line) const
		{
			// https://math.stackexchange.com/questions/2213165/distance-between-two-lines-in-3d-space
			// https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
			Point3Cartesian  p1 = StartPoint();
			Vector3Cartesian d1 = Direction();
			Point3Cartesian  p2 = line.StartPoint();
			Vector3Cartesian d2 = line.Direction();

			Vector3Cartesian n = VectorProduct(d1, d2);

			if (n.IsNullVec())
			{
				// parallel lines
				return Vector3Cartesian(p1, p2) * d1 / d2.NormL2();
			}
			else
			{
				// skew lines
				return Abs(n * Vector3Cartesian(p1, p2)) / n.NormL2();
			}
		}
		// distance between two lines, while also returning nearest points on both lines
		Real Dist(const Line3D& line, Point3Cartesian& out_line1_pnt, Point3Cartesian& out_line2_pnt) const
		{
			// https://math.stackexchange.com/questions/2213165/distance-between-two-lines-in-3d-space
			// https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
			Point3Cartesian  p1 = StartPoint();
			Vector3Cartesian d1 = Direction();
			Point3Cartesian  p2 = line.StartPoint();
			Vector3Cartesian d2 = line.Direction();
			Real dist = 0.0;

			Vector3Cartesian n = VectorProduct(d1, d2);

			if (n.IsNullVec())
			{
				// parallel lines
				dist = Vector3Cartesian(p1, p2) * d1 / d2.NormL2();
				// should we throw exception???
				// because we can't get out_line_pnt's
			}
			else
			{
				// skew lines
				dist = Abs(n * Vector3Cartesian(p1, p2)) / n.NormL2();

				Real t1 = (VectorProduct(d2, n) * Vector3Cartesian(p1, p2)) / POW2(n.NormL2());
				Real t2 = (VectorProduct(d1, n) * Vector3Cartesian(p1, p2)) / POW2(n.NormL2());

				out_line1_pnt = p1 + t1 * d1;
				out_line2_pnt = p2 + t2 * d2;
			}

			return dist;
		}

		// nearest point on line to given point
		Point3Cartesian NearestPointOnLine(const Point3Cartesian& pnt) const
		{
			// https://math.stackexchange.com/questions/1521128/given-a-line-and-a-point-in-3d-how-to-find-the-closest-point-on-the-line         
			Vector3Cartesian line_dir = this->Direction();
			Vector3Cartesian rad_vec_AP(StartPoint(), pnt);

			double t = rad_vec_AP * line_dir / POW2(line_dir.NormL2());

			return StartPoint() + t * line_dir;
		}

		// intersection of two lines
		bool Intersection(const Line3D& line, Point3Cartesian& out_inter_pnt) const
		{
			// https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
			// https://math.stackexchange.com/questions/2213165/distance-between-two-lines-in-3d-space

			return true;
		}

		// perpendicular line that goes through givenpoint
		Line3D PerpendicularLineThroughPoint(const Point3Cartesian& pnt)
		{
			Line3D ret;

			return ret;
		}
	};

	class SegmentLine3D
	{
	private:
		Point3Cartesian _point1;
		Point3Cartesian _point2;

	public:
		SegmentLine3D(Point3Cartesian pnt1, Point3Cartesian pnt2) : _point1(pnt1), _point2(pnt2)
		{
		}
		SegmentLine3D(Point3Cartesian pnt1, Vector3Cartesian direction, Real t)
		{
			_point1 = pnt1;
			_point2 = pnt1 + t * direction;
		}

		Point3Cartesian   StartPoint() const { return _point1; }
		Point3Cartesian&  StartPoint() { return _point1; }

		Point3Cartesian   EndPoint() const { return _point2; }
		Point3Cartesian&  EndPoint() { return _point2; }

		Point3Cartesian		PointOnSegment(Real t)
		{
			if (t < 0.0 || t > 1.0)
				throw std::runtime_error("SegmentLine3D::PointOnSegment - t not in [0,1]");

			return _point1 + t * Direction();
		}

		Real              Length()    const { return _point1.Dist(_point2); }
		Vector3Cartesian  Direction() const { return Vector3Cartesian(_point1, _point2); }
	};

	class Plane3D
	{
	private:
		Real _A, _B, _C, _D;

	public:
		Plane3D(const Point3Cartesian& a, const Vector3Cartesian& normal)
		{
			if (normal.IsNullVec())
				throw std::runtime_error("Plane3D ctor - normal is null vector");

			Vector3Cartesian unitNormal = normal.GetAsUnitVector();

			_A = unitNormal.X();
			_B = unitNormal.Y();
			_C = unitNormal.Z();
			_D = -(a.X() * unitNormal.X() + a.Y() * unitNormal.Y() + a.Z() * unitNormal.Z());
		}
		Plane3D(const Point3Cartesian& a, const Point3Cartesian& b, const Point3Cartesian& c)
			: Plane3D(a, VectorProduct(Vector3Cartesian(a, b), Vector3Cartesian(a, c)))
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

		static Plane3D GetXYPlane() { return Plane3D(Point3Cartesian(0, 0, 0), Vector3Cartesian(0, 0, 1)); }
		static Plane3D GetXZPlane() { return Plane3D(Point3Cartesian(0, 0, 0), Vector3Cartesian(0, 1, 0)); }
		static Plane3D GetYZPlane() { return Plane3D(Point3Cartesian(0, 0, 0), Vector3Cartesian(1, 0, 0)); }

		Real  A() const { return _A; }
		Real& A()				{ return _A; }
		Real  B() const { return _B; }
		Real& B()				{ return _B; }
		Real  C() const { return _C; }
		Real& C()				{ return _C; }
		Real  D() const { return _D; }
		Real& D()				{ return _D; }

		Vector3Cartesian	Normal() const { return Vector3Cartesian(_A, _B, _C); }
		Point3Cartesian		GetPointOnPlane() const {
			if (_A != 0.0)
				return Point3Cartesian(-_D / _A, 0, 0);
			else  if (_B != 0.0)
				return Point3Cartesian(0, -_D / _B, 0);
			else
				return Point3Cartesian(0, 0, -_D / _C);
		}

		void GetCoordAxisSegments(Real& outseg_x, Real& outseg_y, Real& outseg_z)
		{
			if (_A != 0.0)
				outseg_x = -_D / _A;
			else
				outseg_x = Constants::PositiveInf;

			if (_B != 0.0)
				outseg_y = -_D / _B;
			else
				outseg_y = Constants::PositiveInf;

			if (_C != 0.0)
				outseg_z = -_D / _C;
			else
				outseg_z = Constants::PositiveInf;
		}

		// point to plane operations
		bool IsPointOnPlane(const Point3Cartesian& pnt, Real defEps = 1e-15) const
		{
			return std::abs(_A * pnt.X() + _B * pnt.Y() + _C * pnt.Z() + _D) < defEps;
		}
		Real DistToPoint(const Point3Cartesian& pnt) const
		{
			Real a = _A * pnt.X() + _B * pnt.Y() + _C * pnt.Z() + _D;
			Real b = sqrt(_A * _A + _B * _B + _C * _C);

			return std::abs(a / b);
		}

		Point3Cartesian ProjectionToPlane(const Point3Cartesian& pnt) const
		{
			if (IsPointOnPlane(pnt))
				return pnt;

			Line3D line(pnt, Normal());

			// find intersection of this line with plane
			Point3Cartesian ret;
			IntersectionWithLine(line, ret);

			return ret;
		}

		// line to plane operations
		bool IsLineOnPlane(const Line3D& line) const
		{
			// get two points
			const Point3Cartesian pnt1 = line.StartPoint();
			const Point3Cartesian pnt2 = line(1.0);

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
		bool IntersectionWithLine(const Line3D& line, Point3Cartesian& out_inter_pnt) const
		{
			// Calculate the direction vector of the line
			Vector3Cartesian line_dir = line.Direction();

			// Calculate the normal vector of the plane
			Vector3Cartesian plane_normal = Normal();

			// Check if the line is parallel to the plane
			if (ScalarProduct(line_dir, plane_normal) == 0) {
				return false;
			}

			// Calculate the distance between the point on the line and the plane
			double dist = ScalarProduct(Vector3Cartesian(GetPointOnPlane(), line.StartPoint()), plane_normal) / ScalarProduct(line_dir, plane_normal);

			// Calculate the point of intersection
			Point3Cartesian inter_pnt = line.StartPoint() + line_dir * dist;

			// Set the output parameter to the point of intersection
			out_inter_pnt = inter_pnt;

			return true;
		}

		// plane to plane operations
		bool IsParallelToPlane(const Plane3D& plane) const
		{
			Vector3Cartesian norm1(_A, _B, _C);

			Vector3Cartesian norm2(plane._A, plane._B, plane._C);

			return norm1.IsParallelTo(norm2);
		}
		bool IsPerpendicularToPlane(const Plane3D& plane) const
		{
			Vector3Cartesian norm1(_A, _B, _C);
			Vector3Cartesian norm2(plane._A, plane._B, plane._C);

			return norm1.IsPerpendicularTo(norm2);
		}
		Real AngleToPlane(const Plane3D& plane) const
		{
			// to je kut normala
			return this->Normal().AngleToVector(plane.Normal());
		}
		Real DistToPlane(const Plane3D& plane) const
		{
			// TODO finish
			// ili su paralelne, pa imamo neki broj, ili se sijeku pa je 0
			return 0.0;
		}
		bool IntersectionWithPlane(const Plane3D& plane, Line3D& out_inter_line) const
		{
			Vector3Cartesian inter_dir = VectorProduct(Normal(), plane.Normal());

			// Check if the planes are parallel
			if (inter_dir.NormL2() == 0.0) {
				return false;
			}

			// Calculate a point on the intersection line
			Point3Cartesian inter_pnt = GetPointOnPlane();

			// Calculate the distance between the intersection line and the two planes
			double dist1 = ScalarProduct(Vector3Cartesian(inter_pnt, GetPointOnPlane()), plane.Normal());
			double dist2 = ScalarProduct(Vector3Cartesian(inter_pnt, plane.GetPointOnPlane()), Normal());

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
	private:
		Point3Cartesian _pnt1, _pnt2, _pnt3;
	public:
		Triangle3D(Point3Cartesian pnt1, Point3Cartesian pnt2, Point3Cartesian pnt3)
			: _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3)
		{
		}

		Real A() const { return _pnt1.Dist(_pnt2); }
		Real B() const { return _pnt2.Dist(_pnt3); }
		Real C() const { return _pnt3.Dist(_pnt1); }
		
		Point3Cartesian		Pnt1() const  { return _pnt1; }
		Point3Cartesian&	Pnt1()				{ return _pnt1; }
		Point3Cartesian		Pnt2() const  { return _pnt2; }
		Point3Cartesian&	Pnt2()				{ return _pnt2; }
		Point3Cartesian		Pnt3() const  { return _pnt3; }
		Point3Cartesian&	Pnt3()				{ return _pnt3; }
	};

	class TriangleSurface3D : public Triangle3D, IParametricSurfaceRect<3>
	{
	public:
		Real _minX, _maxX, _minY, _maxY;
		Point3Cartesian _center;
		Vector3Cartesian _localX, _localY;

		TriangleSurface3D(Point3Cartesian pnt1, Point3Cartesian pnt2, Point3Cartesian pnt3)
			: Triangle3D(pnt1, pnt2, pnt3)
		{
			// calculate min and max
			_minX = std::min(pnt1.X(), std::min(pnt2.X(), pnt3.X()));
			_maxX = std::max(pnt1.X(), std::max(pnt2.X(), pnt3.X()));
			_minY = std::min(pnt1.Y(), std::min(pnt2.Y(), pnt3.Y()));
			_maxY = std::max(pnt1.Y(), std::max(pnt2.Y(), pnt3.Y()));

			// calculate center
			_center = (pnt1 + pnt2 + pnt3) / 3.0;

			// calculate local coordinate system
			_localX = Vector3Cartesian(pnt1, pnt2).GetAsUnitVector();
			_localY = Vector3Cartesian(pnt1, pnt3).GetAsUnitVector();
		}

		virtual Real getMinX() const { return _minX; }
		virtual Real getMaxX() const { return _maxX; }
		virtual Real getMinY() const { return _minY; }
		virtual Real getMaxY() const { return _maxY; }
	};

	class RectSurface3D : public IParametricSurfaceRect<3>
	{
	public:
		// predstavlja PRAVOKUTNU povrsinu u 3D
		Point3Cartesian _pnt1, _pnt2, _pnt3, _pnt4;
		Real _minX, _maxX, _minY, _maxY;
		Point3Cartesian _center;
		Vector3Cartesian _localX, _localY;

		RectSurface3D() {}
		// zadati i centralnom tockom, vektorom normale, uz dodatni a i b!
		RectSurface3D(Point3Cartesian pnt1, Point3Cartesian pnt2, Point3Cartesian pnt3, Point3Cartesian pnt4)
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
			_localX = Vector3Cartesian(_pnt1, _pnt2).GetAsUnitVector();
			_localY = Vector3Cartesian(_pnt1, _pnt4).GetAsUnitVector();
		}
		virtual Real getMinU() const { return _minX; }
		virtual Real getMaxU() const { return _maxX; }
		virtual Real getMinW() const { return _minY; }
		virtual Real getMaxW() const { return _maxY; }

		Vector3Cartesian getNormal() const {
			return VectorProduct(Vector3Cartesian(_pnt1, _pnt2), Vector3Cartesian(_pnt1, _pnt4)).GetAsUnitVector();
		}
		Point3Cartesian getCenter() const { return _center; }

		Real getArea() const {
			return VectorProduct(Vector3Cartesian(_pnt1, _pnt2), Vector3Cartesian(_pnt1, _pnt4)).NormL2();
		}

		// vraca dva Triangle3D - i orijentacija je parametar! (kako ce odabrati tocke)

		VectorN<Real, 3> operator()(Real u, Real w) const {
			Point3Cartesian ret = _center + u * _localX + w * _localY;
			return VectorN<Real, 3>({ ret.X(), ret.Y(), ret.Z() });
		}
	};

	class IBody
	{
	public:
		virtual bool isInside(const Point3Cartesian& pnt) const = 0;

		// must also have boundary defined

	};

	class BodyWithRectSurfaces : public IBody
	{

	};

	class BodyWithTriangleSurfaces : public IBody
	{

	};

	class BodyWithBoundary : public IBody
	{

	};

	// TODO - IntegrableSolid? koji ima i potrebne funkcije kojima definira granice tijela?
	class IntegrableVolume3D
	{
		// osigurava da se znaju funkcije koje definiraju granice tijela
		// da se moze obaviti volume integracija
	};

	class SolidSurfaces3D : public BodyWithRectSurfaces
	{
	public:
		// solid body in 3D defined by surfaces
		std::vector<RectSurface3D> _surfaces;

	public:
		bool isInside(const Point3Cartesian& pnt) const
		{
			// TODO - implement
			return false;
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
		bool IsInside(const Point3Cartesian& pnt) const
		{
			// TODO - implement
			return false;
		}
	};

	class Cube3D : public SolidSurfaces3D
	{
		// kocka
		Real _a;
		Vector3Cartesian _center;
	public:
		Cube3D(Real a) : _a(a)
		{
			Point3Cartesian pnt1(a / 2, -a / 2, -a / 2);
			Point3Cartesian pnt2(a / 2, a / 2, -a / 2);
			Point3Cartesian pnt3(-a / 2, a / 2, -a / 2);
			Point3Cartesian pnt4(-a / 2, -a / 2, -a / 2);
			Point3Cartesian pnt5(a / 2, -a / 2, a / 2);
			Point3Cartesian pnt6(a / 2, a / 2, a / 2);
			Point3Cartesian pnt7(-a / 2, a / 2, a / 2);
			Point3Cartesian pnt8(-a / 2, -a / 2, a / 2);

			// dodati svih 6 stranica u popis povrsina
			_surfaces.push_back(RectSurface3D(pnt1, pnt4, pnt3, pnt2));     // lower side in xy plane
			_surfaces.push_back(RectSurface3D(pnt5, pnt6, pnt7, pnt8));     // upper side in xy plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt2, pnt6, pnt5));     // front side in yz plane
			_surfaces.push_back(RectSurface3D(pnt4, pnt8, pnt7, pnt3));     // back side in yz plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt5, pnt8, pnt4));     // left side in xz plane
			_surfaces.push_back(RectSurface3D(pnt2, pnt3, pnt7, pnt6));     // right side in xz plane
		}
		Cube3D(Real a, const Vector3Cartesian& center) : Cube3D(a)
		{
			_center = center;
		}
	};
}

#endif
