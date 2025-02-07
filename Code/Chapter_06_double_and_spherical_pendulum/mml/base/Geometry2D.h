#if !defined MML_GEOMETRY_2D_H
#define MML_GEOMETRY_2D_H

#include "MMLBase.h"

#include "base/VectorTypes.h"
#include "base/Geometry.h"

namespace MML
{
	class Line2D
	{
	private:
		Point2Cartesian _point;
		Vector2Cartesian _direction; // unit vector in line direction

	public:
		Line2D(const Point2Cartesian& pnt, const Vector2Cartesian dir)
		{
			_point = pnt;
			_direction = dir.GetUnitVector();
		}

		Line2D(const Point2Cartesian& a, const Point2Cartesian& b)
		{
			Vector2Cartesian dir(a, b);
			_point = a;
			_direction = dir.GetUnitVector();
		}

		Point2Cartesian   StartPoint() const { return _point; }
		Point2Cartesian&  StartPoint() { return _point; }

		Vector2Cartesian  Direction() const { return _direction; }
		Vector2Cartesian& Direction() { return _direction; }

		Point2Cartesian operator()(Real t) const
		{
			Vector2Cartesian dist = t * _direction;
			Point2Cartesian ret = _point + dist;
			return ret;
		}
	};

	class SegmentLine2D
	{
	private:
		Point2Cartesian _point1;
		Point2Cartesian _point2;

	public:
		SegmentLine2D(Point2Cartesian pnt1, Point2Cartesian pnt2) : _point1(pnt1), _point2(pnt2)
		{ }

		SegmentLine2D(const Point2Cartesian& pnt1, const Vector2Cartesian& direction, Real t) : _point1(pnt1)
		{
			_point2 = pnt1 + direction * t;
		}

		Point2Cartesian  StartPoint() const { return _point1; }
		Point2Cartesian& StartPoint() { return _point1; }

		Point2Cartesian  EndPoint()  const { return _point2; }
		Point2Cartesian& EndPoint() { return _point2; }

		Point2Cartesian PointOnSegment(Real t)
		{
			if (t < 0.0 || t > 1.0)
				throw std::invalid_argument("SegmentLine2D::PointOnSegment t must be in [0,1]");

			Vector2Cartesian dist = t * Direction();
			Point2Cartesian ret = _point1 + dist;
			return ret;
		}

		Real                Length()    const { return _point1.Dist(_point2); }
		Vector2Cartesian    Direction() const { return Vector2Cartesian(_point1, _point2); }
	};

	class Triangle2D
	{
	private:
		Point2Cartesian _pnt1, _pnt2, _pnt3;

	public:
		Triangle2D(Point2Cartesian pnt1, Point2Cartesian pnt2, Point2Cartesian pnt3) : _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3)
		{ }

		Point2Cartesian  Pnt1() const { return _pnt1; }
		Point2Cartesian& Pnt1() { return _pnt1; }
		Point2Cartesian  Pnt2() const { return _pnt2; }
		Point2Cartesian& Pnt2() { return _pnt2; }
		Point2Cartesian  Pnt3() const { return _pnt3; }
		Point2Cartesian& Pnt3() { return _pnt3; }

		Real Area() const
		{
			Real a = _pnt1.Dist(_pnt2);
			Real b = _pnt2.Dist(_pnt3);
			Real c = _pnt3.Dist(_pnt1);

			Real s = (a + b + c) / 2.0;

			return sqrt(s * (s - a) * (s - b) * (s - c));
		}
	};

	class Polygon2D
	{
	private:
		std::vector<Point2Cartesian> _points;
	public:
		Polygon2D() {}
		Polygon2D(std::vector<Point2Cartesian> points) : _points(points) {}
		Polygon2D(std::initializer_list<Point2Cartesian> list)
		{
			for (auto element : list)
				_points.push_back(element);
		}

		std::vector<Point2Cartesian>  Points() const { return _points; }
		std::vector<Point2Cartesian>& Points() { return _points; }

		Real Area() const
		{
			Real area = 0.0;
			int n = (int)_points.size();
			for (int i = 0; i < n; i++)
			{
				area += _points[i].X() * _points[(i + 1) % n].Y();
				area -= _points[i].Y() * _points[(i + 1) % n].X();
			}
			area /= 2.0;
			return area;
		}

		bool IsSimple() const
		{
			return false;
		}

		bool IsConvex() const
		{
			int n = (int)_points.size();
			if (n < 3)
				return false;
			bool sign = false;
			for (int i = 0; i < n; i++)
			{
				Vector2Cartesian v1(_points[(i + 1) % n], _points[i]);
				Vector2Cartesian v2(_points[(i + 2) % n], _points[(i + 1) % n]);
				Real cross = v1.X() * v2.Y() - v1.Y() * v2.X();
				if (i == 0)
					sign = cross > 0;
				else if (sign != (cross > 0))
					return false;
			}
			return true;
		}

		std::vector<Triangle2D> Triangularization() const
		{
			std::vector<Triangle2D> triangles;
			int n = (int)_points.size();
			if (n < 3)
				return triangles;
			std::vector<int> indices(n);
			for (int i = 0; i < n; i++)
				indices[i] = i;
			int i = 0;
			while (n > 3)
			{
				int i1 = indices[i];
				int i2 = indices[(i + 1) % n];
				int i3 = indices[(i + 2) % n];
				Triangle2D tri(_points[i1], _points[i2], _points[i3]);
				triangles.push_back(tri);
				indices.erase(indices.begin() + (i + 1) % n);
				n--;
				i = (i + 1) % n;
			}
			Triangle2D tri(_points[indices[0]], _points[indices[1]], _points[indices[2]]);
			triangles.push_back(tri);
			return triangles;
		}

		bool IsInside(Point2Cartesian pnt) const
		{
			int n = (int)_points.size();
			if (n < 3)
				return false;
		
			return false;
		}
	};
}
#endif