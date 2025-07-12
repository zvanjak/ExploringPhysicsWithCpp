#if !defined MML_SURFACE_INTEGRATION_H
#define MML_SURFACE_INTEGRATION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/Geometry3D.h"

#include "core/Derivation.h"
#include "core/Integration.h"
#include "core/FieldOperations.h"

namespace MML
{
	class SurfaceIntegration
	{
	public:
		// povrsinski integral preko
		// - ParametricSurface
		// - ParametricSurfaceRect
		// - SolidSurfaces3D -> 
		//		- SurfacesWithTriangles
		// 		- SurfacesWithRects
		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const IParametricSurfaceRect<3>& surface, const Real x1, const Real x2, const Real y1, const Real y2)
		{
			return 0.0;
		}
		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const BodyWithRectSurfaces& solid)
		{
			Real total = 0.0;
			for (int i = 0; i < solid.getSurfaceCount(); i++)
				total += SurfaceIntegral(vectorField, solid.getSurface(i));

			return total;
		}
		static Real CalcSurfaceContrib(const IVectorFunction<3>& vectorField, const RectSurface3D& surface)
		{
			Real result = 0.0;
			Vector3Cartesian    normal = surface.getNormal();
			Real area = surface.getArea();
			Point3Cartesian    center = surface.getCenter();

			Vector3Cartesian    value = vectorField(Vector3Cartesian(center.X(), center.Y(), center.Z()));

			Real dotProduct = ScalarProduct(normal, value);

			result = dotProduct * area;

			return result;
		}

		static Real SurfaceIntegralImproveRecursively(const IVectorFunction<3>& vectorField, const RectSurface3D& surface, Real prev_value, int level)
		{
			// now we will calculate integral for surface divided to 4 equal parts
			Real result = 0.0;
			Point3Cartesian pnt_mid12 = (surface._pnt1 + surface._pnt2) / 2.0;
			Point3Cartesian pnt_mid23 = (surface._pnt2 + surface._pnt3) / 2.0;
			Point3Cartesian pnt_mid34 = (surface._pnt3 + surface._pnt4) / 2.0;
			Point3Cartesian pnt_mid41 = (surface._pnt4 + surface._pnt1) / 2.0;
			Point3Cartesian center = surface.getCenter();

			RectSurface3D s1(surface._pnt1, pnt_mid12, center, pnt_mid41);
			RectSurface3D s2(pnt_mid12, surface._pnt2, pnt_mid23, center);
			RectSurface3D s3(center, pnt_mid23, surface._pnt3, pnt_mid34);
			RectSurface3D s4(pnt_mid41, center, pnt_mid34, surface._pnt4);

			result += CalcSurfaceContrib(vectorField, s1);
			result += CalcSurfaceContrib(vectorField, s2);
			result += CalcSurfaceContrib(vectorField, s3);
			result += CalcSurfaceContrib(vectorField, s4);

			// compare to prev_value
			if (fabs(result - prev_value) < 0.0001 || level == 0)
			{
				return result;
			}
			else
			{
				Real new_result = 0.0;

				new_result += SurfaceIntegralImproveRecursively(vectorField, s1, result, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, s2, result, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, s3, result, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, s4, result, level - 1);

				return new_result;
			}

			return result;
		}

		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const RectSurface3D& surface)
		{
			Real result = CalcSurfaceContrib(vectorField, surface);

			return SurfaceIntegralImproveRecursively(vectorField, surface, result, 7);
		}
	};
} // end namespace
#endif