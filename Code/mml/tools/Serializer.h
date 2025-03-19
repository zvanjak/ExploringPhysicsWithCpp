#if !defined  MML_SERIALIZER_H
#define MML_SERIALIZER_H

#include "MMLBase.h"

#include "base/Function.h"

namespace MML
{
	class Serializer
	{
	public:
		// Real function serialization 
		static bool SaveRealMultiFunc(std::vector<IRealFunction*> funcs, std::string title, 
																	Real x1, Real x2, int count, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "MULTI_REAL_FUNCTION_VARIABLE_SPACED" << std::endl;

			file << title << std::endl;
			file << funcs.size() << std::endl;
			file << count << std::endl;
			file << x1 << std::endl;
			file << x2 << std::endl;

			for (int i = 0; i < count; i++)
			{
				double x = x1 + (x2 - x1) * i / (count - 1);
				file << x << " ";
				for (int j = 0; j < funcs.size(); j++)
				{
					file << (*funcs[j])(x) << " ";
				}
				file << std::endl;
			}

			file.close();
			return true;
		}

		static bool SaveRealFuncEquallySpaced(const IRealFunction& f, std::string title, 
																					Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "REAL_FUNCTION_EQUALLY_SPACED" << std::endl;
			file << title << std::endl;
			file << "x1: " << x1 << std::endl;
			file << "x2: " << x2 << std::endl;
			file << "NumPoints: " << numPoints << std::endl;

			Real step = (x2 - x1) / (numPoints - 1);
			for (int i = 0; i < numPoints; i++)
			{
				Real x = x1 + i * step;
				file << f(x) << std::endl;
			}
			file.close();
			return true;
		}

		static bool SaveRealFuncEquallySpacedDetailed(const IRealFunction& f, std::string title, 
																									Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "REAL_FUNCTION_EQUALLY_SPACED_DETAILED" << std::endl;
			file << title << std::endl;
			file << "x1: " << x1 << std::endl;
			file << "x2: " << x2 << std::endl;
			file << "NumPoints: " << numPoints << std::endl;

			Real step = (x2 - x1) / (numPoints - 1);
			for (int i = 0; i < numPoints; i++)
			{
				Real x = x1 + i * step;
				file << x << " " << f(x) << std::endl;
			}
			file.close();
			return true;
		}

		static bool SaveRealFuncVariableSpaced(const IRealFunction& f, std::string title, 
																					 Vector<Real> points, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "REAL_FUNCTION_VARIABLE_SPACED" << std::endl;
			file << title << std::endl;

			for (int i = 0; i < points.size(); i++)
			{
				Real x = points[i];
				file << x << " " << f(x) << std::endl;
			}
			file.close();
			return true;
		}

		// Parametric curve serialization
		template<int N>
		static bool SaveParamCurve(const IRealToVectorFunction<N>& f, std::string inType, std::string title, 
															 Real t1, Real t2, int numPoints, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << inType << std::endl;
			file << title << std::endl;
			file << "t1: " << t1 << std::endl;
			file << "t2: " << t2 << std::endl;
			file << "NumPoints: " << numPoints << std::endl;

			Real delta = (t2 - t1) / (numPoints - 1);
			for (Real t = t1; t <= t2; t += delta)
			{
				file << t << " ";
				for (int i = 0; i < N; i++)
					file << f(t)[i] << " ";
				file << std::endl;
			}
			file.close();
			return true;
		}

		template<int N>
		static bool SaveParamCurve(const IRealToVectorFunction<N>& f, std::string inType, std::string title, 
															 Vector<Real> points, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << inType << std::endl;
			file << title << std::endl;
			file << "t1: " << points[0] << std::endl;
			file << "t2: " << points[points.size() - 1] << std::endl;
			file << "NumPoints: " << points.size() << std::endl;
			for (int i = 0; i < points.size(); i++)
			{
				Real t = points[i];
				file << t << " ";
				for (int i = 0; i < N; i++)
					file << f(t)[i] << " ";
				file << std::endl;
			}

			file.close();
			return true;
		}

		template<int N>
		static bool SaveAsParamCurve(std::vector<VectorN<Real, N>> vals, std::string inType, std::string title, 
																 Real t1, Real t2, int numPoints, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << inType << std::endl;
			file << title << std::endl;
			file << "t1: " << t1 << std::endl;
			file << "t2: " << t2 << std::endl;
			file << "NumPoints: " << numPoints << std::endl;

			Real delta = (t2 - t1) / (numPoints - 1);
			for (int i = 0; i < numPoints; i++)
			{
				Real t = t1 + i * delta;
				file << t << " ";
				for (int j = 0; j < N; j++)
					file << vals[i][j] << " ";
				file << std::endl;
			}
			file.close();
			return true;
		}

		static bool SaveAsParamCurve2D(Vector<Real> vec_x, Vector<Real> vec_y)
		{

		}

		//bool SerializeCartesian2D(Real t1, Real t2, int numPoints, std::string fileName) const
		//{
		//    return SaveParamCurve<2>("PARAMETRIC_CURVE_CARTESIAN_2D", t1, t2, numPoints, fileName);
		//}
		//bool SerializeCartesian2DAtPoints(Vector<Real> points, std::string fileName) const
		//{
		//    return SaveParamCurveAtPoints<2>("PARAMETRIC_CURVE_CARTESIAN_2D_AT_POINTS", points, fileName);
		//}
		static bool SaveParamCurveCartesian2D(const IRealToVectorFunction<2>& f, std::string title, 
																					Real t1, Real t2, int numPoints, std::string fileName)
		{
			return SaveParamCurve<2>(f, "PARAMETRIC_CURVE_CARTESIAN_2D", title, t1, t2, numPoints, fileName);
		}
		static bool SaveParamCurveCartesian3D(const IRealToVectorFunction<3>& f, std::string title, 
																					Real t1, Real t2, int numPoints, std::string fileName)
		{
			return SaveParamCurve<3>(f, "PARAMETRIC_CURVE_CARTESIAN_3D", title, t1, t2, numPoints, fileName);
		}
		//bool SerializeCartesian3DAtPoints(Vector<Real> points, std::string fileName) const
		//{
		//    return SaveParamCurveAtPoints("PARAMETRIC_CURVE_CARTESIAN_3D", points, fileName);
		//}
		//bool SerializePolar(Real t1, Real t2, int numPoints, std::string fileName) const
		//{
		//    return SaveParamCurve("PARAMETRIC_CURVE_POLAR", t1, t2, numPoints, fileName);
		//}
		//bool SerializePolarAtPoints(Vector<Real> points, std::string fileName) const
		//{
		//    return SaveParamCurveAtPoints("PARAMETRIC_CURVE_POLAR", points, fileName);
		//}
		//bool SerializeSpherical(Real t1, Real t2, int numPoints, std::string fileName) const
		//{
		//    return SaveParamCurve("PARAMETRIC_CURVE_SPHERICAL", t1, t2, numPoints, fileName);
		//}
		//bool SerializeSphericalAtPoints(Vector<Real> points, std::string fileName) const
		//{
		//    return SaveParamCurveAtPoints("PARAMETRIC_CURVE_SPHERICAL", points, fileName);
		//}

		// Scalar function serialization
		static bool SaveScalarFunc2DCartesian(const IScalarFunction<2>& f, std::string title, 
																					Real x1, Real x2, int numPointsX, 
																					Real y1, Real y2, int numPointsY, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "SCALAR_FUNCTION_CARTESIAN_2D" << std::endl;
			file << title << std::endl;
			file << "x1: " << x1 << std::endl;
			file << "x2: " << x2 << std::endl;
			file << "NumPointsX: " << numPointsX << std::endl;
			file << "y1: " << y1 << std::endl;
			file << "_secDerY: " << y2 << std::endl;
			file << "NumPointsY: " << numPointsY << std::endl;

			Real stepX = (x2 - x1) / (numPointsX - 1);
			Real stepY = (y2 - y1) / (numPointsY - 1);
			for (int i = 0; i < numPointsX; i++)
			{
				for (int j = 0; j < numPointsY; j++)
				{
					Real x = x1 + i * stepX;
					Real y = y1 + j * stepY;
					file << x << " " << y << " " << f(VectorN<Real, 2>{x, y}) << std::endl;
				}
			}
			return true;
		}

		static bool SaveScalarFunc3DCartesian(const IScalarFunction<3>& f, std::string title, 
																					Real x1, Real x2, int numPointsX, 
																					Real y1, Real y2, int numPointsY, 
																					Real z1, Real z2, int numPointsZ, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "SCALAR_FUNCTION_CARTESIAN_3D" << std::endl;
			file << title << std::endl;
			file << "x1: " << x1 << std::endl;
			file << "x2: " << x2 << std::endl;
			file << "NumPointsX: " << numPointsX << std::endl;
			file << "y1: " << y1 << std::endl;
			file << "_secDerY: " << y2 << std::endl;
			file << "NumPointsY: " << numPointsY << std::endl;
			file << "z1: " << z1 << std::endl;
			file << "z2: " << z2 << std::endl;
			file << "NumPointsZ: " << numPointsZ << std::endl;

			Real stepX = (x2 - x1) / (numPointsX - 1);
			Real stepY = (y2 - y1) / (numPointsY - 1);
			Real stepZ = (z2 - z1) / (numPointsZ - 1);
			for (int i = 0; i < numPointsX; i++)
			{
				for (int j = 0; j < numPointsY; j++)
				{
					for (int k = 0; k < numPointsZ; k++)
					{
						Real x = x1 + i * stepX;
						Real y = y1 + j * stepY;
						Real z = z1 + j * stepZ;
						file << x << " " << y << " " << z << " " << f(VectorN<Real, 3>{x, y, z}) << std::endl;
					}
				}
			}
			return true;
		}

		// vector function serialization
		static bool SaveVectorFunc3D(const IVectorFunction<3>& f, std::string inType, std::string title, 
																 Real x1_start, Real x1_end, int numPointsX1, 
																 Real x2_start, Real x2_end, int numPointsX2, 
																 Real x3_start, Real x3_end, int numPointsX3, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << inType << std::endl;
			file << title << std::endl;

			Real stepX = (x1_end - x1_start) / (numPointsX1 - 1);
			Real stepY = (x2_end - x2_start) / (numPointsX2 - 1);
			Real stepZ = (x3_end - x3_start) / (numPointsX3 - 1);
			for (int i = 0; i < numPointsX1; i++)
				for (int j = 0; j < numPointsX2; j++)
					for (int k = 0; k < numPointsX3; k++)
					{
						Real x = x1_start + i * stepX;
						Real y = x2_start + j * stepY;
						Real z = x3_start + k * stepZ;
						auto val = f(VectorN<Real, 3>{x, y, z});
						file << x << " " << y << " " << z << " " << val[0] << " " << val[1] << " " << val[2] << std::endl;
					}

			file.close();
			return true;
		}

		static bool SaveVectorFunc3D(const IVectorFunction<3>& f, std::string inType, std::string title, 
																 Real x1_start, Real x1_end, int numPointsX1, 
																 Real x2_start, Real x2_end, int numPointsX2, 
																 Real x3_start, Real x3_end, int numPointsX3, 
																 std::string fileName, Real upper_threshold)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << inType << std::endl;
			file << title << std::endl;

			Real stepX = (x1_end - x1_start) / (numPointsX1 - 1);
			Real stepY = (x2_end - x2_start) / (numPointsX2 - 1);
			Real stepZ = (x3_end - x3_start) / (numPointsX3 - 1);
			for (int i = 0; i < numPointsX1; i++)
				for (int j = 0; j < numPointsX2; j++)
					for (int k = 0; k < numPointsX3; k++)
					{
						Real x = x1_start + i * stepX;
						Real y = x2_start + j * stepY;
						Real z = x3_start + k * stepZ;
						auto val = f(VectorN<Real, 3>{x, y, z});

						if (val.NormL2() < upper_threshold)
							file << x << " " << y << " " << z << " " << val[0] << " " << val[1] << " " << val[2] << std::endl;
					}

			file.close();
			return true;
		}

		static bool SaveVectorFunc3DCartesian(const IVectorFunction<3>& f, std::string title, 
																					Real x1, Real x2, int numPointsX, 
																					Real y1, Real y2, int numPointsY, 
																					Real z1, Real z2, int numPointsZ, std::string fileName)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_3D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, fileName);
		}
		static bool SaveVectorFunc3DCartesian(const IVectorFunction<3>& f, std::string title, 
																					Real x1, Real x2, int numPointsX, 
																					Real y1, Real y2, int numPointsY, 
																					Real z1, Real z2, int numPointsZ, std::string fileName, Real upper_threshold)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_3D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, fileName, upper_threshold);
		}
		
		
		static bool SaveVectorFuncSpherical(const IVectorFunction<3>& f, std::string title, Real r1, Real r2, int numPointsR, Real theta1, Real theta2, int numPointsTheta, Real phi1, Real phi2, int numPointsPhi, std::string fileName)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_SPHERICAL", title, r1, r2, numPointsR, theta1, theta2, numPointsTheta, phi1, phi2, numPointsPhi, fileName);
		}
		static bool SaveVectorFuncSpherical(const IVectorFunction<3>& f, std::string title, Real r1, Real r2, int numPointsR, Real theta1, Real theta2, int numPointsTheta, Real phi1, Real phi2, int numPointsPhi, std::string fileName, Real upper_threshold)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_SPHERICAL", title, r1, r2, numPointsR, theta1, theta2, numPointsTheta, phi1, phi2, numPointsPhi, fileName, upper_threshold);
		}
	};
}
#endif 