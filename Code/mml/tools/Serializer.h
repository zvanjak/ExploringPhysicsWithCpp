#if !defined  MML_SERIALIZER_H
#define MML_SERIALIZER_H

#include "MMLBase.h"

#include "base/Function.h"

#include "core/ODESystemSolution.h"

namespace MML
{
	class Serializer
	{
	public:
		// Real function serialization 
		// serializing values at equally spaced points in given interval
		static bool SaveRealFunc(const IRealFunction& f, std::string title,
														 Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "REAL_FUNCTION" << std::endl;
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

		// serializing values at given list of points
		static bool SaveRealFunc(const IRealFunction& f, std::string title,
														 Vector<Real> points, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "REAL_FUNCTION" << std::endl;
			file << title << std::endl;
			file << "x1: " << points[0] << std::endl;
			file << "x2: " << points[points.size() - 1] << std::endl;
			file << "NumPoints: " << points.size() << std::endl;

			for (int i = 0; i < points.size(); i++)
			{
				Real x = points[i];
				file << x << " " << f(x) << std::endl;
			}
			file.close();
			return true;
		}

		// serializing multiple functions in a single files
		static bool SaveRealMultiFunc(std::vector<IRealFunction*> funcs, std::string title,
																	Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "MULTI_REAL_FUNCTION" << std::endl;

			file << title << std::endl;
			file << funcs.size() << std::endl;
			file << "x1: " << x1 << std::endl;
			file << "x2: " << x2 << std::endl;
			file << "NumPoints: " << numPoints << std::endl;

			for (int i = 0; i < numPoints; i++)
			{
				double x = x1 + (x2 - x1) * i / (numPoints - 1);
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

		// same as SaveRealFunc, but points are not explicitly written in file
		// (as they can be calculated fomr x1, x2 and numPoints)
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

		template<int N>
		static bool SaveAsParamCurve(std::vector<VectorN<Real, N>> vals, std::string inType, std::string title,
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
				for (int j = 0; j < N; j++)
					file << vals[i][j] << " ";
				file << std::endl;
			}
			file.close();
			return true;
		}

		// save parametric curve in 2D as a list of points (in interval [0, 1.0])
		static bool SaveAsParamCurve2D(Vector<Real> vec_x, Vector<Real> vec_y)
		{
			std::ofstream file("PARAMETRIC_CURVE_CARTESIAN_2D");
			if (!file.is_open())
				return false;
			
			file << "PARAMETRIC_CURVE_CARTESIAN_2D" << std::endl;
			file << "t1: " << 0.0 << std::endl;
			file << "t2: " << 1.0 << std::endl;
			file << "NumPoints: " << vec_x.size() << std::endl;
			for (int i = 0; i < vec_x.size(); i++)
			{
				file << i / (vec_x.size() - 1) << " " << vec_x[i] << " " << vec_y[i] << std::endl;
			}
			file.close();
			return true;
		}

		// Helper/forwarding functions
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
			file << "y2: " << y2 << std::endl;
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

		// 2D vector function serialization
		static bool SaveVectorFunc2D( const IVectorFunction<2>& f, std::string inType, std::string title,
																	Real x1_start, Real x1_end, int numPointsX1,
																	Real x2_start, Real x2_end, int numPointsX2, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;
			file << inType << std::endl;
			file << title << std::endl;
			Real stepX = (x1_end - x1_start) / (numPointsX1 - 1);
			Real stepY = (x2_end - x2_start) / (numPointsX2 - 1);
			for (int i = 0; i < numPointsX1; i++)
				for (int j = 0; j < numPointsX2; j++)
				{
					Real x = x1_start + i * stepX;
					Real y = x2_start + j * stepY;
					auto val = f(VectorN<Real, 2>{x, y});
					file << x << " " << y << " " << val[0] << " " << val[1] << std::endl;
				}
			file.close();
			return true;
		}

		static bool SaveVectorFunc2D( const IVectorFunction<2>& f, std::string inType, std::string title,
																	Real x1_start, Real x1_end, int numPointsX1,
																	Real x2_start, Real x2_end, int numPointsX2,
																	std::string fileName, Real upper_threshold)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;
			file << inType << std::endl;
			file << title << std::endl;
			Real stepX = (x1_end - x1_start) / (numPointsX1 - 1);
			Real stepY = (x2_end - x2_start) / (numPointsX2 - 1);
			for (int i = 0; i < numPointsX1; i++)
				for (int j = 0; j < numPointsX2; j++)
				{
					Real x = x1_start + i * stepX;
					Real y = x2_start + j * stepY;
					auto val = f(VectorN<Real, 2>{x, y});
					if (val.NormL2() < upper_threshold)
						file << x << " " << y << " " << val[0] << " " << val[1] << std::endl;
				}
			file.close();
			return true;
		}

		static bool SaveVectorFunc2DCartesian(const IVectorFunction<2>& f, std::string title,
																					Real x1, Real x2, int numPointsX,
																					Real y1, Real y2, int numPointsY, std::string fileName)
		{
			return SaveVectorFunc2D(f, "VECTOR_FIELD_2D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, fileName);
		}

		static bool SaveVectorFunc2DCartesian(const IVectorFunction<2>& f, std::string title,
																					Real x1, Real x2, int numPointsX,
																					Real y1, Real y2, int numPointsY, std::string fileName, 
																					Real upper_threshold)
		{
			return SaveVectorFunc2D(f, "VECTOR_FIELD_2D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, fileName, upper_threshold);
		}

		// 3D vector function serialization
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
																					Real z1, Real z2, int numPointsZ, std::string fileName, 
																					Real upper_threshold)
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

		// ODESolution serialization
		static bool SaveODESolutionComponentAsFunc(const ODESystemSolution& sol, int compInd, 
																							 std::string title, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;
			
			file << "REAL_FUNCTION" << std::endl;
			file << title << std::endl;
			file << "x1: " << sol.getT1() << std::endl;
			file << "x2: " << sol.getT2() << std::endl;
			file << "NumPoints: " << sol.getTotalSavedSteps() << std::endl;
			
			for (int i = 0; i < sol.getTotalSavedSteps(); i++)
			{
				file << sol.getTValues()[i] << " " << sol.getXValues()[compInd][i] << std::endl;
			}
			
			file.close();
			return true;
		}

		static bool SaveODESolutionAsMultiFunc(const ODESystemSolution& sol, std::string title, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "MULTI_REAL_FUNCTION" << std::endl;

			file << title << std::endl;
			file << sol.getSysDym() << std::endl;
			file << "x1: " << sol.getT1() << std::endl;
			file << "x2: " << sol.getT2() << std::endl;
			file << "NumPoints: " << sol.getTotalSavedSteps() << std::endl;

			for (int i = 0; i < sol.getTotalSavedSteps(); i++)
			{
				file << sol.getTValues()[i] << " ";
				for (int j = 0; j < sol.getSysDym(); j++)
				{
					file << sol.getXValues()[j][i] << " ";
				}
				file << std::endl;
			}
			file.close();
			return true;
		}

		static bool SaveODESolAsParametricCurve2D(const ODESystemSolution& sol, std::string fileName, 
																							int ind1, int ind2, std::string title)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "PARAMETRIC_CURVE_CARTESIAN_2D" << std::endl;
			
			file << title << std::endl;
			file << "t1: " << sol.getT1() << std::endl;
			file << "t2: " << sol.getT2() << std::endl;
			file << "NumPoints: " << sol.getTotalSavedSteps() << std::endl;
			
			for (int i = 0; i < sol.getTotalSavedSteps(); i++)
			{
				file << sol.getTValues()[i] << " ";
				file << sol.getXValues()[ind1][i] << " ";
				file << sol.getXValues()[ind2][i] << " ";
				file << std::endl;
			}
			
			file.close();
			return true;
		}
		
		static bool SaveODESolAsParametricCurve3D(const ODESystemSolution& sol, std::string fileName, 
																							int ind1, int ind2, int ind3, std::string title)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "PARAMETRIC_CURVE_CARTESIAN_3D" << std::endl;

			file << title << std::endl;
			file << "t1: " << sol.getT1() << std::endl;
			file << "t2: " << sol.getT2() << std::endl;
			file << "NumPoints: " << sol.getTotalSavedSteps() << std::endl;

			for (int i = 0; i < sol.getTotalSavedSteps(); i++)
			{
				file << sol.getTValues()[i] << " ";
				file << sol.getXValues()[ind1][i] << " ";
				file << sol.getXValues()[ind2][i] << " ";
				file << sol.getXValues()[ind3][i] << " ";
				file << std::endl;
			}

			file.close();
			return true;
		}

		// particle simulation serialization
		static bool SaveParticleSimulation2D(std::string fileName, int numBalls, 
																				 std::vector<std::vector<Pnt2Cart>> ballPositions, 
																				 std::vector<std::string> ballColors, std::vector<double> ballRadius )
		{
			std::ofstream file(fileName);
			if (file.is_open())
			{
				file << "PARTICLE_SIMULATION_DATA_2D" << std::endl;
				file << "NumBalls: " << numBalls << std::endl;

				for (int i=0; i<numBalls; i++)
				{
					file << "Ball_" << " " << ballColors[i] << " " << ballRadius[i] << std::endl;
				}

				int numSteps = ballPositions[0].size();
				file << "NumSteps: " << numSteps << std::endl;

				for (int i = 0; i < numSteps; i++)
				{
					file << "Step " << i << " 0.1" << std::endl;
					for (int j = 0; j < numBalls; j++)
					{
						file << j << " " << ballPositions[j][i].X() << " " << ballPositions[j][i].Y() << "\n";
					}
				}
				file.close();
			}
			else
			{
				std::cerr << "Unable to open file";
				return false;
			}
			return true;
		}

		static bool SaveParticleSimulation3D(std::string fileName, int numBalls, 
																				 std::vector<std::vector<Pnt3Cart>> ballPositions,
																				 std::vector<std::string> ballColors, std::vector<double> ballRadius)
		{
			std::ofstream file(fileName);
			if (file.is_open())
			{
				file << "PARTICLE_SIMULATION_DATA_3D" << std::endl;
				file << "NumBalls: " << numBalls << std::endl;

				for (int i = 0; i < numBalls; i++)
				{
					file << "Ball_" << i+1 << " " << ballColors[i] << " " << ballRadius[i] << std::endl;
				}

				int numSteps = ballPositions[0].size();
				file << "NumSteps: " << numSteps << std::endl;

				for (int i = 0; i < numSteps; i++)
				{
					file << "Step " << i << " 0.1" << std::endl;
					for (int j = 0; j < numBalls; j++)
					{
						file << j << " " << ballPositions[j][i].X() << " " << ballPositions[j][i].Y() << " " << ballPositions[j][i].Z() << "\n";
					}
				}
				file.close();
			}
			else
			{
				std::cerr << "Unable to open file";
				return false;
			}
			return true;
		}
	};
}
#endif 