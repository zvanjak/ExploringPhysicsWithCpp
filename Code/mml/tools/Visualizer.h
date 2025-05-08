#if !defined MML_VIZUALIZER_H
#define MML_VIZUALIZER_H

#include "MMLBase.h"

#include "base/Function.h"

#include "core/ODESystem.h"
#include "core/ODESystemSolution.h"

#include "tools/Serializer.h"


namespace MML
{
	class Visualizer
	{
		static inline std::string _pathResultFiles{ MML_PATH_ResultFiles };

		static inline std::string _pathRealFuncViz{ MML_PATH_RealFuncViz };
		static inline std::string _pathSurfaceViz{ MML_PATH_SurfaceViz };
		static inline std::string _pathParametricCurve3DViz{ MML_PATH_ParametricCurve3DViz };
		static inline std::string _pathParametricCurve2DViz{ MML_PATH_ParametricCurve2DViz };
		static inline std::string _pathVectorField2DViz{ MML_PATH_VectorField2DViz };
		static inline std::string _pathVectorField3DViz{ MML_PATH_VectorField3DViz };

		static inline std::string _pathParticle2DViz{ MML_PATH_Particle2DViz };
		static inline std::string _pathParticle3DViz{ MML_PATH_Particle3DViz };

	public:
		// visualizations of Real function
		static void VisualizeRealFunction(const IRealFunction& f, std::string title,
																			Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveRealFunc(f, title, x1, x2, numPoints, name);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathRealFuncViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeRealFunction: Not implemented for this OS" << std::endl;
#endif
		}

		static void VisualizeRealFunction(const IRealFunction& f, std::string title,
																			Vector<Real> points, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveRealFunc(f, title, points, name);
#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathRealFuncViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeRealFunction: Not implemented for this OS" << std::endl;
#endif
		}

		static void VisualizeMultiRealFunction(std::vector<IRealFunction*> funcs, std::string title,
																					 std::vector<std::string> func_legend,
																					 Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveRealMultiFunc(funcs, title, x1, x2, numPoints, name);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathRealFuncViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeMultiRealFunction: Not implemented for this OS" << std::endl;
#endif
		}

		// visualizations of Scalar function in 2D
		static void VisualizeScalarFunc2DCartesian(const IScalarFunction<2>& func, std::string title,
																							 Real x1, Real x2, int numPointsX,
																							 Real y1, Real y2, int numPointsY, std::string fileName)
																						{
																							std::string name = _pathResultFiles + fileName;
																							Serializer::SaveScalarFunc2DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, name);

																				#if defined(MML_PLATFORM_WINDOWS)
																							std::string command = _pathSurfaceViz + " " + name;
																							system(command.c_str());
																				#else
																							std::cout << "VisualizeScalarFunc2DCartesian: Not implemented for this OS" << std::endl;
																				#endif
																						}

		// visualizations of Vector fields
		static void VisualizeVectorField2DCartesian(const IVectorFunction<2>& func, std::string title,
																								Real x1, Real x2, int numPointsX,
																								Real y1, Real y2, int numPointsY, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveVectorFunc2DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, name);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathVectorField2DViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeVectorField2DCartesian: Not implemented for this OS" << std::endl;
#endif
		}

		static void VisualizeVectorField3DCartesian(const IVectorFunction<3>& func, std::string title,
																								Real x1, Real x2, int numPointsX,
																								Real y1, Real y2, int numPointsY,
																								Real z1, Real z2, int numPointsZ, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveVectorFunc3DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, name, 3.0);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathVectorField3DViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeVectorField3DCartesian: Not implemented for this OS" << std::endl;
#endif
		}

		// visualizations of Parametric curves
		static void VisualizeParamCurve2D(const IRealToVectorFunction<2>& f, std::string title, 
																			Real t1, Real t2, int numPoints, 
																			std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveParamCurveCartesian2D(f, title, t1, t2, numPoints, name);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathParametricCurve2DViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeParamCurve3D: Not implemented for this OS" << std::endl;
#endif
		}

		static void VisualizeMultiParamCurve2D(std::vector<IRealToVectorFunction<2>*> curves, 
																					 std::string title,
																					 Real t1, Real t2,
																					 int numPoints, std::string fileName)
		{
			// for each function, serialize data to file, with name generated from fileName
			int i = 0;
			std::string params = " ";
			for (auto& curve : curves)
			{
				std::string name = _pathResultFiles + fileName + "_" + std::to_string(i + 1) + ".txt";
				Serializer::SaveParamCurveCartesian2D(*curve, title, t1, t2, numPoints, name);
				i++;

				params = params + name + " ";
			}

			// then call the visualizer with all the file names
#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathParametricCurve2DViz + params;
			system(command.c_str());
#else
			std::cout << "VisualizeMultiParamCurve2D: Not implemented for this OS" << std::endl;
#endif
		}

		static void VisualizeParamCurve3D(const IRealToVectorFunction<3>& f, std::string title, 
																			Real t1, Real t2, int numPoints, 
																			std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveParamCurveCartesian3D(f, title, t1, t2, numPoints, name);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathParametricCurve3DViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeParamCurve3D: Not implemented for this OS" << std::endl;
#endif
		}

		static void VisualizeMultiParamCurve3D(std::vector<IRealToVectorFunction<3>*> curves, 
																					 std::string title,
																					 Real t1, Real t2, int numPoints, 
																					 std::string fileName)
		{
			// for each function, serialize data to file, with name generated from fileName
			int i = 0;
			std::string params = " ";
			for (auto& curve : curves)
			{
				std::string name = _pathResultFiles + fileName + "_" + std::to_string(i + 1);
				Serializer::SaveParamCurveCartesian3D(*curve, title, t1, t2, numPoints, name);
				i++;
				params = params + name + " ";
			}
			// then call the visualizer with all the file names
#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathParametricCurve3DViz + params;
			system(command.c_str());
#else
			std::cout << "VisualizeMultiParamCurve3D: Not implemented for this OS" << std::endl;
#endif
		}

		static void VisualizeMultiParamCurve3D(std::vector<std::string> fileNames)
		{
			std::string params;
			for (auto& name : fileNames)
			{
				params = params + (_pathResultFiles + name + " ");
			}

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathParametricCurve3DViz + " " + params;
			system(command.c_str());
#else
			std::cout << "VisualizeMultiParamCurve3D: Not implemented for this OS" << std::endl;
#endif
		}

		// ODE Solution visualizations
		// Visualizing single variable of ODE system solution as Real function
		static void VisualizeODESysSolCompAsFunc(const ODESystemSolution& sol, int compInd,
																						 std::string title, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveODESolutionComponentAsFunc(sol, compInd, title, name);
#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathRealFuncViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeODESysSolCompAsFunc: Not implemented for this OS" << std::endl;
#endif
		}

		// Visualizing ODE system solution as a multi-function (all variables)
		static void VisualizeODESysSolAsMultiFunc(const ODESystemSolution& sol,
																							std::string title, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveODESolutionAsMultiFunc(sol, title, name);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathRealFuncViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeODESysSolAsMultiFunc: Not implemented for this OS" << std::endl;
#endif
		}

		// Visualizing two variables of ODE system solution as a parametric curve in 2D
		static void VisualizeODESysSolAsParamCurve2(const ODESystemSolution& sol,
																								int ind1, int ind2,
																								std::string title, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveODESolAsParametricCurve2D(sol, name, ind1, ind2, title);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathParametricCurve2DViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeODESysSolAsParamCurve2: Not implemented for this OS" << std::endl;
#endif
		}

		// Visualizing three variables of ODE system solution as a parametric curve in 3D
		static void VisualizeODESysSolAsParamCurve3(const ODESystemSolution& sol,
																								int ind1, int ind2, int ind3,
																								std::string title, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveODESolAsParametricCurve3D(sol, name, ind1, ind2, ind3, title);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathParametricCurve3DViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeODESysSolAsParamCurve3: Not implemented for this OS" << std::endl;
#endif
		}

		// Particle simulation visualizations
		static void VisualizeParticleSimulation2D(std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathParticle2DViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeParticleSimulation2D: Not implemented for this OS" << std::endl;
#endif
		}

		static void VisualizeParticleSimulation3D(std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathParticle3DViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeParticleSimulation3D: Not implemented for this OS" << std::endl;
#endif
		}
	};
}
#endif