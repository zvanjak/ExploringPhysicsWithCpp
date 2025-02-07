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
		static inline std::string _pathParametricCurveViz{ MML_PATH_ParametricCurveViz };
		static inline std::string _pathVectorFieldViz{ MML_PATH_VectorFieldViz };

	public:
		static void VisualizeRealFunction(const IRealFunction& f, std::string title, Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveRealFuncEquallySpacedDetailed(f, title, x1, x2, numPoints, name);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathRealFuncViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeRealFunction: Not implemented for this OS" << std::endl;
#endif
		}

		static void VisualizeMultiRealFunction(std::vector<IRealFunction*> funcs, std::string title, Real x1, Real x2, int numPoints, std::string fileName)
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

		static void VisualizeScalarFunc2DCartesian(const IScalarFunction<2>& func, std::string title, Real x1, Real x2, int numPointsX, Real y1, Real y2, int numPointsY, std::string fileName)
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

		static void VisualizeVectorField3DCartesian(const IVectorFunction<3>& func, std::string title, Real x1, Real x2, int numPointsX, Real y1, Real y2, int numPointsY, Real z1, Real z2, int numPointsZ, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveVectorFunc3DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, name);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathVectorFieldViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeVectorField3DCartesian: Not implemented for this OS" << std::endl;
#endif
		}

		static void VisualizeParamCurve3D(const IRealToVectorFunction<3>& f, std::string title, Real t1, Real t2, int numPoints, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveParamCurveCartesian3D(f, title, t1, t2, numPoints, name);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathParametricCurveViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeParamCurve3D: Not implemented for this OS" << std::endl;
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
			std::string command = _pathParametricCurveViz + " " + params;
			system(command.c_str());
#else
			std::cout << "VisualizeMultiParamCurve3D: Not implemented for this OS" << std::endl;
#endif
		}

		static void VisualizeODESysSolAsMultiFunc(const ODESystemSolution& sol, std::string title, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			sol.Serialize(name, title);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathRealFuncViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeODESysSolAsMultiFunc: Not implemented for this OS" << std::endl;
#endif
		}

		static void VisualizeODESysSolAsParamCurve3(const ODESystemSolution& sol, std::string title, std::string fileName)
		{
			if (sol._sys_dim != 3)
				throw std::runtime_error("VisualizeODESysSolAsParamCurve3: system dimension must be 3");

			std::string name = _pathResultFiles + fileName;
			sol.SerializeAsParametricCurve3D(name, title);

#if defined(MML_PLATFORM_WINDOWS)
			std::string command = _pathParametricCurveViz + " " + name;
			system(command.c_str());
#else
			std::cout << "VisualizeODESysSolAsParamCurve3: Not implemented for this OS" << std::endl;
#endif
		}
	};
}
#endif