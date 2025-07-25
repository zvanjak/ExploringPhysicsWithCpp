#if !defined  MML_ODE_SYSTEM_SOLUTON_H
#define MML_ODE_SYSTEM_SOLUTON_H

#include "MMLBase.h"

#include "base/InterpolatedFunction.h"
#include "base/ODESystem.h"


namespace MML
{
	class ODESystemSolution
	{
		int  _numStepsOK, _numStepsBad;

		int  _sys_dim;
		int  _totalSavedSteps;
	
		Real _t1, _t2;
		Vector<Real> _tval;
		Matrix<Real> _xval;
	
	public:
		ODESystemSolution(Real x1, Real x2, int dim)
			: _t1(x1), _t2(x2), _sys_dim(dim),
				_numStepsOK(0), _numStepsBad(0), _totalSavedSteps(0) {
		}
		ODESystemSolution(Real x1, Real x2, int dim, int maxSteps)
			: _t1(x1), _t2(x2), _sys_dim(dim),
				_numStepsOK(0), _numStepsBad(0), _totalSavedSteps(maxSteps + 1)
		{
			_tval.Resize(maxSteps + 1);
			_xval.Resize(dim, maxSteps + 1);
		}

		void incrementNumStepsOK()	{ _numStepsOK++; }
		void incrementNumStepsBad() { _numStepsBad++; }

		int getSysDym() const { return _sys_dim; }
		
		Real getT1() const { return _t1; }
		Real getT2() const { return _t2; }

		int getNumStepsOK() const { return _numStepsOK; }
		int getNumStepsBad() const { return _numStepsBad; }
		int getTotalNumSteps() const { return _numStepsOK + _numStepsBad; }
		int getTotalSavedSteps() const { return _totalSavedSteps; }

		Vector<Real> getTValues() const { return _tval; }
		Matrix<Real> getXValues() const { return _xval; }
		Vector<Real> getXValues(int component) const 
		{ 
			if (component >= _sys_dim)
				throw("Index out of range in ODESystemSolution::getXValues(component)");
			
			return _xval.VectorFromRow(component); 
		}

		void fillValues(int ind, Real x, Vector<Real>& y)
		{
			if (ind < 0 || ind > _totalSavedSteps)
				throw("Index out of range in ODESystemSolution::fillValues");

			_tval[ind] = x;
			for (int i = 0; i < _sys_dim; i++)
				_xval[i][ind] = y[i];
		}
		void setFinalSize(int numDoneSteps)
		{
			_totalSavedSteps = numDoneSteps + 1;

			Vector<Real> new_xval(numDoneSteps + 1);
			Matrix<Real> new_yval(_sys_dim, numDoneSteps + 1);

			for (int i = 0; i < numDoneSteps + 1; i++)
			{
				new_xval[i] = _tval[i];
				for (int j = 0; j < _sys_dim; j++)
					new_yval[j][i] = _xval[j][i];
			}

			_tval = new_xval;
			_xval = new_yval;
		}

		LinearInterpRealFunc  getSolutionAsLinearInterp(int component) const
		{
			Vector<Real> xsave = _tval;
			Vector<Real> ysave = _xval.VectorFromRow(component);

			return LinearInterpRealFunc(xsave, ysave);
		}
		PolynomInterpRealFunc getSolutionAsPolyInterp(int component, int polyOrder) const
		{
			Vector<Real> xsave = _tval;
			Vector<Real> ysave = _xval.VectorFromRow(component);
			return PolynomInterpRealFunc(xsave, ysave, polyOrder + 1);
		}
		SplineInterpRealFunc  getSolutionAsSplineInterp(int component) const
		{
			Vector<Real> xsave = _tval;
			Vector<Real> ysave = _xval.VectorFromRow(component);
			return SplineInterpRealFunc(xsave, ysave);
		}

		SplineInterpParametricCurve<2> getSolutionAsParametricCurve2D(int ind1, int ind2) const
		{
			Matrix<Real> curve_points(_totalSavedSteps, 2);
			for (int i = 0; i < _totalSavedSteps; i++)
			{
				curve_points(i, 0) = _xval[ind1][i];
				curve_points(i, 1) = _xval[ind2][i];
			}
			return SplineInterpParametricCurve<2>(_t1, _t2, curve_points);
		}

		SplineInterpParametricCurve<3> getSolutionAsParametricCurve3D(int ind1, int ind2, int ind3) const
		{
			Matrix<Real> curve_points(_totalSavedSteps, 3);
			for (int i = 0; i < _totalSavedSteps; i++)
			{
				curve_points(i, 0) = _xval[ind1][i];
				curve_points(i, 1) = _xval[ind2][i];
				curve_points(i, 2) = _xval[ind3][i];
			}
			return SplineInterpParametricCurve<3>(_t1, _t2, curve_points);
		}

		//bool Serialize(std::string fileName, std::string title) const
		//{
		//		std::ofstream file(fileName);
		//		if (!file.is_open())
		//			return false;

		//		file << "MULTI_REAL_FUNCTION" << std::endl;

		//		file << title << std::endl;
		//		file << _sys_dim << std::endl;
		//		file << "x1: " << _t1 << std::endl;
		//		file << "x2: " << _t2 << std::endl;
		//		file << "NumPoints: " << _totalSavedSteps << std::endl;

		//		for (int i = 0; i < _totalSavedSteps; i++)
		//		{
		//			file << _tval[i] << " ";
		//			for (int j = 0; j < _sys_dim; j++)
		//			{
		//				file << _xval[j][i] << " ";
		//			}
		//			file << std::endl;
		//		}
		//		file.close();
		//		return true;
		//}

		//bool SerializeAsParametricCurve2D(std::string fileName, int ind1, int ind2, std::string title) const
		//{
		//	std::ofstream file(fileName);
		//	if (!file.is_open())
		//		return false;
		//	
		//	file << "PARAMETRIC_CURVE_CARTESIAN_2D" << std::endl;
		//	
		//	file << title << std::endl;
		//	file << "t1: " << _t1 << std::endl;
		//	file << "t2: " << _t2 << std::endl;
		//	file << "NumPoints: " << _totalSavedSteps << std::endl;
		//	for (int i = 0; i < _totalSavedSteps; i++)
		//	{
		//		file << _tval[i] << " ";
		//		file << _xval[ind1][i] << " ";
		//		file << _xval[ind2][i] << " ";
		//		file << std::endl;
		//	}
		//	file.close();
		//	return true;
		//}

		//bool SerializeAsParametricCurve3D(std::string fileName, int ind1, int ind2, int ind3, std::string title) const
		//{
		//	std::ofstream file(fileName);
		//	if (!file.is_open())
		//		return false;
		//	
		//	file << "PARAMETRIC_CURVE_CARTESIAN_3D" << std::endl;
		//	
		//	file << title << std::endl;
		//	file << "t1: " << _t1 << std::endl;
		//	file << "t2: " << _t2 << std::endl;
		//	file << "NumPoints: " << _totalSavedSteps << std::endl;
		//	for (int i = 0; i < _totalSavedSteps; i++)
		//	{
		//		file << _tval[i] << " ";

		//		file << _xval[ind1][i] << " ";
		//		file << _xval[ind2][i] << " ";
		//		file << _xval[ind3][i] << " ";
		//		
		//		file << std::endl;
		//	}
		//	file.close();
		//	return true;
		//}

		//bool SerializeAsParametricCurve3D(std::string fileName, std::string title) const
		//{
		//	if (_sys_dim != 3)
		//		return false;

		//	std::ofstream file(fileName);
		//	if (!file.is_open())
		//		return false;

		//	file << "PARAMETRIC_CURVE_CARTESIAN_3D" << std::endl;
		//	file << "t1: " << _t1 << std::endl;
		//	file << "t2: " << _t2 << std::endl;
		//	file << "NumPoints: " << _totalSavedSteps << std::endl;

		//	for (int i = 0; i < _totalSavedSteps; i++)
		//	{
		//		file << _tval[i] << " ";
		//		for (int j = 0; j < _sys_dim; j++)
		//		{
		//			file << _xval[j][i] << " ";
		//		}
		//		file << std::endl;
		//	}

		//	file.close();
		//	return true;
		//}
	};
}

#endif