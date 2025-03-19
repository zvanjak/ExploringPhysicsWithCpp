#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/BaseUtils.h"
#include "base/Vector.h"

#include "tools/ConsolePrinter.h"
#include "tools/Visualizer.h"
#include "tools/Serializer.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "algorithms/RootFinding.h"
#include "algorithms/Statistics.h"

#endif

using namespace MML;
using namespace MML::RootFinding;


class DoublePendulumODE : public IODESystem
{
	Real _m1, _m2;
	Real _l1, _l2;
public:
	DoublePendulumODE(Real m1, Real m2, Real l1, Real l2) : _m1(m1), _m2(m2), _l1(l1), _l2(l2) {}

	int  getDim() const override { return 4; }
	void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
	{
		Real th1 = x[0];
		Real w1  = x[1];
		Real th2 = x[2];
		Real w2  = x[3];

		Real g = 9.81;
		Real divisor = (2 * _m1 + _m2 - _m2 * cos(2*th1 - 2*th2));

		dxdt[0] = w1;
		dxdt[1] = (-g * (2 * _m1 + _m2) * sin(th1) - _m2 * g * sin(th1 - 2 * th2) -
							  2 * sin(th1 - th2) * _m2 * (POW2(w2) * _l2 + POW2(w1) * _l1 * cos(th1 - th2)) 
							) / (_l1 * divisor);
		dxdt[2] = w2;
		dxdt[3] = (2 * sin(th1 - th2) * (POW2(w1) * _l1 * (_m1 + _m2) + 
							 g * (_m1 + _m2) * cos(th1) + 
							 POW2(w2) * _l2 * _m2 * cos(th1 - th2)) 
							) / (_l2 * divisor);
	}
};

void Demo_DoublePendulum()
{
	Real	l1 = 1.0, l2 = 1.0;
	Real  m1 = 0.5, m2 = 1.0;
	DoublePendulumODE		odeSysDoublePend = DoublePendulumODE(m1, m2, l1, l2);

	Real	t1 = 0.0, t2 =50.0;
	int   expectNumSteps = 10000;
	Real	minSaveInterval = (t2 - t1) / expectNumSteps;
	Real	initAngle1 = 0.5;
	Real  initAngle2 = 0.101;
	Vector<Real>	initCond{ initAngle1, 0.0, initAngle2, 0.0 };

	ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(odeSysDoublePend);
	ODESystemSolution sol = adaptSolver.integrate(initCond, t1, t2, minSaveInterval, 1e-08, 0.01);

	Vector<Real> t_vals			 = sol.getTValues();
	Vector<Real> theta1_vals = sol.getXValues(0);
	Vector<Real> theta2_vals = sol.getXValues(2);

	std::cout << "\n\n**** Solving double pendulum  ****\n";
	std::vector<ColDesc>				vecNames{ ColDesc("t", 11, 2, 'F'), 
																				ColDesc("theta 1", 15, 8, 'F'), 
																				ColDesc("theta 2", 15, 8, 'F'), };
	std::vector<Vector<Real>*>	vecVals{ &t_vals, &theta1_vals, &theta2_vals };
	
	VerticalVectorPrinter	vvp(vecNames, vecVals);

	//vvp.Print();

	// getting solutions as polynomials
	PolynomInterpRealFunc 	solAdaptPolyInterp0 = sol.getSolutionAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solAdaptPolyInterp1 = sol.getSolutionAsPolyInterp(2, 3);

	Visualizer::VisualizeRealFunction(solAdaptPolyInterp0, 
																		"Double pendulum - theta1 in time",
																		0.0, 10.0, 200, "double_pendulum_angle1.txt");

	// shown together
	Visualizer::VisualizeMultiRealFunction({ &solAdaptPolyInterp0, &solAdaptPolyInterp1 },
																				 "Double pendulum - both angles", 0.0, 10.0, 200,
																				 "double_pendulum_multi_real_func.txt");

	// form parametric curve 2d from solution
	Matrix<Real> curve_points(t_vals.size(), 2);
	for (int i = 0; i < t_vals.size(); i++)
	{
		curve_points(i, 0) = theta1_vals[i];
		curve_points(i, 1) = theta2_vals[i];
	}

	SplineInterpParametricCurve<2> phaseSpaceTrajectory(0.0, 1.0, curve_points);

	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory, "Double pendulum - phase space trajectory",
																		0.0, 1.0, t_vals.size(), "double_pendulum_phase_space.txt");
}