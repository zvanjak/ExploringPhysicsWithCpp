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
#include "algorithms/FunctionsAnalyzer.h"

#include "Mechanics/Pendulums.h"
#endif

using namespace MML;
using namespace MML::RootFinding;

using namespace MPL;

void Demo_SimplePendulum_Euler()
{
	// alternative way of defining the system with lambda function
	ODESystem pendSys(2, [](Real t, const Vector<Real>& x, Vector<Real>& dxdt)
		{
			Real Length = 1.0;
			dxdt[0] = x[1];
			dxdt[1] = -9.81 / Length * sin(x[0]);
		}
	);

	// setting parameters for our simulation
	Real	t1 = 0.0, t2 = 10.0;
	int   expectNumSteps = 1000;							// means our step = 10 / 1000 = 0.001
	Real	initAngle = Utils::DegToRad(30);		// initial angle set to 30 degrees
	Vector<Real>	initCond{ initAngle, 0.0 };

	// solving and getting solutions
	ODESystemFixedStepSolver	fixedSolver(pendSys, StepCalculators::EulerStepCalc);
	ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);

	Vector<Real> sol_x{ sol.getTValues() }, sol_y1{ sol.getXValues(0) }, sol_y2{ sol.getXValues(1) };

	// console visualization
	std::cout << "\n\n****  Euler method - fixed stepsize  **********\n";
	std::vector<ColDesc> vecNames{ ColDesc("t", 11, 2, 'F'), ColDesc("angle", 15, 8, 'F'), ColDesc("ang.vel.", 15, 8, 'F') };
	std::vector<Vector<Real>*>	vecVals{ &sol_x, &sol_y1, &sol_y2 };

	VerticalVectorPrinter	vvp(vecNames, vecVals);

	vvp.Print();

	// visualizing solutions
	PolynomInterpRealFunc 	solAngle = sol.getSolutionAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solAngleVel = sol.getSolutionAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(solAngle, "Pendulum - angle in time",
		t1, t2, 200, "pendulum_angle_euler.txt");

	Visualizer::VisualizeRealFunction(solAngleVel, "Pendulum - ang.vel. in time",
		t1, t2, 200, "pendulum_ang_vel_euler.txt");

	// shown together as multi-real function
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solAngle, & solAngleVel },
		"Pendulum - both variables", { "Angle", "Ang.vel." },
		t1, t2, 200,
		"pendulum_multi_real_func_euler.txt");

	// visualized together as ODE system solution
	Visualizer::VisualizeODESysSolAsMultiFunc(sol, "Pendulum - Euler method", std::vector<std::string>{"angle", "angle.vel."},
		"pendulum_euler.txt");
}

void Demo_SimplePendulum()
{
	ODESystem pendSys(2, [](Real t, const Vector<Real>& x, Vector<Real>& dxdt)
		{
			Real Length = 1.0;
			dxdt[0] = x[1];
			dxdt[1] = -9.81 / Length * sin(x[0]);
		}
	);

	Real	pendulumLen = 1.0;
	PendulumODE		sys = PendulumODE(pendulumLen);

	Real	t1 = 0.0, t2 = 10.0;
	int   expectNumSteps = 100;
	Real	minSaveInterval = (t2 - t1) / expectNumSteps;

	// initial angle set to 30 degrees in radians, and no pushing at the start
	Real	initAngle = Utils::DegToRad(30), initAngVel = 0.0;
	Vector<Real>	initCond{ initAngle, initAngVel };

	ODESystemFixedStepSolver	fixedSolver(pendSys, StepCalculators::RK4_Basic);
	ODESystemSolution solFixed = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);

	ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(sys);
	ODESystemSolution solAdapt = adaptSolver.integrate(initCond, t1, t2, minSaveInterval / 1.21, 1e-06, 0.01);

	Vector<Real> x_fixed{ solFixed.getTValues() }, y1_fixed{ solFixed.getXValues(0) }, y2_fixed{ solFixed.getXValues(1) };
	Vector<Real> x_adapt{ solAdapt.getTValues() }, y1_adapt{ solAdapt.getXValues(0) }, y2_adapt{ solAdapt.getXValues(1) };

	std::cout << "\n\n****  Runge-Kutta 4th order - fixed stepsize  **********  Runge-Kutta 4th order - adaptive stepper  ****\n";
	std::vector<ColDesc>				vecNames{ ColDesc("t", 11, 2, 'F'), ColDesc("angle", 15, 8, 'F'), ColDesc("ang.vel.", 15, 8, 'F'),
																				ColDesc("t", 22, 2, 'F'), ColDesc("angle", 12, 8, 'F'), ColDesc("ang.vel.", 12, 8, 'F') };
	std::vector<Vector<Real>*>	vecVals{ &x_fixed, &y1_fixed, &y2_fixed,
																			 &x_adapt, &y1_adapt, &y2_adapt };
	VerticalVectorPrinter	vvp(vecNames, vecVals);

	//vvp.Print();

	// getting solutions as polynomials
	PolynomInterpRealFunc 	angleFixedSol = solFixed.getSolutionAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	angleVelFixedSol = solFixed.getSolutionAsPolyInterp(1, 3);

	PolynomInterpRealFunc 	angleAdaptSol = solAdapt.getSolutionAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	angleVelAdaptSol = solAdapt.getSolutionAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(angleAdaptSol, "Pendulum - angle in time",
		t1, t2, 200, "pendulum_angle.txt");

	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &angleAdaptSol, & angleVelAdaptSol },
		"Pendulum - both variables", { "Angle", "Ang.vel" },
		t1, t2, 200,
		"pendulum_multi_real_func.txt");

	Real periodLinear = sys.calcPeriodLinearized();
	Real exactPeriod = sys.calcExactPeriod(initAngle);
	// extracting period T from solutions
	Real simulPeriodFixed = RealFunctionAnalyzer(angleFixedSol).calcRootsPeriod(t1, t2, 1000);
	Real simulPeriodAdapt = RealFunctionAnalyzer(angleAdaptSol).calcRootsPeriod(t1, t2, 1000);

	std::cout << "Pendulum period approx. linear : " << periodLinear << std::endl;
	std::cout << "Pendulum period analytic exact : " << exactPeriod << std::endl;
	std::cout << "Pendulum period RK4 fixed step : " << 2 * simulPeriodFixed << std::endl;
	std::cout << "Pendulum period RK4 adapt.step : " << 2 * simulPeriodAdapt << std::endl;

	// Comparing periods for different initial angles
	Vector<Real> initAnglesDeg{ 1, 2, 5, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0 };
	Vector<Real> initAngles(initAnglesDeg.size());
	for (int i = 0; i < initAngles.size(); i++)
		initAngles[i] = Utils::DegToRad(initAnglesDeg[i]);

	TablePrinter<double, double> data("Angle", 5, 1,
		{ "Linear ", "Exact  ", "Diff %", "Fix.step.sim", "Diff %", "Adapt.step.sim", "Diff % " },
		{ {10,6,'F'},	{10,6,'F'}, {7,3,'F'}, {13,7,'F'}, {8,5,'F'}, {13,7,'F'}, {8,5,'F'} });

	for (int i = 0; i < initAngles.size(); i++)
	{
		initCond[0] = initAngles[i];

		// solving with fixed and adaptive solver
		ODESystemSolution			solF = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);
		ODESystemSolution			solA = adaptSolver.integrate(initCond, t1, t2, minSaveInterval, 1e-06, 0.01);

		// calculating period from solutions
		PolynomInterpRealFunc solFInterp = solF.getSolutionAsPolyInterp(0, 3);
		PolynomInterpRealFunc solAInterp = solA.getSolutionAsPolyInterp(0, 3);
		RealFunctionAnalyzer fa3(solFInterp);
		RealFunctionAnalyzer fa4(solAInterp);
		Real periodSimulFixed = fa3.calcRootsPeriod(t1, t2, 1000);
		Real periodSimulAdapt = fa4.calcRootsPeriod(t1, t2, 1000);

		// analytical formulas
		Real periodLin = sys.calcPeriodLinearized();
		Real periodExact = sys.calcExactPeriod(initAngles[i]);
		Real diffPercent = Abs(periodLin - periodExact) / periodExact * 100;
		Real diffPercentFixed = Abs(2 * periodSimulFixed - periodExact) / periodExact * 100;
		Real diffPercentAdapt = Abs(2 * periodSimulAdapt - periodExact) / periodExact * 100;

		data.addRow(initAnglesDeg[i], { periodLin, periodExact, diffPercent,
																		2 * periodSimulFixed, diffPercentFixed, 2 * periodSimulAdapt, diffPercentAdapt });
	}

	data.Print();

	// Comparing number of adaptive steps needed, for given accuracy
	Vector<Real> acc{ 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10 };
	Vector<Real> numSteps(acc.size());

	std::cout << "\nAccuracy   Num steps   OK steps   Bad steps" << std::endl;
	for (int i = 0; i < acc.size(); i++)
	{
		ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(sys);
		ODESystemSolution solAdapt = adaptSolver.integrate(initCond, t1, t2, minSaveInterval, acc[i], 0.01);
		numSteps[i] = solAdapt.getTotalNumSteps();

		std::cout << std::setw(7) << acc[i] << "        " << std::setw(3) << numSteps[i] << "        "
			<< std::setw(3) << solAdapt.getNumStepsOK() << "        "
			<< std::setw(3) << solAdapt.getNumStepsBad() << std::endl;
	}
}

void Demo_DumpedPendulum()
{
	DumpedPendulumODE odeSys(1.0, 0.5);

	Real	t1 = 0.0, t2 = 20.0;
	int   expectNumSteps = 1000;
	Real	initAngle = 0.5;
	Vector<Real>	initCond{ initAngle, 0.0 };

	ODESystemFixedStepSolver	fixedSolver(odeSys, StepCalculators::RK4_Basic);
	ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);

	Vector<Real> x_fixed = sol.getTValues();
	Vector<Real> y1_fixed = sol.getXValues(0);
	Vector<Real> y2_fixed = sol.getXValues(1);

	std::cout << "\n\n****  Runge-Kutta 4th order - fixed stepsize  **********\n";
	std::vector<ColDesc>				vecNames{ ColDesc("t", 11, 2, 'F'), ColDesc("angle", 15, 8, 'F'), ColDesc("ang.vel.", 15, 8, 'F') };
	std::vector<Vector<Real>*>	vecVals{ &x_fixed, &y1_fixed, &y2_fixed };
	VerticalVectorPrinter	vvp(vecNames, vecVals);
	//vvp.Print();

	PolynomInterpRealFunc 	solInterpY1 = sol.getSolutionAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solInterpY2 = sol.getSolutionAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(solInterpY1, "Dumped Pendulum - angle in time",
		t1, t2, 200, "dumped_pendulum_angle.txt");
	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solInterpY1, & solInterpY2 },
		"Dumped Pendulum - both variables", { "Angle", "Ang.vel." },
		t1, t2, 200,
		"dumped_pendulum_multi_real_func.txt");
}

void Demo_ForcedPendulum()
{
	ForcedPendulumODE forcedPen(1.0, 0.2, 1);

	Real	t1 = 0.0, t2 = 30.0;
	int   expectNumSteps = 3000;
	Real	minSaveInterval = (t2 - t1) / expectNumSteps;
	Real	initAngle = 0.5;
	Vector<Real>	initCond{ initAngle, 0.0 };

	ODESystemSolver<RK5_CashKarp_Stepper> odeSolver(forcedPen);
	ODESystemSolution sol = odeSolver.integrate(initCond, t1, t2, minSaveInterval, 1e-08, 0.001);

	Vector<Real> t_vals = sol.getTValues();
	Vector<Real> y1_fixed = sol.getXValues(0);
	Vector<Real> y2_fixed = sol.getXValues(1);

	std::cout << "\n\n****  Runge-Kutta 4th order - fixed stepsize  **********\n";
	std::vector<ColDesc>				vecNames{ ColDesc("t", 11, 2, 'F'), ColDesc("angle", 15, 8, 'F'), ColDesc("ang.vel.", 15, 8, 'F') };
	std::vector<Vector<Real>*>	vecVals{ &t_vals, &y1_fixed, &y2_fixed };
	VerticalVectorPrinter	vvp(vecNames, vecVals);
	//vvp.Print();

	PolynomInterpRealFunc 	solFixedPolyInterp0 = sol.getSolutionAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solFixedPolyInterp1 = sol.getSolutionAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(solFixedPolyInterp0, "Forced Pendulum - angle in time",
		t1, t2, 1000, "forced_pendulum_angle.txt");
	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solFixedPolyInterp0, & solFixedPolyInterp1 },
		"Forced Pendulum - both variables", { "Angle", "Ang.vel." },
		t1, t2, 1000,
		"forced_pendulum_multi_real_func.txt");

	Matrix<Real> curve_points(t_vals.size(), 2);
	for (int i = 0; i < t_vals.size(); i++)
	{
		curve_points(i, 0) = y1_fixed[i];
		curve_points(i, 1) = y2_fixed[i];
	}

	SplineInterpParametricCurve<2> phaseSpaceTrajectory(0.0, 1.0, curve_points);

	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory, "Forced pendulum - phase space trajectory",
		0.0, 1.0, t_vals.size(), "forced_pendulum_phase_space.txt");
}

int main()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 6 - Pendulum                        ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Demo_SimplePendulum();
	//Demo_SimplePendulumEuler();
	//Demo_DumpedPendulum();
	//Demo_ForcedPendulum();
}