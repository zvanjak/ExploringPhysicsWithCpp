#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/BaseUtils.h"
#include "base/Vector.h"

#include "tools/ConsolePrinter.h"
#include "tools/Visualizer.h"
#include "tools/Serializer.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "Mechanics/Projectiles.h"
#endif

using namespace MML;
using namespace MPL;


void Projectile_in_vacuum()
{
	ProjectileInVacuumODE projSys;

	Real angle = Utils::DegToRad(45);		// launching angle for our projectile
	Real initVel = 700;									// initial velocity in m/s
	Real initHeight = 0;								// initial height in m
	Vector<Real> initCond = projSys.getInitCond(angle, initHeight, initVel);

	Real t1 = 0, t2 = 1.01 * projSys.TimeOfFlight(angle, initHeight, initVel); // time of flight for the projectile motion in vacuum
	int expectNumSteps = 200;			// expected number of steps in the solution

	ODESystemFixedStepSolver	odeSolver(projSys, StepCalculators::EulerStepCalc);
	ODESystemSolution sol = odeSolver.integrate(initCond, t1, t2, expectNumSteps);

	//ODESystemSolver<RK5_CashKarp_Stepper> odeSolver(projSys);
	//ODESystemSolution sol = odeSolver.integrate(initCond, t1, t2, 0.01, 1e-08, 0.01);

	Vector<Real> sol_t{ sol.getTValues() }, sol_x{ sol.getXValues(0) }, sol_y{ sol.getXValues(1) }, sol_vx{ sol.getXValues(2) }, sol_vy{ sol.getXValues(3) };

	// Visualizing the solution
	PolynomInterpRealFunc 	sol_Y_of_t = sol.getSolutionAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(sol_Y_of_t, "Projectile in vacuum - height in time",
		t1, t2, 100, "projectile_height_in_time.txt");

	LinearInterpRealFunc s1(sol_x, sol_y);			// forming y(x) function

	//Real x2 = projSys.CalcRange(angle, initHeight, initVel); // range of the projectile motion in vacuum
	//Visualizer::VisualizeRealFunction(s1, "Projectile in vacuum - y(x)", 0, x2, 100, "projectile_y(x).txt");

	//Visualizer::VisualizeODESysSolAsMultiFunc(sol, "Projectile in vacuum - all variables", std::vector<std::string>{"x", "y", "vx", "vy"},
	//																					"projectile_launch_ode_sol.txt");

	// visualize solution for a given set of initial angles
	Vector<Real> angles{ Utils::DegToRad(20), Utils::DegToRad(30), Utils::DegToRad(45), Utils::DegToRad(50), Utils::DegToRad(55), Utils::DegToRad(60), Utils::DegToRad(65) };
	Vector<Real> listT2;		// time of flight for each angle
	for (auto angle : angles) {
		Real flightTime = projSys.TimeOfFlight(angle, initHeight, initVel); // time of flight for the projectile motion in vacuum
		listT2.push_back(flightTime);
	}

	Vector<LinearInterpRealFunc> solutions;
	std::vector<IRealFunction*> solYs;
	std::vector<std::string> strAngles;
	for (int i = 0; i < angles.size(); i++)
	{
		Vector<Real> initCond = projSys.getInitCond(angles[i], initHeight, initVel);
		ODESystemSolution sol = odeSolver.integrate(initCond, t1, listT2[i], expectNumSteps);
		//ODESystemSolution sol = odeSolver.integrate(initCond, t1, listT2[i], 0.01, 1e-06, 0.01);

		Vector<Real> sol_x1{ sol.getXValues(0) }, sol_y1{ sol.getXValues(1) };
		LinearInterpRealFunc y_of_x(sol_x1, sol_y1); // y(x) function

		solutions.push_back(y_of_x);
	}

	// get vector of RealFunctions for y-component from each solution
	// init with pointers from solutions vector
	int i = 0;
	for (auto& fun : solutions)
	{
		solYs.push_back(&fun);
		strAngles.push_back("Angle - " + std::to_string(Utils::RadToDeg(angles[i++])));
	}

	// now visualize all solutions
	Visualizer::VisualizeMultiRealFunction(solYs,
		"Projectile in vacuum - different angles", strAngles,
		0, 50000, 100, "projectile_vacuum_multi.txt");
}

void Projectile_with_air_resistance()
{
	ProjectileWithAirResistanceODE projSys(4e-5); // air resistance coefficient

	Real angle = Utils::DegToRad(45);		// launching angle for our projectile
	Real initVel = 700;									// initial velocity in m/s
	Real initHeight = 0;								// initial height in m
	Vector<Real> initCond = projSys.getInitCond(angle, initHeight, initVel);

	Real t1 = 0, t2 = 100;
	int expectNumSteps = 300;			// expected number of steps in the solution

	ODESystemFixedStepSolver	odeSolver(projSys, StepCalculators::EulerStepCalc);
	ODESystemSolution sol = odeSolver.integrate(initCond, t1, t2, expectNumSteps);

	//ODESystemSolver<RK5_CashKarp_Stepper> odeSolver(projSys);
	//ODESystemSolution sol = odeSolver.integrate(initCond, t1, t2, 0.01, 1e-08, 0.01);

	Vector<Real> sol_t{ sol.getTValues() }, sol_x{ sol.getXValues(0) }, sol_y{ sol.getXValues(1) }, sol_vx{ sol.getXValues(2) }, sol_vy{ sol.getXValues(3) };

	// Visualizing the solution
	PolynomInterpRealFunc 	solY = sol.getSolutionAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(solY, "Projectile height in time with air resistance",
		t1, t2, 300, "projectile_height_in_time_with_air_resistance.txt");
	Visualizer::VisualizeODESysSolAsMultiFunc(sol, "Projectile launch with air resistance - Euler method", std::vector<std::string>{"x", "y", "vx", "vy"},
		"projectile_launch_ode_sol_with_air_resistance.txt");

	// visualize solution for a given set of initial angles
	Vector<Real> angles{ Utils::DegToRad(20), Utils::DegToRad(30), Utils::DegToRad(45), Utils::DegToRad(50), Utils::DegToRad(55), Utils::DegToRad(60), Utils::DegToRad(65) };

	Vector<Real> listT2; // time of flight for each angle
	for (auto angle : angles)
	{
		Real flightTime = t2; // projSys.TimeOfFlight(angle, initHeight, initVel); // time of flight for the projectile motion in vacuum
		listT2.push_back(flightTime);
	}
	//	Vector<Real> listT2{ 50, 70, 100, 100, 100, 90, 90 }; // time of flight for each angle

	Vector<LinearInterpRealFunc> solutions;
	std::vector<IRealFunction*> solYs;
	std::vector<std::string> strAngles;
	for (int i = 0; i < angles.size(); i++)
	{
		Vector<Real> initCond = projSys.getInitCond(angles[i], initHeight, initVel);
		ODESystemSolution sol = odeSolver.integrate(initCond, t1, listT2[i], expectNumSteps);
		//ODESystemSolution sol = odeSolver.integrate(initCond, t1, listT2[i], 0.01, 1e-06, 0.01);

		Vector<Real> sol_x1{ sol.getXValues(0) }, sol_y1{ sol.getXValues(1) };
		LinearInterpRealFunc y_of_x(sol_x1, sol_y1); // y(x) function

		solutions.push_back(y_of_x);
	}

	// get vector of RealFunctions for y-component from each solution
	// init with pointers from solutions vector
	int i = 0;
	for (auto& fun : solutions)
	{
		solYs.push_back(&fun);
		strAngles.push_back("Angle - " + std::to_string(Utils::RadToDeg(angles[i++])));
	}

	// now visualize all solutions
	Visualizer::VisualizeMultiRealFunction(solYs,
		"Projectile with air resistance - different angles", strAngles,
		0, 22000, 300, "projectile_air_res_multi.txt");
}

void Projectile_with_changing_air_density()
{
	ProjectileChangingAirDensityODE projSys(4e-5); // air resistance coefficient

	Real angle = Utils::DegToRad(45);		// launching angle for our projectile
	Real initVel = 700;									// initial velocity in m/s
	Real initHeight = 0;								// initial height in m
	Vector<Real> initCond = projSys.getInitCond(angle, initHeight, initVel);

	Real t1 = 0, t2 = 100;
	int expectNumSteps = 300;			// expected number of steps in the solution

	ODESystemFixedStepSolver	odeSolver(projSys, StepCalculators::EulerStepCalc);
	ODESystemSolution sol = odeSolver.integrate(initCond, t1, t2, expectNumSteps);

	Vector<Real> sol_t{ sol.getTValues() }, sol_x{ sol.getXValues(0) }, sol_y{ sol.getXValues(1) }, sol_vx{ sol.getXValues(2) }, sol_vy{ sol.getXValues(3) };

	// Visualizing the solution
	PolynomInterpRealFunc 	solY = sol.getSolutionAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(solY, "Projectile height in time with air resistance dependent on air density",
		t1, t2, 300, "projectile_height_in_time_with_air_resistance_air_density.txt");
	Visualizer::VisualizeODESysSolAsMultiFunc(sol, "Projectile launch with air resistance - Euler method", std::vector<std::string>{"x", "y", "vx", "vy"},
		"projectile_height_in_time_with_air_resistance_air_density.txt");

	// visualize solution for a given set of initial angles
	Vector<Real> angles{ Utils::DegToRad(20), Utils::DegToRad(30), Utils::DegToRad(45), Utils::DegToRad(50), Utils::DegToRad(55), Utils::DegToRad(60), Utils::DegToRad(65) };

	Vector<Real> listT2; // time of flight for each angle
	for (auto angle : angles)
	{
		Real flightTime = t2; // projSys.TimeOfFlight(angle, initHeight, initVel); // time of flight for the projectile motion in vacuum
		listT2.push_back(flightTime);
	}
	//	Vector<Real> listT2{ 50, 70, 100, 100, 100, 90, 90 }; // time of flight for each angle

	Vector<LinearInterpRealFunc> solutions;
	std::vector<IRealFunction*> solYs;
	std::vector<std::string> strAngles;
	for (int i = 0; i < angles.size(); i++)
	{
		Vector<Real> initCond = projSys.getInitCond(angles[i], initHeight, initVel);
		ODESystemSolution sol = odeSolver.integrate(initCond, t1, listT2[i], expectNumSteps);
		//ODESystemSolution sol = odeSolver.integrate(initCond, t1, listT2[i], 0.01, 1e-06, 0.01);

		Vector<Real> sol_x1{ sol.getXValues(0) }, sol_y1{ sol.getXValues(1) };
		LinearInterpRealFunc y_of_x(sol_x1, sol_y1); // y(x) function

		solutions.push_back(y_of_x);
	}

	// get vector of RealFunctions for y-component from each solution
	// init with pointers from solutions vector
	int i = 0;
	for (auto& fun : solutions)
	{
		solYs.push_back(&fun);
		strAngles.push_back("Angle - " + std::to_string(Utils::RadToDeg(angles[i++])));
	}

	// now visualize all solutions
	Visualizer::VisualizeMultiRealFunction(solYs,
		"Projectile with air resistance dependent on air density - different angles", strAngles,
		0, 22000, 300, "projectile_air_res_air_density_multi.txt");
}

void Projectile_with_drag_coeff_dependent_on_speed()
{
	BaseballWithDragCoeffDependentOnSpeedODE projSys; // air resistance coefficient

	Real angle = Utils::DegToRad(45);		// launching angle for our projectile
	Real initVel = 35;									// initial velocity in m/s
	Real initHeight = 0;								// initial height in m
	Vector<Real> initCond = projSys.getInitCond(angle, initHeight, initVel);

	Real t1 = 0, t2 = 5;
	int expectNumSteps = 300;			// expected number of steps in the solution

	ODESystemFixedStepSolver	odeSolver(projSys, StepCalculators::EulerStepCalc);
	ODESystemSolution sol = odeSolver.integrate(initCond, t1, t2, expectNumSteps);

	Vector<Real> sol_t{ sol.getTValues() }, sol_x{ sol.getXValues(0) }, sol_y{ sol.getXValues(1) }, sol_vx{ sol.getXValues(2) }, sol_vy{ sol.getXValues(3) };

	// Visualizing the solution
	PolynomInterpRealFunc 	solY = sol.getSolutionAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(solY, "Baseball height in time with drag.coeff. dependent on speed",
		t1, t2, 300, "baseball_height_in_time.txt");
	Visualizer::VisualizeODESysSolAsMultiFunc(sol, "Projectile launch with air resistance - Euler method", std::vector<std::string>{"x", "y", "vx", "vy"},
		"baseball_height_in_time.txt");

	// visualize solution for a given set of initial angles
	Vector<Real> angles{ Utils::DegToRad(20), Utils::DegToRad(30), Utils::DegToRad(45), Utils::DegToRad(50), Utils::DegToRad(55), Utils::DegToRad(60), Utils::DegToRad(65) };

	Vector<Real> listT2; // time of flight for each angle
	for (auto angle : angles)
	{
		Real flightTime = t2; // projSys.TimeOfFlight(angle, initHeight, initVel); // time of flight for the projectile motion in vacuum
		listT2.push_back(flightTime);
	}
	//	Vector<Real> listT2{ 50, 70, 100, 100, 100, 90, 90 }; // time of flight for each angle

	Vector<LinearInterpRealFunc> solutions;
	std::vector<IRealFunction*> solYs;
	std::vector<std::string> strAngles;
	for (int i = 0; i < angles.size(); i++)
	{
		Vector<Real> initCond = projSys.getInitCond(angles[i], initHeight, initVel);
		ODESystemSolution sol = odeSolver.integrate(initCond, t1, listT2[i], expectNumSteps);
		//ODESystemSolution sol = odeSolver.integrate(initCond, t1, listT2[i], 0.01, 1e-06, 0.01);

		Vector<Real> sol_x1{ sol.getXValues(0) }, sol_y1{ sol.getXValues(1) };
		LinearInterpRealFunc y_of_x(sol_x1, sol_y1); // y(x) function

		solutions.push_back(y_of_x);
	}

	// get vector of RealFunctions for y-component from each solution
	// init with pointers from solutions vector
	int i = 0;
	for (auto& fun : solutions)
	{
		solYs.push_back(&fun);
		strAngles.push_back("Angle - " + std::to_string(Utils::RadToDeg(angles[i++])));
	}

	// now visualize all solutions
	Visualizer::VisualizeMultiRealFunction(solYs,
		"Projectile with drag.coeff dependent on speed - different angles", strAngles,
		0, 100, 300, "projectile_air_res_air_density_multi.txt");
}

void Compare_solutions_for_set_of_angles(const IProjectileODE& sys1, std::string label1, const IProjectileODE& sys2, std::string label2,
	Real initVel, Real initHeight, int numSteps, Vector<Real> angles, Vector<Real> timeOfFlight)
{

	std::vector<IRealFunction*> solYs;
	std::vector<std::string> strAngles;
	Real maxX = 0.0;

	// first, solve the system in vacuum
	Vector<LinearInterpRealFunc> solutionsVacuum;
	ODESystemFixedStepSolver	odeSolver(sys1, StepCalculators::EulerStepCalc);
	for (int i = 0; i < angles.size(); i++)
	{
		Vector<Real> initCond = sys1.getInitCond(angles[i], initHeight, initVel);
		ODESystemSolution sol = odeSolver.integrate(initCond, 0.0, timeOfFlight[i], numSteps);

		Vector<Real> sol_x1{ sol.getXValues(0) }, sol_y1{ sol.getXValues(1) };
		LinearInterpRealFunc y_of_x(sol_x1, sol_y1); // y(x) function

		if (y_of_x.MaxX() > maxX)
			maxX = y_of_x.MaxX(); // find maximum x value for all solutions

		solutionsVacuum.push_back(y_of_x);
	}
	int i = 0;
	for (auto& fun : solutionsVacuum)
	{
		solYs.push_back(&fun);
		strAngles.push_back(label1 + ", angle - " + std::to_string(Utils::RadToDeg(angles[i++])));
	}

	// now, solve the system with air resistance
	Vector<LinearInterpRealFunc> solutionsAir;

	ODESystemFixedStepSolver	odeSolverAir(sys2, StepCalculators::EulerStepCalc);
	for (int i = 0; i < angles.size(); i++)
	{
		Vector<Real> initCond = sys2.getInitCond(angles[i], initHeight, initVel);
		ODESystemSolution sol = odeSolverAir.integrate(initCond, 0.0, timeOfFlight[i], numSteps);

		Vector<Real> sol_x1{ sol.getXValues(0) }, sol_y1{ sol.getXValues(1) };
		LinearInterpRealFunc y_of_x(sol_x1, sol_y1); // y(x) function

		solutionsAir.push_back(y_of_x);
	}
	i = 0;
	for (auto& fun : solutionsAir)
	{
		solYs.push_back(&fun);
		strAngles.push_back(label2 + ", angle - " + std::to_string(Utils::RadToDeg(angles[i++])));
	}

	// now visualize all solutions
	Visualizer::VisualizeMultiRealFunction(solYs,
		"Comparing " + label1 + " and " + label2 + " solutions", strAngles,
		0, maxX * 1.1, 100, "projectile_comparison_" + label1 + "_vs_" + label2 + ".txt");
}

void Compare_solutions_vacuum_vs_air_resistance()
{
	ProjectileInVacuumODE						projSysVacuum;
	ProjectileWithAirResistanceODE	projSysAir(0.01); // air resistance coefficient

	Real initVel = 20;									// initial velocity in m/s
	Real initHeight = 0.0;
	int	 expectNumSteps = 200;					// expected number of steps in the solution

	// compare solutions for a given set of initial angles
	Vector<Real> angles{ Utils::DegToRad(20), Utils::DegToRad(30), Utils::DegToRad(40),Utils::DegToRad(45), Utils::DegToRad(50), Utils::DegToRad(55),Utils::DegToRad(60) };
	Vector<Real> listT2;		// time of flight for each angle
	for (auto angle : angles) {
		Real flightTime = projSysVacuum.TimeOfFlight(angle, initHeight, initVel); // time of flight for the projectile motion in vacuum
		listT2.push_back(flightTime);
	}

	Compare_solutions_for_set_of_angles(projSysVacuum, "Vacuum", projSysAir, "Air resist.", initVel, initHeight, expectNumSteps, angles, listT2);
}

void Compare_solutions_air_resistance_vs_changing_air_density()
{
	ProjectileWithAirResistanceODE	projSysAir(4e-5); // air resistance coefficient
	ProjectileChangingAirDensityODE projSysChangingDensity(4e-5); // air resistance coefficient

	Real initVel = 500;									// initial velocity in m/s
	Real initHeight = 0.0;
	int	 expectNumSteps = 200;					// expected number of steps in the solution

	// compare solutions for a given set of initial angles
	Vector<Real> angles{ Utils::DegToRad(20), Utils::DegToRad(30), Utils::DegToRad(40),Utils::DegToRad(45), Utils::DegToRad(50), Utils::DegToRad(55),Utils::DegToRad(60) };
	Vector<Real> listT2;		// time of flight for each angle
	for (auto angle : angles) {
		Real flightTime = 100; //  projSysVacuum.TimeOfFlight(angle, initHeight, initVel); // time of flight for the projectile motion in vacuum
		listT2.push_back(flightTime);
	}

	Compare_solutions_for_set_of_angles(projSysAir, "Air resist", projSysChangingDensity, "Chang.air dens.", initVel, initHeight, expectNumSteps, angles, listT2);
}

void Compare_solutions_air_resistance_vs_drag_dependent_on_speed()
{
	ProjectileWithAirResistanceODE						projSysAir(4e-5); // air resistance coefficient
	BaseballWithDragCoeffDependentOnSpeedODE	projSysDragDepOnSpeed;

	Real initVel = 35;									// initial velocity in m/s
	Real initHeight = 0.0;
	int	 expectNumSteps = 300;					// expected number of steps in the solution

	// compare solutions for a given set of initial angles
	Vector<Real> angles{ Utils::DegToRad(20), Utils::DegToRad(30), Utils::DegToRad(40),Utils::DegToRad(45), Utils::DegToRad(50), Utils::DegToRad(55),Utils::DegToRad(60) };
	Vector<Real> listT2;		// time of flight for each angle
	for (auto angle : angles) {
		Real flightTime = 10; //  projSysVacuum.TimeOfFlight(angle, initHeight, initVel); // time of flight for the projectile motion in vacuum
		listT2.push_back(flightTime);
	}

	Compare_solutions_for_set_of_angles(projSysAir, "Air resist", projSysDragDepOnSpeed, "Drag.dep.on speed.", initVel, initHeight, expectNumSteps, angles, listT2);
}

int main()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****               EXAMPLE 7 - Throwing things in the air          ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	//Projectile_in_vacuum();
	//Projectile_with_air_resistance();
	//Compare_solutions_vacuum_vs_air_resistance();
	//Compare_solutions_air_resistance_vs_changing_air_density();
	Compare_solutions_air_resistance_vs_drag_dependent_on_speed();

	//Projectile_with_changing_air_density();
	//Projectile_with_drag_coeff_dependent_on_speed();
}