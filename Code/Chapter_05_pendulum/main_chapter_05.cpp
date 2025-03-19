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


using namespace MML;
using namespace MML::RootFinding;

class PendulumODE : public IODESystem
{
	Real _Length;
public:
	PendulumODE(Real length) : _Length(length) {}

	int  getDim() const override { return 2; }
	void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
	{
		dxdt[0] = x[1];
		dxdt[1] = -9.81 / _Length * sin(x[0]);
	}
};

double calculatePendulumPeriod(double length, double initialAngle)
{
	const double g = 9.81;
	double k = std::sin(initialAngle / 2);
	double ellipticIntegral = std::comp_ellint_1(k);

	double period = 4 * std::sqrt(length / g) * ellipticIntegral;

	return period;
}

Real calcFunctionPeriod(const IRealFunction& func, Real t1, Real t2)
{
	int		numFoundRoots;
	Vector<Real> root_brack_x1(10), root_brack_x2(10);
	FindRootBrackets(func, t1, t2, 4000, root_brack_x1, root_brack_x2, numFoundRoots);

	Vector<Real> roots(numFoundRoots);
	Vector<Real> rootDiffs(numFoundRoots - 1);
	for (int i = 0; i < numFoundRoots; i++)
	{
		roots[i] = FindRootBisection(func, root_brack_x1[i], root_brack_x2[i], 1e-7);

		if (i > 0)
			rootDiffs[i - 1] = roots[i] - roots[i - 1];
	}

	return Statistics::Avg(rootDiffs);
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
	Real	initAngle = 0.5;
	Vector<Real>	initCond{ initAngle, 0.0 };

	ODESystemFixedStepSolver	fixedSolver(pendSys, StepCalculators::RK4_Basic);
	ODESystemSolution solFixed = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);

	ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(sys);
	ODESystemSolution solAdapt = adaptSolver.integrate(initCond, t1, t2, minSaveInterval / 1.21, 1e-06, 0.01);

	Vector<Real> x_fixed = solFixed.getTValues();
	Vector<Real> y1_fixed = solFixed.getXValues(0);
	Vector<Real> y2_fixed = solFixed.getXValues(1);

	Vector<Real> x_adapt = solAdapt.getTValues();
	Vector<Real> y1_adapt = solAdapt.getXValues(0);
	Vector<Real> y2_adapt = solAdapt.getXValues(1);

	std::cout << "\n\n****  Runge-Kutta 4th order - fixed stepsize  **********  Runge-Kutta 4th order - adaptive stepper  ****\n";
	std::vector<ColDesc>				vecNames{ ColDesc("t", 11, 2, 'F'), ColDesc("angle", 15, 8, 'F'), ColDesc("ang.vel.", 15, 8, 'F'),
																				ColDesc("t", 22, 2, 'F'), ColDesc("angle", 12, 8, 'F'), ColDesc("ang.vel.", 12, 8, 'F') };
	std::vector<Vector<Real>*>	vecVals{ &x_fixed, &y1_fixed, &y2_fixed,
																			 &x_adapt, &y1_adapt, &y2_adapt };
	VerticalVectorPrinter	vvp(vecNames, vecVals);

	vvp.Print();

	// getting solutions as polynomials
	PolynomInterpRealFunc 	solFixedPolyInterp0 = solFixed.getSolutionAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solFixedPolyInterp1 = solFixed.getSolutionAsPolyInterp(1, 3);

	PolynomInterpRealFunc 	solAdaptPolyInterp0 = solAdapt.getSolutionAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solAdaptPolyInterp1 = solAdapt.getSolutionAsPolyInterp(1, 3);

	SplineInterpRealFunc		solSplineInterp0 = solAdapt.getSolutionAsSplineInterp(0);

	Visualizer::VisualizeRealFunction(solAdaptPolyInterp0, "Pendulum - angle in time",
		0.0, 10.0, 200, "pendulum_angle.txt");

	// shown together
	Visualizer::VisualizeMultiRealFunction({ &solAdaptPolyInterp0, &solAdaptPolyInterp1 },
		"Pendulum - both variables", 0.0, 10.0, 200,
		"pendulum_multi_real_func.txt");

	// extracting period T from solutions
	Real periodLinear = 2.0 * Constants::PI * sqrt(pendulumLen / 9.81);
	std::cout << "Pendulum period approx. linear : " << periodLinear << std::endl;

	Real simulPeriodFixed = calcFunctionPeriod(solFixedPolyInterp0, t1, t2);
	std::cout << "Pendulum period RK4 fixed step : " << 2 * simulPeriodFixed << std::endl;

	Real simulPeriodAdapt = calcFunctionPeriod(solAdaptPolyInterp0, t1, t2);
	std::cout << "Pendulum period RK4 adapt.step : " << 2 * simulPeriodAdapt << std::endl;

	// calculate exact period for pendulum of length L, and given initial angle phi
	Real exactPeriod = calculatePendulumPeriod(pendulumLen, initAngle);
	std::cout << "Pendulum period analytic exact : " << exactPeriod << std::endl;


	// Comparing periods for different initial angles

	Vector<Real> initAnglesDeg{ 1, 2, 5, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0 };
	Vector<Real> initAngles(initAnglesDeg.size());				// convert to radians
	for (int i = 0; i < initAngles.size(); i++)
		initAngles[i] = initAnglesDeg[i] * Constants::PI / 180.0;

	Vector<Real> periods(initAngles.size());
	Real periodLin, periodSimulFixed, periodSimulAdapt, periodExact;

	std::cout << "\nAngle      Linear.      Exact     Diff %   Fix.step.sim Adapt.step.sim" << std::endl;
	for (int i = 0; i < initAngles.size(); i++)
	{
		initCond[0] = initAngles[i];

		// calculate period from fixed solution
		ODESystemSolution			solF = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);
		PolynomInterpRealFunc solFInterp = solF.getSolutionAsPolyInterp(0, 3);
		periodSimulFixed = calcFunctionPeriod(solFInterp, t1, t2);

		// calculate period from adaptive solution
		ODESystemSolution			solA = adaptSolver.integrate(initCond, t1, t2, minSaveInterval, 1e-06, 0.01);
		PolynomInterpRealFunc solAInterp = solA.getSolutionAsPolyInterp(0, 3);
		periodSimulAdapt = calcFunctionPeriod(solAInterp, t1, t2);

		// analytical formulas
		periodLin = 2.0 * Constants::PI * sqrt(pendulumLen / 9.81);
		periodExact = calculatePendulumPeriod(pendulumLen, initAngles[i]);

		double diffPercent = Abs(periodLin - periodExact) / periodExact * 100;
		std::cout << std::setw(2) << initAnglesDeg[i] << " deg:  " << periodLin << "    " 
			<< periodExact << "   " << std::setw(5) << diffPercent << "    " << 2 * periodSimulFixed << "    " << 2 * periodSimulAdapt << std::endl;
	}

	// Comparing number of adaptive steps needed, for given accuracy

	Vector<Real> acc{ 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8 };
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

int main()
{
	Demo_SimplePendulum();

  return 0;
}

