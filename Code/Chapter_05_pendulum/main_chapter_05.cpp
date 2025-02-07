#include "MMLBase.h"

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Vector.h"

#include "tools/ConsolePrinter.h"
#include "tools/Visualizer.h"
#include "tools/Serializer.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "algorithms/RootFinding.h"



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

double calculatePendulumPeriod(double length, double initialAngle) {
	const double g = 9.81; // acceleration due to gravity
	double k = std::sin(initialAngle / 2);
	double ellipticIntegral = std::comp_ellint_1(k);
	double period = 4 * std::sqrt(length / g) * ellipticIntegral;
	return period;
}

void Demo_ExactPendulum()
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

	Real	t1 = 0.0, t2 = 1000.0;
	int   expectNumSteps = 10000;
	Real initAngle = 0.5;
	Vector<Real>	initCond{ initAngle, 0.0 };

	ODESystemFixedStepSolver	sysSolver(pendSys, StepCalculators::RK4_Basic);
	ODESystemSolution sol = sysSolver.integrate(initCond, t1, t2, expectNumSteps);

	const double h1 = 0.01, hmin = 0.0;
	ODESystemSolver<RK4_CashKarp_Stepper> solver(sys);
	Real minSaveInterval = (t2 - t1) / expectNumSteps;
	ODESystemSolution sol3 = solver.integrate(initCond, t1, t2, minSaveInterval, 1e-06, h1, hmin);

	Vector<Real> x_fixed = sol.getXValues();
	Vector<Real> y1_fixed = sol.getYValues(0);
	Vector<Real> y2_fixed = sol.getYValues(1);

	Vector<Real> x_adapt = sol3.getXValues();
	Vector<Real> y1_adapt = sol3.getYValues(0);
	Vector<Real> y2_adapt = sol3.getYValues(1);

	std::cout << "\n\n****  Runge-Kutta 4th order - fixed stepsize  **********  Runge-Kutta 4th order - adaptive stepper  ****\n";
	std::vector<ColDesc>				vecNames{ ColDesc("t", 11, 2, 'F'), ColDesc("angle", 15, 8, 'F'), ColDesc("ang.vel.", 15, 8, 'F'),
																				ColDesc("t", 22, 2, 'F'), ColDesc("angle", 12, 8, 'F'), ColDesc("ang.vel.", 12, 8, 'F') };
	std::vector<Vector<Real>*>	vecVals{ &x_fixed, &y1_fixed, &y2_fixed,
																			 &x_adapt, &y1_adapt, &y2_adapt };
	VerticalVectorPrinter	vvp(vecNames, vecVals);

	//vvp.Print();

	PolynomInterpRealFunc 	solPolyInterp0 = sol3.getSolutionAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solPolyInterp1 = sol3.getSolutionAsPolyInterp(1, 3);
	SplineInterpRealFunc		solSplineInterp0 = sol3.getSolutionAsSplineInterp(0);

	SplineInterpParametricCurve<2> spline(sol3.getYValues());

	Visualizer::VisualizeRealFunction(solSplineInterp0, "Pendulum - angle in time", 0.0, 10.0, 200, "pendulum_angle.txt");

	// shown together
	//Visualizer::VisualizeMultiRealFunction({ &solPolyInterp0, &solPolyInterp1 }, 
	//																			"Pendulum - both variables", 0.0, 10.0, 200, "pendulum_multi_real_func.txt");

	// extracting period T from solutions
	int numFoundRoots;
	Vector<Real> root_brack_x1(10), root_brack_x2(10);
	FindRootBrackets(solPolyInterp0, t1, t2, 4000, root_brack_x1, root_brack_x2, numFoundRoots);

	std::cout << "Found " << numFoundRoots << " roots" << std::endl;
	Vector<Real> roots(numFoundRoots);
	for (int i = 0; i < numFoundRoots; i++)
	{
		roots[i] = FindRootBisection(solPolyInterp0, root_brack_x1[i], root_brack_x2[i], 1e-7);

		std::cout << "Root " << i << " : " << roots[i];
		if (i > 0)
			std::cout << "      Difference: " << roots[i] - roots[i - 1];
		std::cout << std::endl;
	}
	Real period = 2.0 * Constants::PI * sqrt(pendulumLen / 9.81);
	std::cout << "Period of pendulum: " << period / 2 << std::endl;

	// calculate exact period for pendulum of length L, and given initial angle phi
	Real exactPeriod = calculatePendulumPeriod(pendulumLen, initAngle);
	std::cout << "Exact period of pendulum: " << exactPeriod / 2.0 << std::endl;



	// comparison for different initial angles

	// usporediti za različite eps 1e-2, 1e-3, ... , 1e-8 kako izgledauju stepovi i vizualizirati ih
}

int main()
{
	Demo_ExactPendulum();

  return 0;
}

