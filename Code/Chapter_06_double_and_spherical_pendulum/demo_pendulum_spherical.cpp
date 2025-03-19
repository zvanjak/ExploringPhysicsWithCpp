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


class SphericalPendulumODE : public IODESystem
{
	Real _l;
public:
	SphericalPendulumODE() : _l(1.0) {}
	SphericalPendulumODE(Real l) : _l(l) {}

	int  getDim() const override { return 4; }
	void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
	{
		Real tht = x[0];
		Real w1  = x[1];
		Real phi = x[2];
		Real w2  = x[3];

		Real g = 9.81;

		dxdt[0] = w1;
		dxdt[1] = sin(tht) * cos(tht) * POW2(w2) - g / _l * sin(tht);
		dxdt[2] = w2;
		dxdt[3] = -2.0 * w1 * w2 / tan(tht);
	}
};


void Demo_SphericalPendulum()
{
	Real	l = 1.0;
	SphericalPendulumODE		odeSysSpherPend = SphericalPendulumODE(l);

	Real	t1 = 0.0, t2 = 5.0;
	int   expectNumSteps = 100;
	Real	minSaveInterval = (t2 - t1) / expectNumSteps;
	Real	initAngleTheta = 0.5;
	Real  initAnglePhi = 0.1;
	Vector<Real>	initCond{ initAngleTheta, 0.0, initAnglePhi, 0.05 };

	ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(odeSysSpherPend);
	ODESystemSolution solAdapt = adaptSolver.integrate(initCond, t1, t2, minSaveInterval / 1.21, 1e-06, 0.01);

	Vector<Real> t_vals  = solAdapt.getTValues();
	Vector<Real> theta_vals = solAdapt.getXValues(0);
	Vector<Real> theta_dot_vals = solAdapt.getXValues(1);
	Vector<Real> phi_vals = solAdapt.getXValues(2);
	Vector<Real> phi_dot_vals = solAdapt.getXValues(3);

	std::cout << "\n\n****  Runge-Kutta 4th order - fixed stepsize  **********  Runge-Kutta 4th order - adaptive stepper  ****\n";
	std::vector<ColDesc>				vecNames{ ColDesc("t", 11, 2, 'F'), ColDesc("theta", 15, 8, 'F'), ColDesc("theta dot", 15, 8, 'F'),
																				ColDesc("phi", 12, 8, 'F'), ColDesc("phi dot", 12, 8, 'F') };
	std::vector<Vector<Real>*>	vecVals{ &t_vals, &theta_vals, &theta_dot_vals, &phi_vals, &phi_dot_vals  };
	VerticalVectorPrinter	vvp(vecNames, vecVals);

	vvp.Print();

	// getting solutions as polynomials
	PolynomInterpRealFunc 	solAdaptPolyInterp0 = solAdapt.getSolutionAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solAdaptPolyInterp1 = solAdapt.getSolutionAsPolyInterp(2, 3);

	SplineInterpRealFunc		solSplineInterp0 = solAdapt.getSolutionAsSplineInterp(0);
	SplineInterpParametricCurve<2> spline(solAdapt.getXValues());

	Visualizer::VisualizeRealFunction(solAdaptPolyInterp0, "Spherical pendulum - theta in time",
		t1, t2, 200, "spherical_pendulum_theta.txt");

	// shown together
	Visualizer::VisualizeMultiRealFunction({ &solAdaptPolyInterp0, &solAdaptPolyInterp1 },
		"Spherical pendulum - both angles", t1, t2, 200,
		"spherical_pendulum_multi_real_func.txt");

}