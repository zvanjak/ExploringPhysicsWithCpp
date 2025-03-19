#include "MMLBase.h"
#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "core/Derivation.h"
#include "base/Function.h"

#include "base/VectorTypes.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "tools/Serializer.h"
#include "tools/Visualizer.h"


using namespace MML;


class TwoBodyGravityConfig
{
public:
	Real _G = 100;
	Real _m1, _m2;
	Vector3Cartesian _r1, _r2;
	Vector3Cartesian _v1, _v2;

	TwoBodyGravityConfig(Real m1, Real m2, Vector3Cartesian r1, Vector3Cartesian r2, Vector3Cartesian v1, Vector3Cartesian v2)
		: _m1(m1), _m2(m2), _r1(r1), _r2(r2), _v1(v1), _v2(v2)
	{	}

	TwoBodyGravityConfig(Real G, Real m1, Real m2, Vector3Cartesian r1, Vector3Cartesian r2, Vector3Cartesian v1, Vector3Cartesian v2)
		: _G(G), _m1(m1), _m2(m2), _r1(r1), _r2(r2), _v1(v1), _v2(v2)
	{	}

};

class TwoBodyGravitySystemODE : public IODESystem
{
	TwoBodyGravityConfig _config;

public:
	TwoBodyGravitySystemODE(const TwoBodyGravityConfig& config)
		: _config(config)
	{	}

	int getDim() const override { return 12; }
	void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
	{
		Vector3Cartesian r1{ x[0], x[1], x[2] };
		Vector3Cartesian r2{ x[3], x[4], x[5] };
		Vector3Cartesian v1{ x[6], x[7], x[8] };
		Vector3Cartesian v2{ x[9], x[10], x[11] };

		Vector3Cartesian r12 = r2 - r1;
		Real dist = r12.NormL2();
		
		Vector3Cartesian f12 = _config._G * _config._m1 * _config._m2 / POW3(dist) * r12;
		Vector3Cartesian f21 = -f12;
		
		dxdt[0] = v1[0];
		dxdt[1] = v1[1];
		dxdt[2] = v1[2];
		dxdt[3] = v2[0];
		dxdt[4] = v2[1];
		dxdt[5] = v2[2];
		dxdt[6] = f12[0] / _config._m1;
		dxdt[7] = f12[1] / _config._m1;
		dxdt[8] = f12[2] / _config._m1;
		dxdt[9] = f21[0] / _config._m2;
		dxdt[10] = f21[1] / _config._m2;
		dxdt[11] = f21[2] / _config._m2;
	}
};

class TwoBodiesGravitySimulator
{

public:

	// solve and visualize trajectory
	void SolveAndShowTrajectories(const TwoBodyGravityConfig &config, Real t)
	{
		TwoBodyGravitySystemODE ode(config);
		
		ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(ode);
		ODESystemSolution solAdapt = adaptSolver.integrate(Vector<Real>{ config._r1[0], config._r1[1], config._r1[2], 
																																		 config._r2[0], config._r2[1], config._r2[2], 
																																		 config._v1[0], config._v1[1], config._v1[2], 
																																		 config._v2[0], config._v2[1], config._v2[2] }, 
																											 0, t, 0.01, 1e-06, 0.01);
		
		Vector<Real> t_vals = solAdapt.getTValues();
		Vector<Real> r1_x_vals = solAdapt.getXValues(0);
		Vector<Real> r1_y_vals = solAdapt.getXValues(1);
		Vector<Real> r1_z_vals = solAdapt.getXValues(2);
		Vector<Real> r2_x_vals = solAdapt.getXValues(3);
		Vector<Real> r2_y_vals = solAdapt.getXValues(4);
		Vector<Real> r2_z_vals = solAdapt.getXValues(5);
		
		// form ParametricCurve from these 3 vectors
		std::vector<VectorN<Real, 3>> res1;
		for (int i = 0; i < t_vals.size(); i++)
			res1.push_back(VectorN<Real, 3>{r1_x_vals[i], r1_y_vals[i], r1_z_vals[i]});
		Serializer::SaveAsParamCurve<3>(res1, "PARAMETRIC_CURVE_CARTESIAN_3D", "gravity_example_body1",
																		0, t, t_vals.size(),
																		MML_PATH_ResultFiles + "gravity_example_body1.txt");
		
		std::vector<VectorN<Real, 3>> res2;
		for (int i = 0; i < t_vals.size(); i++)
			res2.push_back(VectorN<Real, 3>{r2_x_vals[i], r2_y_vals[i], r2_z_vals[i]});
		Serializer::SaveAsParamCurve<3>(res2, "PARAMETRIC_CURVE_CARTESIAN_3D", "gravity_example_body2",
																		0, t, t_vals.size(),
																		MML_PATH_ResultFiles + "gravity_example_body2.txt");

		Visualizer::VisualizeMultiParamCurve3D({ "gravity_example_body1.txt", "gravity_example_body2.txt" });
	}

	// calculate center of mass
	Vec3Cart CenterOfMass(const TwoBodyGravityConfig &config) const
	{
		return (config._m1 * config._r1 + config._m2 * config._r2) / (config._m1 + config._m2);
	}

	// calculate reduced mass
	Real ReducedMass(const TwoBodyGravityConfig &config) const
	{
		return config._m1 * config._m2 / (config._m1 + config._m2);
	}

	// calculate total potential energy
	Real TotalPotentialEnergy(const TwoBodyGravityConfig &config) const
	{
		return -config._G * config._m1 * config._m2 / (config._r1 - config._r2).NormL2();
	}

	// calculate total kinetic energy
	Real TotalKineticEnergy(const TwoBodyGravityConfig &config) const
	{
		return 0.5 * config._m1 * config._v1.NormL2() + 0.5 * config._m2 * config._v2.NormL2();
	}
};

// SolarSystem gravity simulator

class SolarSystemGravitySimulator
{
	// ima definiranu putanju za sve planete
	ParametricCurve<3> _planetPaths[9];
};

void Demo_TwoMasses()
{
	Real G = 200;
	Real m1 = 100;
	Real m2 = 200;
	Vec3Cart r1{ -100, -50, -50 };
	Vec3Cart r2{ 50, -30, -40 };
	Vec3Cart v1{ -5, 0, 5.0 };
	Vec3Cart v2{ 5, 2, -5.0 };

	TwoBodyGravityConfig config(G, m1, m2, r1, r2, v1, v2);

	TwoBodiesGravitySimulator sim;

	std::cout << "Center of mass   = " << sim.CenterOfMass(config) << std::endl;
	std::cout << "Kinetic energy   = " << sim.TotalKineticEnergy(config) << std::endl;
	std::cout << "Potential energy = " << sim.TotalPotentialEnergy(config) << std::endl;

	sim.SolveAndShowTrajectories(config, 1000);
}

int main()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                       EXAMPLE 8 - Gravity                     ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Demo_TwoMasses();

	// to je small body simulator, u zadanom polju tijela, koje se gibaju po zadanim putanjama

	// imamo zadane putanje Sunca i svih planeta kao ellipse curve

	// imammo, u prikladnim heliocentricnim koordinatam, i putanju Voyagera 1 i 2 kroz vrijeme

	// CILJ je, na osnovu dane putanje, odsmilurati precizno, kad je Voyager morao paliti motore
	// odnosno kad je bilo dodatne silen na njega, da bi se zadrzao na putanji

	return 0;
}

