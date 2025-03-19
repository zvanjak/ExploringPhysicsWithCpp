#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Geometry3D.h"
#include "base/Function.h"
#include "base/InterpolatedFunction.h"

#include "core/Derivation.h"

#include "tools/Serializer.h"
#include "tools/Visualizer.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

//#include "../test_data/Fields.h"


using namespace MML;

// Simulating gravity for N bodies in 3D space
//   - multiple ways - simple, solving ODE
// Calculating potential lineas

class GravityMass
{
public:
	Real	_mass;
	Vector3Cartesian _position;
	Vector3Cartesian _velocity;

	GravityMass(const Real& mass, const Vector3Cartesian& position)
		: _mass(mass), _position(position) {
	}
	GravityMass(const Real& mass, const Vector3Cartesian& position, const Vector3Cartesian& velocity)
		: _mass(mass), _position(position), _velocity(velocity) {
	}
};

class NBodyGravityConfig
{
	std::vector<GravityMass> _masses;
public:
	int NumBodies() const { return (int)_masses.size(); }

	void AddBody(Real mass, Vector3Cartesian position, Vector3Cartesian velocity)
	{
		_masses.push_back(GravityMass(mass, position, velocity));
	}

	Real						 Mass(int i)		 const { return _masses[i]._mass; }
	Vector3Cartesian Position(int i) const { return _masses[i]._position; }
	Vector3Cartesian Velocity(int i) const { return _masses[i]._velocity; }

	void SetPosition(int i, Vector3Cartesian pos) { _masses[i]._position = pos; }
	void SetVelocity(int i, Vector3Cartesian vel) { _masses[i]._velocity = vel; }
};

class NBodyGravityPotentialField : public IScalarFunction<3>
{
private:
	Real _G;
	NBodyGravityConfig _config;
public:
	NBodyGravityPotentialField(NBodyGravityConfig config) : _G(6.67430e-11), _config(config) {}
	NBodyGravityPotentialField(Real G, NBodyGravityConfig config) : _G(G), _config(config) {}

	Real operator()(const VectorN<Real, 3>& x) const
	{
		Real pot = 0.0;
		for (int i = 0; i < _config.NumBodies(); i++)
		{
			// FIX!
			// pot += Fields::InverseRadialPotentialFieldCart(_G * _config.Mass(i), _config.Position(i) - x);
		}
		return pot;
	}
};
class NBodyGravityForceField : public IVectorFunction<3>
{
private:
	Real _G;
	NBodyGravityConfig _config;
public:
	NBodyGravityForceField(NBodyGravityConfig config) : _G(6.67430e-11), _config(config) {}
	NBodyGravityForceField(Real G, NBodyGravityConfig config) : _G(G), _config(config) {}

	VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const
	{
		VectorN<Real, 3> force(0.0);
		for (int i = 0; i < _config.NumBodies(); i++)
		{
			// FIX!
			// force = force + Fields::InverseRadialPotentialForceFieldCart(_G * _config.Mass(i), _config.Position(i) - x);
		}
		return force;
	}
};

// set up ODESystem za simuliranje sustava tijela u gravitacijskom polju
// a ima i druga opcija - u DANOM sustavu tijela, simulirati neko trece tijelo!
class NBodyMotionSolverEuler
{
protected:
	NBodyGravityConfig _config;

public:
	NBodyMotionSolverEuler(NBodyGravityConfig inConfig) : _config(inConfig) {}

	void simulateOneStep(const Real dt)
	{
		std::vector<Vector3Cartesian> force(_config.NumBodies(), Vector3Cartesian(0.0, 0.0, 0.0));

		for (int i = 0; i < _config.NumBodies(); i++) {
			// calculate force on body i
			for (int j = 0; j < _config.NumBodies(); j++) {
				if (i != j) {
					Real gravityConstant = 30;
					Vec3Cart radialVec = _config.Position(i) - _config.Position(j);
					force[i] = force[i] - gravityConstant * _config.Mass(i) * _config.Mass(j) / POW3(radialVec.NormL2()) * radialVec;
				}
			}
		}
		// advancing velocities and positions,  using the most simple Euler method of first order
		for (int i = 0; i < _config.NumBodies(); i++) {
			_config.SetVelocity(i, _config.Velocity(i) + force[i] * dt / _config.Mass(i));
			_config.SetPosition(i, _config.Position(i) + _config.Velocity(i) * dt);
		}
	}
	std::vector<std::vector<VectorN< Real, 3 >>> simulate(const Real dt, const int steps)
	{
		std::vector<std::vector<VectorN<Real, 3>>> trajectories(_config.NumBodies());

		// saving initial positions
		for (int i = 0; i < _config.NumBodies(); i++)
			trajectories[i].push_back(_config.Position(i));

		for (int i = 0; i < steps; i++) {
			simulateOneStep(dt);

			for (int i = 0; i < _config.NumBodies(); i++)
				trajectories[i].push_back(_config.Position(i));
		}
		return trajectories;
	}
};

void Example9_N_body_problem_Euler()
{
	NBodyGravityConfig config;
	config.AddBody(10000, Vector3Cartesian{ 0.0,    0.0,    0.0 }, Vector3Cartesian{ 0.0,   2.0,   5.0 });
	config.AddBody(20, Vector3Cartesian{ -110.0,  -50.0,   10.0 }, Vector3Cartesian{ 0.0,   50,   0.0 });
	config.AddBody(10, Vector3Cartesian{ 130.0,   50.0,   70.0 }, Vector3Cartesian{ 0.0,   -50,   0 });
	config.AddBody(20, Vector3Cartesian{ -20.0,  100.0, -110.0 }, Vector3Cartesian{ 50,   0.0,   0.0 });
	config.AddBody(10, Vector3Cartesian{ 70.0, -110.0,   70.0 }, Vector3Cartesian{ -50,   50,   50.0 });

	NBodyMotionSolverEuler solver(config);

	Real t1 = 0.0, t2 = 20.0;
	const int  steps = 2000;
	const Real dt = (t2 - t1) / steps;

	auto res = solver.simulate(dt, steps);

	for (int i = 0; i < config.NumBodies(); i++)
	{
		Serializer::SaveAsParamCurve<3>(res[i], "PARAMETRIC_CURVE_CARTESIAN_3D", "Body" + std::to_string(i + 1),
			0.0, dt * steps, steps + 1,
			MML_PATH_ResultFiles + "body" + std::to_string(i) + ".txt");
	}

	Visualizer::VisualizeMultiParamCurve3D({ "body0.txt", "body1.txt", "body2.txt", "body3.txt", "body4.txt" });
}


class NBodySystemODE : public IODESystem
{
	NBodyGravityConfig _config;
	const Real G = 30; //6.67430e-11;

public:
	NBodySystemODE(NBodyGravityConfig inConfig) : _config(inConfig) {}

	int getDim() const { return 6 * _config.NumBodies(); }
	void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const
	{
		// filling in dxdt vector with N * 6 derivations of our variables (x, y, z, vx, vy, vz)
		for (int i = 0; i < _config.NumBodies(); i++)
		{
			// calculating force on body i
			Vector3Cartesian force(0, 0, 0);
			for (int j = 0; j < _config.NumBodies(); j++)
			{
				if (i != j)      // checking for self-force 
				{
					Vector3Cartesian vec_dist(x[6 * j] - x[6 * i],
																		x[6 * j + 2] - x[6 * i + 2],
																		x[6 * j + 4] - x[6 * i + 4]);

					force = force + G * _config.Mass(i) * _config.Mass(j) / POW3(vec_dist.NormL2()) * vec_dist;
				}
			}

			// x coord
			dxdt[6 * i]			= x[6 * i + 1];
			dxdt[6 * i + 1] = 1 / _config.Mass(i) * force.X();
			// y coord
			dxdt[6 * i + 2] = x[6 * i + 3];
			dxdt[6 * i + 3] = 1 / _config.Mass(i) * force.Y();
			// z coord
			dxdt[6 * i + 4] = x[6 * i + 5];
			dxdt[6 * i + 5] = 1 / _config.Mass(i) * force.Z();
		}
	}
};

void Example9_N_body_problem_RK5Adaptive()
{
	NBodyGravityConfig config;
	config.AddBody(10000, Vector3Cartesian{ 0.0,    0.0,    0.0 }, Vector3Cartesian{ 0.0,   2.0,   5.0 });
	config.AddBody(20, Vector3Cartesian{ -110.0,  -50.0,   10.0 }, Vector3Cartesian{ 0.0,   50,   0.0 });
	config.AddBody(10, Vector3Cartesian{ 130.0,   50.0,   70.0 }, Vector3Cartesian{ 0.0,   -50,   0 });
	config.AddBody(20, Vector3Cartesian{ -20.0,  100.0, -110.0 }, Vector3Cartesian{ 50,   0.0,   0.0 });
	config.AddBody(10, Vector3Cartesian{ 70.0, -110.0,   70.0 }, Vector3Cartesian{ -50,   50,   50.0 });

	NBodySystemODE odeSys(config);

	// set-up initial conditions, according to config
	Vector<Real>	initCond(6 * config.NumBodies());
	for(int i=0; i<config.NumBodies(); i++ )
	{
		initCond[6 * i] = config.Position(i).X();
		initCond[6 * i + 1] = config.Velocity(i).X();
		initCond[6 * i + 2] = config.Position(i).Y();
		initCond[6 * i + 3] = config.Velocity(i).Y();
		initCond[6 * i + 4] = config.Position(i).Z();
		initCond[6 * i + 5] = config.Velocity(i).Z();
	}

	Real	t1 = 0.0, t2 = 20.0;
	int   expectNumSteps = 1000;
	Real	minSaveInterval = (t2 - t1) / expectNumSteps;

	ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(odeSys);
	ODESystemSolution solAdapt = adaptSolver.integrate(initCond, t1, t2, minSaveInterval, 1e-07, 0.05, 0.0);

	std::cout << "Done steps: OK - " << solAdapt.getNumStepsOK() << "  Bad - " << solAdapt.getNumStepsBad() << std::endl;

	Vector<Real> t_vals = solAdapt.getTValues();
	for( int i=0; i<config.NumBodies(); i++ )
	{
		Vector<Real> body_x_vals = solAdapt.getXValues(6 * i);
		Vector<Real> body_y_vals = solAdapt.getXValues(6 * i + 2);
		Vector<Real> body_z_vals = solAdapt.getXValues(6 * i + 4);

		// form ParametricCurve from these 3 vectors
		std::vector<VectorN<Real, 3>> res;
		for (int i = 0; i < t_vals.size(); i++)
			res.push_back(VectorN<Real, 3>{body_x_vals[i], body_y_vals[i], body_z_vals[i]});

		Serializer::SaveAsParamCurve<3>(res, "PARAMETRIC_CURVE_CARTESIAN_3D", "body_RK4_" + std::to_string(i + 1),
																		t1, t2, t_vals.size(),
																		MML_PATH_ResultFiles + "body_rk4_" + std::to_string(i) + ".txt");
	}

	Visualizer::VisualizeMultiParamCurve3D({ "body_rk4_0.txt", "body_rk4_1.txt", "body_rk4_2.txt", "body_rk4_3.txt", "body_rk4_4.txt" });
}

int main()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****               EXAMPLE 9 - solving N-body problem              ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	//Example9_N_body_problem_Euler();
	Example9_N_body_problem_RK5Adaptive();

	return 0;
}

