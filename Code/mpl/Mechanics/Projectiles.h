#if !defined MPL_PROJECTILES_H
#define MPL_PROJECTILES_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

using namespace MML;

namespace MPL
{
	class IProjectileODE : public IODESystem
	{
	public:
		virtual Vector<Real> getInitCond(Real angle, Real initHeight, Real velocity) const = 0;
	};

	// representing ODE system for projectile motion in vacuum
	class ProjectileInVacuumODE : public IProjectileODE
	{
	public:
		int  getDim() const override { return 4; }
		void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[2];
			dxdt[1] = x[3];
			dxdt[2] = 0;
			dxdt[3] = -9.81;

			if (x[1] < 0) // if the projectile hits the ground, stop the motion
			{
				dxdt[0] = 0; // stop x motion
				dxdt[1] = 0; // stop y motion
				dxdt[2] = 0; // stop x velocity
				dxdt[3] = 0; // stop y velocity

				//std::cout << "Projectile hit the ground at t = " << t << std::endl;
			}
		}
		Vector<Real> getInitCond(Real angle, Real initHeight, Real velocity) const
		{
			Vector<Real> initCond(4);
			initCond[0] = 0;											// initial x position
			initCond[1] = initHeight;							// initial y position
			initCond[2] = velocity * cos(angle);	// initial x velocity (m/s)
			initCond[3] = velocity * sin(angle);	// initial y velocity (m/s)
			return initCond;
		}

		Real CalcRange(Real angle, Real initHeight, Real velocity) const
		{
			// calculating range of the projectile motion in vacuum
			Real g = 9.81; // gravitational acceleration
			return (velocity * cos(angle) / g) * (velocity * sin(angle) + sqrt(velocity * velocity * sin(angle) * sin(angle) + 2 * g * initHeight));
		}
		Real TimeOfFlight(Real angle, Real initHeight, Real velocity) const
		{
			// calculating time of flight for the projectile motion in vacuum
			Real g = 9.81;
			Real v0y = velocity * sin(angle);
			return (v0y + sqrt(v0y * v0y + 2 * g * initHeight)) / g;
		}
	};

	class ProjectileWithAirResistanceODE : public IProjectileODE
	{
		Real _dragCoefficient; // air resistance coefficient
	public:
		ProjectileWithAirResistanceODE(Real dragCoefficient = 0.1) : _dragCoefficient(dragCoefficient) {}

		int  getDim() const override { return 4; }
		void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
		{
			Real v = sqrt(x[2] * x[2] + x[3] * x[3]);			// speed of the projectile

			dxdt[0] = x[2];
			dxdt[1] = x[3];
			dxdt[2] = -_dragCoefficient * v * x[2];
			dxdt[3] = -9.81 - _dragCoefficient * v * x[3];

			if (x[1] < 0) // if the projectile hits the ground, stop the motion
			{
				dxdt[0] = 0; // stop x motion
				dxdt[1] = 0; // stop y motion
				dxdt[2] = 0; // stop x velocity
				dxdt[3] = 0; // stop y velocity
				//std::cout << "Projectile hit the ground at t = " << t << std::endl;
			}
		}
		Vector<Real> getInitCond(Real angle, Real initHeight, Real velocity) const
		{
			Vector<Real> initCond(4);
			initCond[0] = 0;											// initial x position
			initCond[1] = initHeight;							// initial y position
			initCond[2] = velocity * cos(angle);	// initial x velocity (m/s)
			initCond[3] = velocity * sin(angle);	// initial y velocity (m/s)
			return initCond;
		}
	};

	// representing ODE system for projectile motion with air resistance, 
	// but with included changes of air density with altitude
	enum class AirDensityModel
	{
		Isothermal,		// isothermal model of air density
		Adiabatic			// adiabatic model of air density
	};

	class ProjectileChangingAirDensityODE : public IProjectileODE
	{
		Real _dragCoefficient0;					// air resistance coefficient at sea level
		Real _airDensity0 = 1.225;			// air density at sea level in kg/m^3
		
		AirDensityModel _airDensityModel = AirDensityModel::Isothermal; // model of air density

	public:
		ProjectileChangingAirDensityODE(Real dragCoefficient0 = 0.1) : _dragCoefficient0(dragCoefficient0) {}

		int  getDim() const override { return 4; }
		void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
		{
			Real v = sqrt(x[2] * x[2] + x[3] * x[3]);				// speed of the projectile

			dxdt[0] = x[2];
			dxdt[1] = x[3];
			dxdt[2] = -dragCoefficient(x[3], v) * v * x[2];
			dxdt[3] = -9.81 - dragCoefficient(x[3], v) * v * x[3];
		}
		Vector<Real> getInitCond(Real angle, Real initHeight, Real velocity) const
		{
			Vector<Real> initCond(4);
			initCond[0] = 0;											// initial x position
			initCond[1] = initHeight;							// initial y position
			initCond[2] = velocity * cos(angle);	// initial x velocity (m/s)
			initCond[3] = velocity * sin(angle);	// initial y velocity (m/s)
			return initCond;
		}

		// function to calculate air density based on altitude (isothermal)
		Real airDensity(Real altitude) const
		{
			// Simplified model: air density decreases exponentially with altitude
			// at sea level, air density is approximately 1.225 kg/m^3
			return 1.225 * exp(-altitude / 8000); // 8000 m is a rough scale height for the atmosphere
		}

		// function to calculate air density based on altitude (adiabatic)
		Real airDensityAdiabatic(Real altitude) const
		{
			// Adiabatic model: air density decreases with altitude, but with a different rate
			// at sea level, air density is approximately 1.225 kg/m^3
			return 1.225 * pow((1 - (0.0065 * altitude / 288.15)), 5.2561); // using the adiabatic lapse rate
		}

		// function to calculate drag coefficient based on air density and speed
		Real dragCoefficient(Real altitude, Real speed) const
		{
			if( _airDensityModel == AirDensityModel::Adiabatic )
			{
				Real density = airDensityAdiabatic(altitude);
				return _dragCoefficient0 * density / _airDensity0; // simplified drag coefficient formula
			}
			
			// Default to isothermal model
			Real density = airDensity(altitude);
			return _dragCoefficient0 * density / _airDensity0; // simplified drag coefficient formula
		}
	};

	// representing ODE system for projectile motion with air resistance,
	// but with drag coefficient dependent on speed
	class BaseballWithDragCoeffDependentOnSpeedODE : public IProjectileODE
	{
	public:
		BaseballWithDragCoeffDependentOnSpeedODE()	{ }

		int  getDim() const override { return 4; }
		void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
		{
			Real v = sqrt(x[2] * x[2] + x[3] * x[3]);				// speed of the projectile

			dxdt[0] = x[2];
			dxdt[1] = x[3];
			dxdt[2] = -dragCoefficient(v) * v * x[2];
			dxdt[3] = -9.81 - dragCoefficient(v) * v * x[3];
		}
		Vector<Real> getInitCond(Real angle, Real initHeight, Real velocity) const
		{
			Vector<Real> initCond(4);
			initCond[0] = 0;											// initial x position
			initCond[1] = initHeight;							// initial y position
			initCond[2] = velocity * cos(angle);	// initial x velocity (m/s)
			initCond[3] = velocity * sin(angle);	// initial y velocity (m/s)
			return initCond;
		}

		Real dragCoefficient(Real v) const
		{
			Real vd = 35;
			Real delta = 5;

			return 0.0039 + 0.0058 / (1 + exp((v - vd) / delta));
		}
	};
}

#endif