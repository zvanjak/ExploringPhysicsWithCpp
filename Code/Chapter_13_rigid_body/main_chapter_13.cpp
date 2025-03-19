#include "MMLBase.h"

#include "base/VectorN.h"

#include "core/Derivation.h"
#include "core/Integration.h"
#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransf3D.h"


using namespace MML;


struct DiscreteMass {
    Vector3Cartesian _position;
    double _mass;

    DiscreteMass(const Vector3Cartesian &position, const double& mass)
        : _position(position), _mass(mass)  { }
};
struct DiscreteMassesConfig {
    std::vector<DiscreteMass> _masses;

    DiscreteMassesConfig(const std::vector<DiscreteMass>& masses)
        : _masses(masses)  { }
};

class DiscreteMassMomentOfInertiaTensorCalculator
{
    DiscreteMassesConfig _massesConfig;
public:
    DiscreteMassMomentOfInertiaTensorCalculator(const DiscreteMassesConfig& massesConfig)
        : _massesConfig(massesConfig) { }

    Tensor2<3> calculate()
    {
        Tensor2<3> tensor(2,0);  // can be (0,2) or (1,1) as well (it is a Cartesian tensor)
        for (const auto& mass : _massesConfig._masses)
        {
            Vector3Cartesian pos = mass._position;
            tensor(0,0) += mass._mass * (pos.Y() * pos.Y() + pos.Z() * pos.Z());
            tensor(1,1) += mass._mass * (pos.X() * pos.X() + pos.Z() * pos.Z());
            tensor(2,2) += mass._mass * (pos.X() * pos.X() + pos.Y() * pos.Y());
            
            tensor(0,1) -= mass._mass * pos.X() * pos.Y();
            tensor(0,2) -= mass._mass * pos.X() * pos.Z();
            tensor(1,2) -= mass._mass * pos.Y() * pos.Z();
        }
        tensor(1,0) = tensor(0,1);
        tensor(2,0) = tensor(0,2);
        tensor(2,1) = tensor(1,2);

        return tensor;
    }
};

// continuous mass
class IContinuousMass
{
public:
    Real _x1, _x2;

    Real (*_y1)(Real);
    Real (*_y2)(Real);

    Real (*_z1)(Real, Real);
    Real (*_z2)(Real, Real);

    IContinuousMass(Real x1, Real x2, Real (*y1)(Real), Real (*y2)(Real), Real (*z1)(Real, Real), Real (*z2)(Real, Real))
        : _x1(x1), _x2(x2), _y1(y1), _y2(y2), _z1(z1), _z2(z2)
    { }  

    virtual Real getDensity(const VectorN<Real, 3> &x) = 0;
};

class ContinuousMass : public IContinuousMass
{
    Real (*_density)(const VectorN<Real, 3> &x);

public:
    ContinuousMass(Real x1, Real x2, Real (*y1)(Real), Real (*y2)(Real), Real (*z1)(Real, Real), Real (*z2)(Real, Real), Real (*density)(const VectorN<Real, 3> &x))
        : IContinuousMass(x1, x2, y1, y2, z1, z2), _density(density)
    { } 

    virtual Real getDensity(const VectorN<Real, 3> &x)
    {
        return _density(x);
    }
};

class ContinuousMassConstDensity : public IContinuousMass
{
    Real _density;
public:
    ContinuousMassConstDensity(Real x1, Real x2, Real (*y1)(Real), Real (*y2)(Real), Real (*z1)(Real, Real), Real (*z2)(Real, Real), Real density)
        : IContinuousMass(x1, x2, y1, y2, z1, z2), _density(density)
    { } 

    virtual Real getDensity(const VectorN<Real, 3> &x)
    {
        return _density;
    }
};

class ContMassScalarFuncBase : public IScalarFunction<3>
{
protected:
  IContinuousMass &_mass;
public:
	ContMassScalarFuncBase(IContinuousMass& mass) : _mass(mass) {}
};

struct Func11 : public ContMassScalarFuncBase
{
	Func11(IContinuousMass& mass) : ContMassScalarFuncBase(mass) {}
	Real operator()(const VectorN<Real, 3>& x) const
	{
		return _mass.getDensity(x) * (x[1]*x[1] + x[2]*x[2]);
	}
};
struct Func22 : public ContMassScalarFuncBase
{
	Func22(IContinuousMass& mass) : ContMassScalarFuncBase(mass) {}
	Real operator()(const VectorN<Real, 3>& x) const
	{
		return _mass.getDensity(x) * (x[0] * x[0] + x[2] * x[2]);
	}
};
struct Func33 : public ContMassScalarFuncBase
{
	Func33(IContinuousMass& mass) : ContMassScalarFuncBase(mass) {}
	Real operator()(const VectorN<Real, 3>& x) const
	{
		return _mass.getDensity(x) * (x[0] * x[0] + x[1] * x[1]);
	}
};
struct Func12 : public ContMassScalarFuncBase
{
	Func12(IContinuousMass& mass) : ContMassScalarFuncBase(mass) {}
	Real operator()(const VectorN<Real, 3>& x) const
	{
		return -_mass.getDensity(x) * x[0] * x[1];
	}
};
struct Func13 : public ContMassScalarFuncBase
{
	Func13(IContinuousMass& mass) : ContMassScalarFuncBase(mass) {}
	Real operator()(const VectorN<Real, 3>& x) const
	{
		return -_mass.getDensity(x) * x[0] * x[2];
	}
};
struct Func23 : public ContMassScalarFuncBase
{
	Func23(IContinuousMass& mass) : ContMassScalarFuncBase(mass) {}
	Real operator()(const VectorN<Real, 3>& x) const
	{
		return -_mass.getDensity(x) * x[1] * x[2];
	}
};

class ContinuousMassMomentOfInertiaTensorCalculator
{
    IContinuousMass &_mass;
public:
    ContinuousMassMomentOfInertiaTensorCalculator(IContinuousMass &mass)
        : _mass(mass)
    { }

		// calculating tensor of inertia for continuous mass around z axis
    Tensor2<3> calculate()
    {
        Tensor2<3> tensor(2,0);

				Func11 _f11(_mass);
				Func22 _f22(_mass);
				Func33 _f33(_mass);
				Func12 _f12(_mass);
				Func13 _f13(_mass);
				Func23 _f23(_mass);
				tensor(0, 0) = Integrate3D(_f11, _mass._x1, _mass._x2, _mass._y1, _mass._y2, _mass._z1, _mass._z2);
				tensor(1, 1) = Integrate3D(_f22, _mass._x1, _mass._x2, _mass._y1, _mass._y2, _mass._z1, _mass._z2);
				tensor(2, 2) = Integrate3D(_f33, _mass._x1, _mass._x2, _mass._y1, _mass._y2, _mass._z1, _mass._z2);
				tensor(0, 1) = Integrate3D(_f12, _mass._x1, _mass._x2, _mass._y1, _mass._y2, _mass._z1, _mass._z2);
				tensor(0, 2) = Integrate3D(_f13, _mass._x1, _mass._x2, _mass._y1, _mass._y2, _mass._z1, _mass._z2);
				tensor(1, 2) = Integrate3D(_f23, _mass._x1, _mass._x2, _mass._y1, _mass._y2, _mass._z1, _mass._z2);
				tensor(1, 0) = tensor(0, 1);
				tensor(2, 0) = tensor(0, 2);
				tensor(2, 1) = tensor(1, 2);

        ScalarFunctionFromStdFunc<3> fDensity( std::function<Real(const VectorN<Real, 3>&)>{ std::bind( &IContinuousMass::getDensity, &_mass, std::placeholders::_1) } );

        Real vol = Integrate3D( fDensity, 
                                _mass._x1, _mass._x2,              
                                _mass._y1, _mass._y2,
                                _mass._z1, _mass._z2);
                                    
        return tensor;
    }
};

void Example13_tensor_of_inertia_discrete_masses()
{
    // define set of discrete masses
    double a = 1;
    Vector3Cartesian pos1(a, a, 0);
    Vector3Cartesian pos2(-a, a, 0);
    Vector3Cartesian pos3(-a, -a, 0);
    Vector3Cartesian pos4(a, -a, 0);
    Vector3Cartesian pos5(0, 0, 4*a);
    double m1 = 2;
    double m2 = 1;
    double m3 = 4;
    double m4 = 1;
    double m5 = 1;

    DiscreteMass mass1(pos1, m1);
    DiscreteMass mass2(pos2, m2);
    DiscreteMass mass3(pos3, m3);
    DiscreteMass mass4(pos4, m4);
    DiscreteMass mass5(pos5, m5);

    std::vector<DiscreteMass> masses = {mass1, mass2, mass3, mass4, mass5};
    DiscreteMassesConfig massesConfig(masses);

    DiscreteMassMomentOfInertiaTensorCalculator calculator(massesConfig);
    Tensor2<3> tensor_orig = calculator.calculate();

    std::cout << "Tensor of inertia: " << std::endl;
    std::cout << tensor_orig << std::endl;

    // investigating what happens if we change coord.system, in two cases:
    // 1. using coord.system transform we calculate TRANSFORMED (original) tensor
    // 2. using coord.system transform we calculate NEW set of masses and then calculate tensor

    // new coord.system is rotated around x axis for 30 degrees
    CoordTransfCart3DRotationXAxis coord_transf(30.0 * Constants::PI / 180.0);

    // 1) - calculated tensor transformation
    Tensor2<3> tensor_transf = coord_transf.transfTensor2(tensor_orig, Vector3Cartesian(1, 1, 1));

    std::cout << "Tensor of inertia transformed: " << std::endl;
    std::cout << tensor_transf << std::endl;

    // 2) - change masses position and calculate new tensor 
    DiscreteMassesConfig massesTransConfig(masses);
    for (auto& mass : massesTransConfig._masses)
		mass._position = coord_transf.transf(mass._position);

    DiscreteMassMomentOfInertiaTensorCalculator calculator2(massesTransConfig);
    Tensor2<3> tensor_changed = calculator2.calculate();

    std::cout << "Tensor of inertia rotated masses: " << std::endl;
    std::cout << tensor_changed << std::endl;
}

void Example13_tensor_of_intertia_continuous_mass()
{
  // create ContinuousMass representation for cube of constant density
	ContinuousMassConstDensity cube(-0.5, 0.5, 
                                  [](Real x) { return -0.5; }, [](Real x) { return 0.5; }, 
                                  [](Real x, Real y) { return -0.5; }, [](Real x, Real y) { return 0.5; }, 
                                  1.0);

	ContinuousMassMomentOfInertiaTensorCalculator calculator(cube);

	Tensor2<3> tensor = calculator.calculate();

	std::cout << "Tensor of inertia for cube: " << std::endl;
	std::cout << tensor << std::endl;

  // let's do the same for sphere
	ContinuousMassConstDensity sphere(-1.0, 1.0,
		[](Real x) { return -sqrt(1 - x * x); }, [](Real x) { return sqrt(1 - x * x); },
		[](Real x, Real y) { return -sqrt(1 - x * x - y * y); }, [](Real x, Real y) { return sqrt(1 - x * x - y * y); },
		1.0);

	ContinuousMassMomentOfInertiaTensorCalculator calculator2(sphere);

	Tensor2<3> tensor2 = calculator2.calculate();

	std::cout << "Tensor of inertia for sphere: " << std::endl;
	std::cout << tensor2 << std::endl;

	std::cout << "Theoretical value for sphere: " << std::endl;
	std::cout << 2.0 / 5.0 * (4.0 / 3.0 * Constants::PI)  << std::endl;
}

void Example13_tensor_of_inertia()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                   EXAMPLE 13 - tensor of inertia               ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Example13_tensor_of_inertia_discrete_masses();
    Example13_tensor_of_intertia_continuous_mass();
}

int main()
{
	Example13_tensor_of_inertia();
    
  return 0;
}

