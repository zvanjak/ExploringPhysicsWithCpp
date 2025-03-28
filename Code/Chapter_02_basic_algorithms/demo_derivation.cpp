#include "MMLBase.h"

#include "base/Function.h"
#include "core/Derivation.h"

using namespace MML;


void Example2_RealFunction_Derivation()
{
	RealFunction       f1{ [](Real x) { return (Real)(sin(x) * (1.0 + 0.5 * x * x)); } };

	// numerical derivation of real function (available orders - 1, 2, 4, 6, 8)
	double der_f1 = Derivation::NDer1(f1, 0.5);
	double der_f2 = Derivation::NDer2(f1, 0.5, 1e-6);   // setting explicit step size
	Real err;
	double der_f6 = Derivation::NDer6(f1, 0.5, &err);   // if we need error estimate    
	double der_f8 = Derivation::NDer8(f1, 0.5);
	
	// we can use default Derive routine (set to NDer4)
	double der_f4 = Derivation::Derive(f1, 0.5);
	// available also with error estimate
	double der_f41 = Derivation::DeriveErr(f1, 0.5, &err);

	// second and third derivatives
	double sec_der_f1		= Derivation::NSecDer2(f1, 0.5);
	double third_der_f1 = Derivation::NThirdDer2(f1, 0.5);

	// creating new function that is derivation of existing function
	//RealFuncDerived4    f1_der4(f1);        // 4th order derivation
}

void Example2_ScalarVectorFunction_Derivation()
{
	// scalar and vector functions
	ScalarFunction<3>   f2Scal([](const VectorN<Real, 3>& x) { return (Real)(1.0 / pow(x.NormL2(), 2)); });
	VectorFunction<3>   f3Vec([](const VectorN<Real, 3>& x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });

	VectorN<Real, 3> der_point{ 1.0, 1.0, 1.0 };

	double					 der_f2Scal			= Derivation::NDer1Partial(f2Scal, 1, der_point);
	VectorN<Real, 3> der_f2Scal_all = Derivation::NDer1PartialByAll(f2Scal, der_point);

	double							 der_f3Vec				= Derivation::NDer1Partial(f3Vec, 1, 1, der_point);
	VectorN<Real, 3>     der_f3Vec_by1		= Derivation::NDer2PartialByAll(f3Vec, 1, der_point);
	MatrixNM<Real, 3, 3> der_f3Vec_by_all = Derivation::NDer4PartialAllByAll(f3Vec, der_point);
}

void Example2_ParametricCurve_Derivation()
{
	// parametric curves
	ParametricCurve<3>  f4Curve([](Real x) { return VectorN<Real, 3>{x, 2 * x, 3 * x}; });

	VectorN<Real, 3>    der_f4Curve   = Derivation::NDer1(f4Curve, 1.0);
	VectorN<Real, 3>    sec_f4Curve   = Derivation::NSecDer1(f4Curve, 1.0);
	VectorN<Real, 3>    third_f4Curve = Derivation::NThirdDer1(f4Curve, 1.0);
}

void Demo_Derivation()
{
	Example2_RealFunction_Derivation();
	Example2_ScalarVectorFunction_Derivation();
	Example2_ParametricCurve_Derivation();
}