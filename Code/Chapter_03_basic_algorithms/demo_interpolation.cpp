#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "base/Function.h"
#include "base/InterpolatedFunction.h"

#include "tools/Visualizer.h"

using namespace MML;


void Demo_Interpolation()
{
	// we are manually creating data for interpolation in two vectors
	// but generally you would have some data from some source
	Real x1 = 0.0;
	Real x2 = 10.0;
	Vector<Real> vec_x{ 0.0, 1.0, 2.0,  3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
	Vector<Real> vec_y{ 0.0, 1.0, 4.0, -2.0, -16, -13,   6, 6.2, 5.5, 12.0, 3.0 };

	// create InterpolatedFunction object
	LinearInterpRealFunc    lin_interp(vec_x, vec_y);
	PolynomInterpRealFunc   poly_interp3(vec_x, vec_y, 3);
	PolynomInterpRealFunc   poly_interp4(vec_x, vec_y, 4);
	SplineInterpRealFunc		spline_interp(vec_x, vec_y);

	// visualize interpolations
  Visualizer::VisualizeRealFunction(lin_interp, "linear interpolation", 
																		x1, x2, 100, "example2_lin_interp.txt");
	Visualizer::VisualizeRealFunction(poly_interp3, "polynomial interpolation 3rd order", 
																		x1, x2, 100, "example2_poly_interp3.txt");
	Visualizer::VisualizeRealFunction(poly_interp4, "polynomial interpolation 4th order", 
																		x1, x2, 100, "example2_poly_interp4.txt");
	Visualizer::VisualizeRealFunction(spline_interp, "spline interpolation", 
																		x1, x2, 100, "example2_spline_interp.txt");

	// visualize all together
	Visualizer::VisualizeMultiRealFunction({ &lin_interp, &poly_interp3, &poly_interp4, &spline_interp },
																					"Interpolation comparison", { "Linear", "Poly 3", "Poly 4", "Spline" },
																					0.0, 10.0, 100, "example2_interpolation.txt");
}