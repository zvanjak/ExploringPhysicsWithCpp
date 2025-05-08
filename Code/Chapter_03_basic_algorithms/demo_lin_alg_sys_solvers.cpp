#include "MMLBase.h"

#include "core/LinAlgEqSolvers.h"

using namespace MML;


void Demo_Lin_alg_sys_solvers()
{
	// initialize a 4x4 matrix, defining the system of linear equations
	Matrix<Real> mat(4, 4, { -2,  4,  3,  5,
														3,  8, -7, -2,
														1, -3,  7,  4,
														5, -4,  9,  1 });
	Matrix<Real> matcopy(mat);
	std::cout << "Initial matrix:\n";  mat.Print(std::cout, 10, 3);

	// right side vector
	Vector<Real> rhs({ 3, -2, 1, 5 }), rhscopy(rhs);
	std::cout << "\nRight side: ";     rhs.Print(std::cout, 10, 3);

	// solution is returned in rhs vector
	GaussJordanSolver<Real>::Solve(mat, rhs);

	Vector<Real> sol2 = GaussJordanSolver<Real>::SolveConst(matcopy, rhscopy);

	std::cout << "\nSolution: " << rhs << std::endl;
	std::cout << "mat * sol = "; (matcopy * rhs).Print(std::cout, 10, 3);
}