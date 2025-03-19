#if !defined  MML_LINEAR_ALG_EQ_SOLVERS_H
#define MML_LINEAR_ALG_EQ_SOLVERS_H

#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/BaseUtils.h"

namespace MML
{
	///////////////////////   GAUSS-JORDAN SOLVER    /////////////////////////////
	template<class Type>
	class GaussJordanSolver
	{
	public:
		// solving with Matrix RHS (ie. solving simultaneously for multiple RHS)
		static bool Solve(Matrix<Type>& a, Matrix<Type>& b)
		{
			int i, icol, irow, j, k, l, ll;
			Real big;
			Type dum, pivinv;

			int n = a.RowNum();
			int m = b.ColNum();
			std::vector<int> indxc(n), indxr(n), ipiv(n);
			for (j = 0; j < n; j++) ipiv[j] = 0;
			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j < n; j++)
					if (ipiv[j] != 1)
						for (k = 0; k < n; k++) {
							if (ipiv[k] == 0) {
								if (Abs(a[j][k]) >= big) {
									big = Abs(a[j][k]);
									irow = j;
									icol = k;
								}
							}
						}
				++(ipiv[icol]);
				if (irow != icol) {
					for (l = 0; l < n; l++) std::swap(a[irow][l], a[icol][l]);
					for (l = 0; l < m; l++) std::swap(b[irow][l], b[icol][l]);
				}
				indxr[i] = irow;
				indxc[i] = icol;

				if (a[icol][icol] == Real{ 0.0 })
					return false;
				// throw SingularMatrixError("GaussJordanSolver::Solve - Singular Matrix");

				pivinv = Real{ 1.0 } / a[icol][icol];
				a[icol][icol] = 1.0;
				for (l = 0; l < n; l++) a[icol][l] *= pivinv;
				for (l = 0; l < m; l++) b[icol][l] *= pivinv;
				for (ll = 0; ll < n; ll++)
					if (ll != icol) {
						dum = a[ll][icol];
						a[ll][icol] = 0.0;
						for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
						for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
					}
			}
			for (l = n - 1; l >= 0; l--) {
				if (indxr[l] != indxc[l])
					for (k = 0; k < n; k++)
						std::swap(a[k][indxr[l]], a[k][indxc[l]]);
			}

			return true;
		}
		
		// solving for a given RHS vector
		static bool Solve(Matrix<Type>& a, Vector<Type>& b)
		{
			auto bmat = Utils::ColumnMatrixFromVector(b);

			bool ret = Solve(a, bmat);
			b = bmat.VectorFromColumn(0);
			return ret;
		}

		// solving for a given RHS vector, but with return value
		// (in case of singular matrix 'a', exception is thrown)
		static Vector<Type> SolveConst(Matrix<Type>& a, const Vector<Type>& b)
		{
			Vector<Type> ret(b.size());
			Matrix<Type> bmat = Utils::ColumnMatrixFromVector(b);
			bool success = Solve(a, bmat);
			if (!success)
				throw SingularMatrixError("GaussJordanSolver::Solve - Singular Matrix");
			ret = bmat.VectorFromColumn(0);
			return ret;
		}
	};
}
// end namespace
#endif