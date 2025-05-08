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
		static bool SolveInPlace(Matrix<Type>& a, Matrix<Type>& b)
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
		static bool SolveInPlace(Matrix<Type>& a, Vector<Type>& b)
		{
			auto bmat = Utils::ColumnMatrixFromVector(b);

			bool ret = SolveInPlace(a, bmat);
			b = bmat.VectorFromColumn(0);
			return ret;
		}

		// solving for a given RHS vector, but with return value
		// (in case of singular matrix 'a', exception is thrown)
		static Vector<Type> Solve(Matrix<Type>& a, const Vector<Type>& b)
		{
			Matrix<Type> bmat = Utils::ColumnMatrixFromVector(b);

			bool success = SolveInPlace(a, bmat);
			
			if (!success)
				throw SingularMatrixError("GaussJordanSolver::Solve - Singular Matrix");
			
			return bmat.VectorFromColumn(0);
		}

		// solving for a given RHS vector, but with return value
		// both matrix and vector are const and are not changed
		// (in case of singular matrix 'a', exception is thrown)
		static Vector<Type> SolveConst(const Matrix<Type>& a, const Vector<Type>& b)
		{
			if (a.RowNum() != a.ColNum())
				throw MatrixDimensionError("GaussJordanSolver::SolveConst - matrix must be square", a.RowNum(), a.ColNum(), -1, -1);
			if (a.RowNum() != b.size())
				throw MatrixDimensionError("GaussJordanSolver::SolveConst - matrix and vector must have same size", a.RowNum(), a.ColNum(), b.size(), -1);
			if (b.size() == 0)
				throw MatrixDimensionError("GaussJordanSolver::SolveConst - vector must be non-empty", a.RowNum(), a.ColNum(), b.size(), -1);
			if (a.RowNum() == 0 || a.ColNum() == 0 )
				throw MatrixDimensionError("GaussJordanSolver::SolveConst - matrix must be non-empty", a.RowNum(), a.ColNum(), -1, -1);

			Matrix<Real> mat(a);
			Matrix<Type> bmat = Utils::ColumnMatrixFromVector(b);

			bool success = SolveInPlace(mat, bmat);

			if (!success)
				throw SingularMatrixError("GaussJordanSolver::SolveConst - Singular Matrix");

			return bmat.VectorFromColumn(0);
		}
	};
}
// end namespace
#endif