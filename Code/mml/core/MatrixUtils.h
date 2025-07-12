#if !defined MML_MATRIXUTILS_H
#define MML_MATRIXUTILS_H

#include "MMLBase.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/Matrix.h"
#include "base/MatrixNM.h"
#include "base/Polynom.h"

namespace MML
{
	namespace Utils
	{
		void FaddeevAlg(const Matrix<Real>& A, PolynomRealFunc& outCharPoly, Real& outDet, Matrix<Real>& outInv)
		{
			int N = A.RowNum();

			Matrix<Real> nullMat(N, N);		// null matrix
			auto identMat = Matrix<Real>::GetUnitMatrix(N);

			Vector<Matrix<Real>> M(N + 1, nullMat);

			outCharPoly.SetDegree(N);

			M[0] = nullMat;
			outCharPoly[N] = 1.0;
			for (int k = 1; k <= A.RowNum(); k++)
			{
				M[k] = A * M[k - 1] + outCharPoly[N - k + 1] * identMat;
				outCharPoly[N - k] = -1.0 / k * (A * M[k]).Trace();
			}

			outDet = std::pow(-1.0, N) * outCharPoly[0];
			outInv = M[N] / outDet; 
		}
	}
} // namespace MML
#endif


