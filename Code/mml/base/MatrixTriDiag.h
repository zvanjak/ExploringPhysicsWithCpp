#if !defined MML_MATRIX_TRIDIAG_H
#define MML_MATRIX_TRIDIAG_H

#include "MMLBase.h"

namespace MML
{
	template<class Type>
	class TridiagonalMatrix
	{
	private:
		int _dim;
		Vector<Type> _belowDiag;
		Vector<Type> _diag;
		Vector<Type> _aboveDiag;

	public:
		TridiagonalMatrix(int dim, const Vector<Type>& a, const Vector<Type>& diag, const Vector<Type>& c) : _belowDiag(a), _diag(diag), _aboveDiag(c), _dim(dim)
		{
			if( dim < 3)
				throw("Error in TridiagonalMatrix constructor: dim < 3");

			if (a.size() != dim || diag.size() != dim || c.size() != dim)
				throw("Error in TridiagonalMatrix constructor: wrong dimensions");
		}

		TridiagonalMatrix(int dim, std::initializer_list<Type> values) : _dim(dim), _belowDiag(dim), _diag(dim), _aboveDiag(dim)
		{
			if (values.size() != dim * 3 - 2)
				throw("Error in TridiagonalMatrix constructor: wrong dimensions");

			auto val = values.begin();
			_belowDiag[0] = 0.0;
			_diag[0] = *val++;
			_aboveDiag[0] = *val++;
			for (int i = 1; i < dim - 1; ++i)
			{
				_belowDiag[i] = *val++;
				_diag[i] = *val++;
				_aboveDiag[i] = *val++;
			}
			_belowDiag[dim - 1] = *val++;
			_diag[dim - 1] = *val++;
			_aboveDiag[dim - 1] = 0.0;
		}

		int RowNum() const { return _dim; }
		int ColNum() const { return _dim; }

		Type  operator()(int i, int j) const {
			if (i == j)
				return _diag[i];
			else if (i == j - 1)
				return _aboveDiag[i];
			else if (i == j + 1 && j < _dim - 1)
				return _belowDiag[i];
			else
				return 0.0;
		}
		Type& operator()(int i, int j) {
			if (i == j)
				return _diag[i];
			else if (i == j - 1)
				return _aboveDiag[i];
			else if (i == j + 1 && j < _dim - 1)
				return _belowDiag[i];
			else
				throw MatrixAccessBoundsError("TridiagonalMatrix::operator()", i, j, _dim, _dim);
		}

		// TODO 0.9 - dodati IsEqual
		// TODO 1.0 - imaju li smisla operacije s regularnim matricama?
		TridiagonalMatrix operator+(const TridiagonalMatrix& b) const
		{
			if (_dim != b._dim)
				throw MatrixDimensionError("TridiagonalMatrix::operator+() - must be same dim", _dim, _dim, b._dim, b._dim);

			TridiagonalMatrix temp(*this);
			for (int i = 0; i < _belowDiag.size(); i++)
				temp._belowDiag[i] += b._belowDiag[i];
			for (int i = 0; i < _diag.size(); i++)
				temp._diag[i] += b._diag[i];
			for (int i = 0; i < _aboveDiag.size(); i++)
				temp._aboveDiag[i] += b._aboveDiag[i];

			return temp;
		}
		TridiagonalMatrix operator-(const TridiagonalMatrix& b) const
		{
			if (_dim != b._dim)
				throw MatrixDimensionError("TridiagonalMatrix::operator+() - must be same dim", _dim, _dim, b._dim, b._dim);

			TridiagonalMatrix temp(*this);
			for (int i = 0; i < _belowDiag.size(); i++)
				temp._belowDiag[i] -= b._belowDiag[i];
			for (int i = 0; i < _diag.size(); i++)
				temp._diag[i] -= b._diag[i];
			for (int i = 0; i < _aboveDiag.size(); i++)
				temp._aboveDiag[i] -= b._aboveDiag[i];

			return temp;
		}

		friend TridiagonalMatrix operator*(const TridiagonalMatrix& a, Type b)
		{
			int	i, j;
			TridiagonalMatrix	ret(a);

			for (int i = 0; i < ret._belowDiag.size(); i++)
				ret._belowDiag[i] *= b;
			for (int i = 0; i < ret._diag.size(); i++)
				ret._diag[i] *= b;
			for (int i = 0; i < ret._aboveDiag.size(); i++)
				ret._aboveDiag[i] *= b;

			return ret;
		}
		friend TridiagonalMatrix operator*(Type a, const TridiagonalMatrix& b)
		{
			int	i, j;
			TridiagonalMatrix	ret(b);

			for (int i = 0; i < ret._belowDiag.size(); i++)
				ret._belowDiag[i] *= a;
			for (int i = 0; i < ret._diag.size(); i++)
				ret._diag[i] *= a;
			for (int i = 0; i < ret._aboveDiag.size(); i++)
				ret._aboveDiag[i] *= a;

			return ret;
		}
		friend TridiagonalMatrix operator/(const TridiagonalMatrix& a, Type b)
		{
			int	i, j;
			TridiagonalMatrix	ret(a);

			for (int i = 0; i < ret._belowDiag.size(); i++)
				ret._belowDiag[i] /= b;
			for (int i = 0; i < ret._diag.size(); i++)
				ret._diag[i] /= b;
			for (int i = 0; i < ret._aboveDiag.size(); i++)
				ret._aboveDiag[i] /= b;

			return ret;
		}

		// TODO - invert
		// TODO - 0.9 transpose

		void Solve(Vector<Type>& rhs, Vector<Type>& sol)
		{
			int j, n = _belowDiag.size();
			Real bet;
			Vector<Type> gam(n);

			if (_diag[0] == 0.0)
				throw("Error 1 in tridag");

			sol[0] = rhs[0] / (bet = _diag[0]);

			for (j = 1; j < n; j++) {
				gam[j] = _aboveDiag[j - 1] / bet;
				bet = _diag[j] - _belowDiag[j] * gam[j];

				if (bet == 0.0)
					throw("Error 2 in tridag");

				sol[j] = (rhs[j] - _belowDiag[j] * sol[j - 1]) / bet;
			}
			for (j = (n - 2); j >= 0; j--)
				sol[j] -= gam[j + 1] * sol[j + 1];
		}
		Vector<Type> Solve(Vector<Type>& rhs)
		{
			Vector<Type> sol(rhs.size());
			Solve(rhs, sol);
			return sol;
		}

		void   Print(std::ostream& stream, int width, int precision) const
		{
			stream << "Dim: " << _dim << std::endl;
			for (int i = 0; i < _dim; i++)
			{
				stream << "[ ";
				for (int j = 0; j < _dim; j++) {
					stream << std::setw(width) << std::setprecision(precision) << (*this)(i, j) << ", ";
				}
				stream << " ]" << std::endl;
			}
		}
	};
}

#endif