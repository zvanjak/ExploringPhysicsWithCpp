#if !defined MML_MATRIX_SYM_H
#define MML_MATRIX_SYM_H

#include "MMLBase.h"

#include "Matrix.h"

namespace MML
{
	template<class Type>
	class MatrixSym
	{
	private:
		int  _dim;
		Type** _ptrData;

		void Init(int dim)
		{
			if (dim <= 0)
				throw MatrixDimensionError("MatrixSym::Init - dimension must be positive", dim, -1, -1, -1);

			_dim = dim;

			_ptrData = new Type * [dim];

			int numElem = (dim * dim + dim) / 2;
			if (_ptrData)
			{
				_ptrData[0] = numElem > 0 ? new Type[numElem] : nullptr;

				if (numElem > 0 && _ptrData[0] == nullptr)
					throw MatrixAllocationError("MatrixSym::Init - allocation error", dim, -1);

				for (int i = 1; i < dim; i++)
					_ptrData[i] = _ptrData[i - 1] + i;
			}
			else
				throw MatrixAllocationError("Matrix::Init - allocation error", dim, -1);
		}

	public:
		typedef Type value_type;      // make Type available externally

		///////////////////////          Constructors and destructor       //////////////////////
		MatrixSym() : _dim(0), _ptrData{ nullptr } {}
		MatrixSym(int dim) : _dim(dim)
		{
			Init(dim);
			for (int i = 0; i < _dim; ++i)
				for (int j = 0; j < _dim; ++j)
					_ptrData[i][j] = 0;
		}
		MatrixSym(int dim, Type val) : _dim(dim)
		{
			Init(dim);
			for (int i = 0; i < _dim; ++i)
				for (int j = 0; j <= i; ++j)
					_ptrData[i][j] = val;
		}
		MatrixSym(int dim, std::initializer_list<Type> values) : _dim(dim)
		{
			Init(dim);
			if (values.size() != (dim * dim + dim) / 2)
				throw MatrixDimensionError("Error in MatrixSym constructor: wrong dimensions", dim, -1, -1, -1);

			auto val = values.begin();
			for (int i = 0; i < _dim; ++i)
				for (int j = 0; j <= i; ++j)
					if (val != values.end())
						_ptrData[i][j] = *val++;
		}
		MatrixSym(const MatrixSym& m) : _dim(m._dim)
		{
			Init(m._dim);

			for (int i = 0; i < _dim; ++i)
				for (int j = 0; j <= i; ++j)
					_ptrData[i][j] = m._ptrData[i][j];
		}
		MatrixSym(MatrixSym&& m)
		{
			_ptrData = m._ptrData;

			_dim = m._dim;

			m._dim = 0;
			m._ptrData = nullptr;
		}
		~MatrixSym()
		{
			if (_ptrData != NULL) {
				delete[](_ptrData[0]);
				delete[](_ptrData);
			}
		}

		////////////////////////               Standard stuff             ////////////////////////
		int Dim() const { return (int)_dim; }
		int RowNum() const { return (int)_dim; }
		int ColNum() const { return (int)_dim; }

		bool IsEqual(const MatrixSym& b, Type eps = Defaults::MatrixIsEqualTolerance) const
		{
			if (Dim() != b.Dim())
				return false;

			for (int i = 0; i < Dim(); i++)
				for (int j = 0; j <= i; j++)
				{
					if (Abs(_ptrData[i][j] - b._ptrData[i][j]) > eps)
						return false;
				}

			return true;
		}
		static bool AreEqual(const MatrixSym& a, const MatrixSym& b, Type eps = Defaults::MatrixIsEqualTolerance)
		{
			return a.IsEqual(b, eps);
		}

		Matrix<Type> GetAsMatrix() const
		{
			Matrix<Type> ret(Dim(), Dim());

			for (int i = 0; i < Dim(); i++)
				for (int j = 0; j < Dim(); j++)
					ret[i][j] = (*this)(i, j);

			return ret;
		}

		/////////////////////          Init from regular Matrix           /////////////////////
		MatrixSym InitFromLower(Matrix<Type>& b)
		{
			if (b.RowNum() != b.ColNum())
				throw MatrixDimensionError("MatrixSym::InitFromLower - must be square matrix", b.RowNum(), b.ColNum(), -1, -1);

			MatrixSym ret(b.RowNum());
			for (int i = 0; i < b.RowNum(); i++)
				for (int j = 0; j <= i; j++)
					ret._ptrData[i][j] = b[i][j];

			return ret;
		}
		MatrixSym InitFromUpper(Matrix<Type>& b)
		{
			if (b.RowNum() != b.ColNum())
				throw MatrixDimensionError("MatrixSym::InitFromUpper - must be square matrix", b.RowNum(), b.ColNum(), -1, -1);

			MatrixSym ret(b.RowNum());
			for (int i = 0; i < b.RowNum(); i++)
				for (int j = i; j < b.RowNum(); j++)
					ret._ptrData[i][j] = b[i][j];

			return ret;
		}

		/////////////////////          Vector-Matrix conversion           /////////////////////
		Vector<Type> VectorFromRow(int rowInd)
		{
			if (rowInd >= RowNum())
				throw MatrixAccessBoundsError("VectorFromRow - row index must be less then a.RowNum()", rowInd, 0, RowNum(), ColNum());

			Vector<Type> ret(ColNum());
			for (int i = 0; i < ColNum(); i++)
				ret[i] = this->a(rowInd, i);

			return ret;
		}
		Vector<Type> VectorFromColumn(int colInd)
		{
			if (colInd >= ColNum())
				throw MatrixAccessBoundsError("VectorFromColumn - column index must be less then a.ColNum()", 0, colInd, RowNum(), ColNum());

			Vector<Type> ret(RowNum());
			for (int i = 0; i < RowNum(); i++)
				ret[i] = this->a(i, colInd);

			return ret;
		}
		Vector<Type> VectorFromDiagonal()
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("VectorFromDiagonal - must be square matrix", RowNum(), ColNum(), -1, -1);

			Vector<Type> ret(RowNum());
			for (int i = 0; i < RowNum(); i++)
				ret[i] = this->a(i, i);

			return ret;
		}

		///////////////////////////            Operators             ///////////////////////////
		MatrixSym& operator=(const MatrixSym& m)
		{
			if (this == &m)
				return *this;

			if (_dim != m._dim)
			{
				delete[](_ptrData[0]);
				delete[](_ptrData);

				Init(m.Dim());
			}

			for (size_t i = 0; i < m.Dim(); ++i)
				for (size_t j = 0; j <= i; ++j)
					_ptrData[i][j] = m._ptrData[i][j];

			return *this;
		}
		MatrixSym& operator=(MatrixSym&& m)
		{
			if (this == &m)
				return *this;

			std::swap(_ptrData, m._ptrData);
			std::swap(_dim, m._dim);

			return *this;
		}

		////////////////////////           Access operators             ///////////////////////
		Type  operator()(int i, int j) const {
			if (i < j)
				return _ptrData[j][i];
			else
				return _ptrData[i][j];
		}
		Type& operator()(int i, int j) {
			if (i < j)
				return _ptrData[j][i];
			else
				return _ptrData[i][j];
		}

		// version with checking bounds
		Type  ElemAt(int i, int j) const
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("MatrixSym::ElemAt", i, j, RowNum(), ColNum());

			return operator()(i, j);
		}
		Type& ElemAt(int i, int j)
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("MatrixSym::ElemAt", i, j, RowNum(), ColNum());

			return operator()(i, j);
		}

		///////////////////////           Arithmetic operators             //////////////////////
		MatrixSym operator+(const MatrixSym& b) const
		{
			if (_dim != b._dim)
				throw MatrixDimensionError("MatrixSym::operator+() - must be same dim", _dim, -1, b._dim, -1);

			MatrixSym temp(_dim);
			for (size_t i = 0; i < Dim(); i++)
				for (size_t j = 0; j <= i; j++)
					temp._ptrData[i][j] = b._ptrData[i][j] + _ptrData[i][j];

			return temp;
		}
		MatrixSym operator-()
		{
			MatrixSym temp(_dim);
			for (size_t i = 0; i < Dim(); i++)
				for (size_t j = 0; j <= i; j++)
					temp._ptrData[i][j] = -_ptrData[i][j];

			return temp;
		}
		MatrixSym operator-(const MatrixSym& b) const
		{
			if (_dim != b._dim)
				throw MatrixDimensionError("MatrixSym::operator-() - must be same dim", _dim, -1, b._dim, -1);

			MatrixSym temp(_dim);
			for (size_t i = 0; i < Dim(); i++)
				for (size_t j = 0; j <= i; j++)
					temp._ptrData[i][j] = b._ptrData[i][j] - _ptrData[i][j];

			return temp;
		}

		Matrix<Type> operator*(const MatrixSym& b) const
		{
			if (Dim() != b.Dim())
				throw MatrixDimensionError("MatrixSym::operator*(MatrixSym &) - a.colNum must be equal to b.rowNum", _dim, _dim, b._dim, b._dim);

			Matrix<Type>	ret(RowNum(), b.ColNum());
			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++)
				{
					ret[i][j] = 0;
					for (int k = 0; k < ColNum(); k++)
						ret[i][j] += (*this)(i, k) * b(k, j);
				}

			return	ret;
		}

		Matrix<Type> operator*(const Matrix<Type>& b) const
		{
			if (Dim() != b.RowNum())
				throw MatrixDimensionError("MatrixSym::operator*(Matrix &) - a.colNum must be equal to b.rowNum", _dim, _dim, b.RowNum(), b.ColNum());

			Matrix<Type>	ret(RowNum(), b.ColNum());
			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++)
				{
					ret[i][j] = 0;
					for (int k = 0; k < ColNum(); k++)
						ret[i][j] += (*this)(i, k) * b(k, j);
				}

			return	ret;
		}

		friend MatrixSym operator*(const MatrixSym& a, Type b)
		{
			int	i, j;
			MatrixSym	ret(a.Dim());

			for (i = 0; i < a.Dim(); i++)
				for (j = 0; j <= i; j++)
					ret[i][j] = a._ptrData[i][j] * b;

			return ret;
		}
		friend MatrixSym operator*(Type a, const MatrixSym& b)
		{
			int	i, j;
			MatrixSym	ret(a.Dim());

			for (i = 0; i < a.Dim(); i++)
				for (j = 0; j <= i; j++)
					ret[i][j] = a * b._ptrData[i][j];

			return ret;
		}
		friend MatrixSym operator/(const MatrixSym& a, Type b)
		{
			int	i, j;
			MatrixSym	ret(a.Dim());

			for (i = 0; i < a.Dim(); i++)
				for (j = 0; j <= i; j++)
					ret[i][j] = a._ptrData[i][j] / b;

			return ret;
		}

		friend Vector<Type> operator*(const MatrixSym& a, const Vector<Type>& b)
		{
			if (a.Dim() != b.size())
				throw MatrixDimensionError("operator*(MatSym a, Vec b) - a.Dim must be equal to vector size", a.Dim(), a.Dim(), (int)b.size(), -1);

			Vector<Type>	ret(a.RowNum());
			for (int i = 0; i < a.Dim(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < a.Dim(); j++)
					ret[i] += a(i, j) * b[j];
			}

			return ret;
		}
		friend Vector<Type> operator*(const Vector<Type>& a, const MatrixSym& b)
		{
			if (a.size() != b.Dim())
			{
				//std::string error = std::format("Hello {}!\n", "world");
				throw MatrixDimensionError("operator*(Vec a, MatSym b) - vector size must be equal to b.Dim", (int)a.size(), -1, b.Dim(), -1);
			}

			Vector<Type>	ret(b.Dim());
			for (int i = 0; i < b.Dim(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < b.Dim(); j++)
					ret[i] += a[j] * b(i, j);
			}

			return ret;
		}

		////////////////////////            Inverse              ///////////////////////
		Matrix<Type> GetInverse() const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixSym::GetInverse - must be square matrix", _dim, _dim, -1, -1);

			Matrix<Type> a = this->GetAsMatrix();              // making a copy, where inverse will be stored at the end

			a.Invert();

			return a;
		}

		///////////////////////////               I/O                 ///////////////////////////
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}
		void   Print(std::ostream& stream, int width, int precision) const
		{
			stream << "Rows: " << Dim() << std::endl;

			for (int i = 0; i < Dim(); i++)
			{
				stream << "[ ";
				for (int j = 0; j < Dim(); j++)
				{
					stream << std::setw(width) << std::setprecision(precision) << (*this)(i, j) << ", ";
				}
				stream << " ]" << std::endl;
			}
		}
		friend std::ostream& operator<<(std::ostream& stream, const MatrixSym& a)
		{
			a.Print(stream, 10, 3);

			return stream;
		}
	};
}
#endif