#if !defined  MML_MATRIX_H
#define MML_MATRIX_H

#include "MMLBase.h"

#include "Vector.h"

namespace MML
{
	template<class Type>
	class Matrix
	{
	private:
		int  _rows;
		int  _cols;
		Type** _data;

		void Init(int rows, int cols)
		{
			if (rows <= 0 || cols < 0)
				throw MatrixDimensionError("Matrix::Init - rowNum and colNum must be positive", rows, cols, -1, -1);

			_rows = rows;
			_cols = cols;
			int numElem = rows * cols;

			_data = new Type * [rows];
			if (_data)
			{
				_data[0] = numElem > 0 ? new Type[numElem] : nullptr;

				if (numElem > 0 && _data[0] == nullptr)
					throw MatrixAllocationError("Matrix::Init - allocation error", rows, cols);

				for (int i = 1; i < rows; i++)
					_data[i] = _data[i - 1] + cols;
			}
			else
				throw MatrixAllocationError("Matrix::Init - allocation error", rows, cols);
		}

	public:
		typedef Type value_type;      // make Type available externally

		///////////////////////          Constructors and destructor       //////////////////////
		explicit Matrix() : _rows(0), _cols(0), _data{ nullptr } {}
		explicit Matrix(int rows, int cols) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = 0;
		}
		explicit Matrix(int rows, int cols, const Type& val) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = val;
		}
		// useful if you have a pointer to continuous 2D array (can be in row-, or column-wise memory layout)
		explicit Matrix(int rows, int cols, Type* val, bool isRowWise = true) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);

			if (isRowWise)
				for (int i = 0; i < _rows; ++i)
					for (int j = 0; j < _cols; ++j)
						_data[i][j] = *val++;
			else
				for (int j = 0; j < _cols; ++j)
					for (int i = 0; i < _rows; ++i)
						_data[i][j] = *val++;
		}
		// in strict mode, you must supply ALL necessary values for complete matrix initialization
		explicit Matrix(int rows, int cols, std::initializer_list<Type> values, bool strictMode = true) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);

			auto val = values.begin();
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					if (val != values.end())
						_data[i][j] = *val++;
					else {
						if (strictMode)
							throw MatrixDimensionError("Matrix::Matrix - not enough values in initializer list", _rows, _cols, -1, -1);
						else
							_data[i][j] = Type{ 0 };
					}
		}

		Matrix(const Matrix& m) : _rows(m._rows), _cols(m._cols)
		{
			Init(m._rows, m._cols);

			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = m._data[i][j];
		}
		// creating submatrix from given matrix 'm'
		Matrix(const Matrix& m, int ind_row, int ind_col, int row_num, int col_num)
		{
			if (ind_row < 0 || ind_row >= m._rows || ind_col < 0 || ind_col >= m._cols)
				throw MatrixDimensionError("Matrix::Matrix - invalid row or column index", m._rows, m._cols, ind_row, ind_col);

			if (row_num <= 0 || col_num <= 0)
				throw MatrixDimensionError("Matrix::Matrix - rowNum and colNum must be positive", row_num, col_num, -1, -1);

			if (ind_row + row_num > m._rows || ind_col + col_num > m._cols)
				throw MatrixDimensionError("Matrix::Matrix - submatrix out of bounds", m._rows, m._cols, ind_row, ind_col);

			_rows = row_num;
			_cols = col_num;

			Init(row_num, col_num);

			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = m._data[ind_row + i][ind_col + j];
		}
		Matrix(Matrix&& m)
		{
			_data = m._data;

			_rows = m._rows;
			_cols = m._cols;

			m._rows = 0;
			m._cols = 0;
			m._data = nullptr;
		}
		~Matrix()
		{
			if (_data != NULL) {
				delete[](_data[0]);
				delete[](_data);
			}
		}

		void Resize(int rows, int cols, bool preserveElements = false)
		{
			if (preserveElements == true)
			{
				if (rows <= 0 || cols <= 0)
					throw MatrixDimensionError("Matrix::Resize - rowNum and colNum must be positive", rows, cols, -1, -1);
				
				Type** oldData = _data;
				int oldRows = _rows;
				int oldCols = _cols;
				
				Init(rows, cols);
				
				for (int i = 0; i < std::min(_rows, oldRows); ++i)
					for (int j = 0; j < std::min(_cols, oldCols); ++j)
						_data[i][j] = oldData[i][j];
				
				delete[] oldData[0];
				delete[] oldData;
			}
			else
			{
				if (rows <= 0 || cols <= 0)
					throw MatrixDimensionError("Matrix::Resize - rowNum and colNum must be positive", rows, cols, -1, -1);

				if (rows == RowNum() && cols == ColNum())      // nice :)
					return;

				if (_data != NULL) {
					delete[](_data[0]);
					delete[](_data);
				}

				Init(rows, cols);

				for (int i = 0; i < _rows; ++i)
					for (int j = 0; j < _cols; ++j)
						_data[i][j] = 0;
			}
		}   

		///////////////////////              Standard stuff                //////////////////////
		inline int RowNum() const { return _rows; }
		inline int ColNum() const { return _cols; }

		static Matrix GetUnitMatrix(int dim)
		{
			if (dim <= 0)
				throw MatrixDimensionError("Matrix::GetUnitMatrix - dimension must be positive", dim, dim, -1, -1);

			Matrix unitMat(dim, dim);
			unitMat.MakeUnitMatrix();

			return unitMat;
		}
		void   MakeUnitMatrix(void)
		{
			if (_rows == _cols)
			{
				for (int i = 0; i < _rows; i++)
					for (int j = 0; j < _cols; j++)
						if (i == j)
							_data[i][j] = 1;
						else
							_data[i][j] = 0;
			}
			else
				throw MatrixDimensionError("Matrix::MakeUnitMatrix - must be square matrix", _rows, _cols, -1, -1);
		}

		Matrix GetLower(bool includeDiagonal = true) const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::GetLower - must be square matrix", _rows, _cols, -1, -1);

			Matrix ret(RowNum(), ColNum());
			for (int i = 0; i < RowNum(); i++)
			{
				if (includeDiagonal)
					for (int j = 0; j <= i; j++)
						ret[i][j] = _data[i][j];
				else
					for (int j = 0; j < i; j++)
						ret[i][j] = _data[i][j];
			}

			return ret;
		}
		Matrix GetUpper(bool includeDiagonal = true) const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::GetUpper - must be square matrix", _rows, _cols, -1, -1);

			Matrix ret(RowNum(), ColNum());
			for (int i = 0; i < RowNum(); i++)
			{
				if (includeDiagonal)
					for (int j = i; j < ColNum(); j++)
						ret[i][j] = _data[i][j];
				else
					for (int j = i + 1; j < ColNum(); j++)
						ret[i][j] = _data[i][j];
			}

			return ret;
		}

		void InitRowWithVector(int rowInd, const Vector<Type>& vec)
		{
			if (rowInd < 0 || rowInd >= RowNum())
				throw MatrixAccessBoundsError("Matrix::InitRowWithVector - invalid row index", rowInd, 0, RowNum(), ColNum());
			if (vec.size() != ColNum())
				throw MatrixDimensionError("Matrix::InitRowWithVector - vector size must match column number", RowNum(), ColNum(), -1, -1);
			
			for (int j = 0; j < ColNum(); j++)
				_data[rowInd][j] = vec[j];
		}
		void InitColWithVector(int colInd, const Vector<Type>& vec)
		{
			if (colInd < 0 || colInd >= ColNum())
				throw MatrixAccessBoundsError("Matrix::InitColumnWithVector - invalid column index", 0, colInd, RowNum(), ColNum());
			if (vec.size() != RowNum())
				throw MatrixDimensionError("Matrix::InitColumnWithVector - vector size must match row number", RowNum(), ColNum(), -1, -1);
			
			for (int i = 0; i < RowNum(); i++)
				_data[i][colInd] = vec[i];
		}
		
		void SwapRows(int k, int l)
		{
			if (k < 0 || k >= RowNum() || l < 0 || l >= RowNum())
				throw MatrixDimensionError("Matrix::SwapRows - invalid row index", RowNum(), ColNum(), k, l);
			if (k == l)
				return;
			for (int j = 0; j < ColNum(); j++)
				std::swap(_data[k][j], _data[l][j]);
		}
		void SwapCols(int k, int l)
		{
			if (k < 0 || k >= ColNum() || l < 0 || l >= ColNum())
				throw MatrixDimensionError("Matrix::SwapCols - invalid column index", RowNum(), ColNum(), k, l);
			if (k == l)
				return;
			for (int i = 0; i < RowNum(); i++)
				std::swap(_data[i][k], _data[i][l]);
		}

		///////////////////////          Matrix to Vector conversions      //////////////////////
		Vector<Type> VectorFromRow(int rowInd) const
		{
			if (rowInd < 0 || rowInd >= RowNum())
				throw MatrixAccessBoundsError("VectorFromRow - invalid row index", rowInd, 0, RowNum(), ColNum());

			Vector<Type> ret(ColNum());
			for (int i = 0; i < ColNum(); i++)
				ret[i] = (*this)(rowInd, i);

			return ret;
		}
		Vector<Type> VectorFromColumn(int colInd) const
		{
			if (colInd < 0 || colInd >= ColNum())
				throw MatrixAccessBoundsError("VectorFromColumn - invalid column index", 0, colInd, RowNum(), ColNum());

			Vector<Type> ret(RowNum());
			for (int i = 0; i < RowNum(); i++)
				ret[i] = (*this)(i, colInd);

			return ret;
		}
		Vector<Type> VectorFromDiagonal() const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("VectorFromDiagonal - must be square matrix", RowNum(), ColNum(), -1, -1);

			Vector<Type> ret(RowNum());
			for (int i = 0; i < RowNum(); i++)
				ret[i] = (*this)(i, i);

			return ret;
		}

		///////////////////////               Matrix properties            //////////////////////
		bool IsUnit(double eps = Defaults::IsMatrixUnitPrecision) const
		{
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					if (i != j && Abs((*this)[i][j]) > eps)
						return false;

			for (int i = 0; i < RowNum(); i++)
				if (Abs((*this)[i][i] - Real{ 1.0 }) > eps)
					return false;

			return true;
		}
		bool IsDiagonal(double eps = Defaults::IsMatrixDiagonalPrecision) const
		{
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					if (i != j && Abs((*this)[i][j]) > eps)
						return false;

			return true;
		}
		bool IsDiagDominant() const
		{
			for (int i = 0; i < RowNum(); i++)
			{
				Type sum = 0.0;
				for (int j = 0; j < ColNum(); j++)
					if (i != j)
						sum += Abs((*this)[i][j]);

				if (Abs((*this)[i][i]) < sum)
					return false;
			}
			return true;
		}
		bool IsSymmetric() const
		{
			if (RowNum() != ColNum())
				return false;

			for (int i = 0; i < RowNum(); i++)
				for (int j = i + 1; j < ColNum(); j++)
					if ((*this)[i][j] != (*this)[j][i])
						return false;

			return true;
		}
		bool IsAntiSymmetric() const
		{
			if (RowNum() != ColNum())
				return false;
			for (int i = 0; i < RowNum(); i++)
				for (int j = i + 1; j < ColNum(); j++)
					if ((*this)[i][j] != -(*this)[j][i])
						return false;
			return true;
		}

		Real NormL1() const
		{
			Real norm = 0;
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					norm += Abs(_data[i][j]);
			return norm;
		}
		Real NormL2() const
		{
			Real norm = 0;
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					norm += _data[i][j] * _data[i][j];
			return std::sqrt(norm);
		}
		Real NormLInf() const
		{
			Real norm = 0;
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					norm = std::max(norm, Abs(_data[i][j]));
			return norm;
		}

		///////////////////////             Assignment operators           //////////////////////
		Matrix& operator=(const Matrix& m)
		{
			if (this == &m)
				return *this;

			if (_rows != m._rows || _cols != m._cols)
			{
				delete[](_data[0]);
				delete[](_data);

				Init(m._rows, m._cols);
			}

			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = m._data[i][j];

			return *this;
		}
		Matrix& operator=(Matrix&& m) noexcept
		{
			if (this == &m)
				return *this;

			std::swap(_data, m._data);
			std::swap(_rows, m._rows);
			std::swap(_cols, m._cols);

			return *this;
		}

		///////////////////////               Access operators             //////////////////////
		inline Type* operator[](int i) { return _data[i]; }
		inline const Type* operator[](const int i) const { return _data[i]; }

		inline Type  operator()(int i, int j) const { return _data[i][j]; }
		inline Type& operator()(int i, int j) { return _data[i][j]; }

		// version with checked access
		Type  at(int i, int j) const
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("Matrix::at", i, j, RowNum(), ColNum());

			return _data[i][j];
		}
		Type& at(int i, int j)
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("Matrix::at", i, j, RowNum(), ColNum());

			return _data[i][j];
		}

		///////////////////////              Equality operations           //////////////////////
		bool operator==(const Matrix& b) const
		{
			if (_rows != b._rows || _cols != b._cols)
				return false;

			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					if (_data[i][j] != b._data[i][j])
						return false;

			return true;
		}
		bool operator!=(const Matrix& b) const
		{
			return !(*this == b);
		}

		bool IsEqualTo(const Matrix<Type>& b, Type eps = Defaults::MatrixEqualityPrecision) const
		{
			if (RowNum() != b.RowNum() || ColNum() != b.ColNum())
				return false;

			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
				{
					if (Abs(_data[i][j] - b._data[i][j]) > eps)
						return false;
				}

			return true;
		}
		static bool AreEqual(const Matrix& a, const Matrix& b, Type eps = Defaults::MatrixEqualityPrecision)
		{
			return a.IsEqualTo(b, eps);
		}

		///////////////////////              Arithmetic operators          //////////////////////
		Matrix  operator-() const            // unary minus
		{
			Matrix temp(_rows, _cols);
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					temp._data[i][j] = -_data[i][j];

			return temp;
		}

		Matrix  operator+(const Matrix& b) const
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator+() - must be same dim", _rows, _cols, b._rows, b._cols);

			Matrix temp(_rows, _cols);
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					temp._data[i][j] = b._data[i][j] + _data[i][j];

			return temp;
		}
		Matrix& operator+=(const Matrix& b)
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator+=() - must be same dim", _rows, _cols, b._rows, b._cols);
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					_data[i][j] += b._data[i][j];
			return *this;
		}
		Matrix  operator-(const Matrix& b) const
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator-() - must be same dim", _rows, _cols, b._rows, b._cols);

			Matrix temp(_rows, _cols);
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					temp._data[i][j] = _data[i][j] - b._data[i][j];

			return temp;
		}
		Matrix& operator-=(const Matrix& b)
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator-=() - must be same dim", _rows, _cols, b._rows, b._cols);
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					_data[i][j] -= b._data[i][j];
			return *this;
		}
		Matrix  operator*(const Matrix& b) const
		{
			if (ColNum() != b.RowNum())
				throw MatrixDimensionError("Matrix::operator*() - a.colNum must be equal to b.rowNum", _rows, _cols, b._rows, b._cols);

			Matrix	ret(RowNum(), b.ColNum());
			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++)
				{
					ret._data[i][j] = 0;
					for (int k = 0; k < ColNum(); k++)
						ret._data[i][j] += _data[i][k] * b._data[k][j];
				}

			return	ret;
		}

		Matrix  operator*(const Type& b) const
		{
			int	i, j;
			Matrix	ret(RowNum(), ColNum());

			for (i = 0; i < RowNum(); i++)
				for (j = 0; j < ColNum(); j++)
					ret[i][j] = _data[i][j] * b;

			return ret;
		}
		Matrix& operator*=(const Type& b)
		{
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					_data[i][j] *= b;
			return *this;
		}
		Matrix  operator/(const Type& b) const
		{
			Matrix	ret(RowNum(), ColNum());

			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					ret[i][j] = _data[i][j] / b;

			return ret;
		}
		Matrix& operator/=(const Type& b)
		{
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					_data[i][j] /= b;
			return *this;
		}

		Vector<Type> operator*(const Vector<Type>& b) const
		{
			if (ColNum() != b.size())
				throw MatrixDimensionError("operator*(Mat a, Vec b) - a.colNum must be equal to vector size", _rows, _cols, (int)b.size(), -1);

			Vector<Type>	ret(RowNum());
			for (int i = 0; i < RowNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < ColNum(); j++)
					ret[i] += _data[i][j] * b[j];
			}

			return ret;
		}

		friend Matrix operator*(const Type& a, const Matrix<Type>& b)
		{
			int	i, j;
			Matrix	ret(b.RowNum(), b.ColNum());

			for (i = 0; i < b.RowNum(); i++)
				for (j = 0; j < b.ColNum(); j++)
					ret[i][j] = a * b._data[i][j];

			return ret;
		}
		friend Vector<Type> operator*(const Vector<Type>& a, const Matrix<Type>& b)
		{
			if (a.size() != b.RowNum())
			{
				//std::string error = std::format("Hello {}!\n", "world");
				throw MatrixDimensionError("operator*(Vec a, Mat b) - vector size must be equal to b.rowNum", (int)a.size(), -1, b._rows, b._cols);
			}

			Vector<Type>	ret(b.ColNum());
			for (int i = 0; i < b.ColNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < b.RowNum(); j++)
					ret[i] += a[j] * b(j, i);
			}

			return ret;
		}

		///////////////////////            Trace, Inverse & Transpose      //////////////////////
		Type   Trace() const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::Trace - must be square matrix", _rows, _cols, -1, -1);

			Type sum = 0;
			for (int i = 0; i < RowNum(); i++)
				sum += _data[i][i];

			return sum;
		}

		void   Invert()
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::Invert - must be square matrix", _rows, _cols, -1, -1);

			Matrix& a = *this;
			Matrix  b(RowNum(), 1);      // dummy rhs

			b(0, 0) = 1.0;

			int i, icol, irow, j, k, l, ll;
			Type dum, pivinv;
			Real big;

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

				if (a[icol][icol] == 0.0)
					throw SingularMatrixError("Matrix::Invert, Singular Matrix");

				pivinv = 1.0 / a[icol][icol];
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
		}
		Matrix GetInverse() const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::GetInverse - must be square matrix", _rows, _cols, -1, -1);

			Matrix a(*this);              // making a copy, where inverse will be stored at the end
			a.Invert();

			return a;
		}

		void   Transpose()
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::Transpose - in-place Transpose possible only for square matrix", _rows, _cols, -1, -1);

			for (int i = 0; i < RowNum(); i++)
				for (int j = i + 1; j < ColNum(); j++)
					std::swap(_data[i][j], _data[j][i]);
		}
		Matrix GetTranspose() const
		{
			Matrix ret(ColNum(), RowNum());

			for (int i = 0; i < ColNum(); i++)
				for (int j = 0; j < RowNum(); j++)
					ret[i][j] = _data[j][i];

			return ret;
		}

		///////////////////////                    I/O                    //////////////////////
		void   Print(std::ostream& stream, int width, int precision) const
		{
			stream << "Rows: " << RowNum() << " Cols: " << ColNum();

			if (RowNum() == 0 || ColNum() == 0) {
				stream << " - Empty matrix" << std::endl;
				return;
			}
			else
				stream << std::endl;

			for (int i = 0; i < RowNum(); i++)
			{
				stream << "[ ";
				for (int j = 0; j < ColNum(); j++)
				{
					if (j == ColNum() - 1)
						stream << std::setw(width) << std::setprecision(precision) << _data[i][j];
					else
						stream << std::setw(width) << std::setprecision(precision) << _data[i][j] << ", ";
				}
				if (i == RowNum() - 1)
					stream << " ]";
				else
					stream << " ]" << std::endl;
			}
		}
		void   Print(std::ostream& stream, int width, int precision, Real zeroThreshold) const
		{
			stream << "Rows: " << RowNum() << " Cols: " << ColNum();

			if (RowNum() == 0 || ColNum() == 0) {
				stream << " - Empty matrix" << std::endl;
				return;
			}
			else
				stream << std::endl;

			for (int i = 0; i < RowNum(); i++)
			{
				stream << "[ ";
				for (int j = 0; j < ColNum(); j++)
				{
					Type value{ 0 };
					if (Abs(_data[i][j]) > zeroThreshold)
						value = _data[i][j];

					if (j == ColNum() - 1)
						stream << std::setw(width) << std::setprecision(precision) << value;
					else
						stream << std::setw(width) << std::setprecision(precision) << value << ", ";
				}
				if (i == RowNum() - 1)
					stream << " ]";
				else
					stream << " ]" << std::endl;
			}
		}

		friend std::ostream& operator<<(std::ostream& stream, const Matrix& a)
		{
			a.Print(stream, 10, 3);

			return stream;
		}
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}

		static bool LoadFromFile(std::string inFileName, Matrix& outMat)
		{
			std::ifstream file(inFileName);

			if (file.is_open())
			{
				int rows, cols;
				file >> rows >> cols;

				outMat.Resize(rows, cols);
				for (int i = 0; i < outMat.RowNum(); i++)
					for (int j = 0; j < outMat.ColNum(); j++)
						file >> outMat[i][j];

				file.close();
			}
			else {
				std::cerr << "Error: could not open file " << inFileName << " for reading." << std::endl;
				return false;
			}

			return true;
		}
		static bool SaveToFile(const Matrix& mat, std::string inFileName)
		{
			std::ofstream file(inFileName);

			if (file.is_open())
			{
				file << mat.RowNum() << " " << mat.ColNum() << std::endl;
				for (int i = 0; i < mat.RowNum(); i++)
				{
					for (int j = 0; j < mat.ColNum(); j++)
						file << mat(i, j) << " ";
					file << std::endl;
				}
				file.close();
			}
			else {
				std::cerr << "Error: could not create file " << inFileName << " for writing." << std::endl;
				return false;
			}

			return true;
		}
	};

	//////////////////////               Default Matrix typdefs         ////////////////((////
	typedef Matrix<int>     MatrixInt;
	typedef Matrix<float>   MatrixFlt;
	typedef Matrix<double>  MatrixDbl;
	typedef Matrix<Complex> MatrixComplex;

	typedef Matrix<int>     MatI;
	typedef Matrix<float>   MatF;
	typedef Matrix<double>  MatD;
	typedef Matrix<Complex> MatC;
}
#endif