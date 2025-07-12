#if !defined  MML_MATRIX_H
#define MML_MATRIX_H

#include "MMLBase.h"

#include "Vector.h"

namespace MML
{
	template<class Type>
	class MatrixView {
		Type* data;
		int rows, cols, stride;
	public:
		MatrixView(Type* data_, int rows_, int cols_, int stride_)
			: data(data_), rows(rows_), cols(cols_), stride(stride_) {
		}

		// Element access
		Type& operator()(int i, int j) { return data[i * stride + j]; }
		const Type& operator()(int i, int j) const { return data[i * stride + j]; }

		int RowNum() const { return rows; }
		int ColNum() const { return cols; }

		// Optional: iterators for flat traversal
		// ... (can be added as needed)
	};

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
		void ReleaseMemory()
		{
			if (_data != NULL) {
				delete[](_data[0]);
				delete[](_data);
			}
		}

		bool IsMatrixTypeComplex() const
		{
			return std::is_same_v<Type, std::complex<double>> ||
				std::is_same_v<Type, std::complex<float>> ||
				std::is_same_v<Type, std::complex<long double>>;
		}

	public:
		typedef Type value_type;      // make Type available externally

		///////////////////////          Constructors and destructor       //////////////////////
		explicit Matrix() : _rows(0), _cols(0), _data{ nullptr } {}
		explicit Matrix(int rows, int cols) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);
			if constexpr (is_MML_simple_numeric<Type>) {
				for (int i = 0; i < _rows; ++i)
					for (int j = 0; j < _cols; ++j)
						_data[i][j] = 0;
			}
		}
		explicit Matrix(int rows, int cols, const Type& val) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = val;
		}
		// Constructor from std::vector<std::vector<Type>>
		explicit Matrix(const std::vector<std::vector<Type>>& values) {
			if (values.empty() || values[0].empty()) {
				_rows = 0;
				_cols = 0;
				_data = nullptr;
				return;
			}
			_rows = static_cast<int>(values.size());
			_cols = static_cast<int>(values[0].size());

			Init(_rows, _cols);

			for (int i = 0; i < _rows; ++i) {
				if (static_cast<int>(values[i].size()) != _cols)
					throw MatrixDimensionError("Matrix::Matrix(std::vector<std::vector<Type>>): inconsistent row sizes", _rows, _cols, i, static_cast<int>(values[i].size()));
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = values[i][j];
			}
		}
		// Constructor from std::array<std::array<Type, Cols>, Rows>
		template <std::size_t Rows, std::size_t Cols>
		explicit Matrix(const std::array<std::array<Type, Cols>, Rows>& arr)
			: _rows(static_cast<int>(Rows)), _cols(static_cast<int>(Cols))
		{
			Init(_rows, _cols);
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = arr[i][j];
		}
		// Constructor from a pointer to continuous 2D array (can be in row-, or column-wise memory layout)
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
		Matrix(Matrix&& m) noexcept : _data(nullptr), _rows(0), _cols(0)
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
			ReleaseMemory();
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

				ReleaseMemory();

				Init(rows, cols);

				if constexpr (is_MML_simple_numeric<Type>) {
					for (int i = 0; i < _rows; ++i)
						for (int j = 0; j < _cols; ++j)
							_data[i][j] = 0;
				}
			}
		}
		void MakeUnitMatrix(void)
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

		// Static factory methods
		static Matrix GetUnitMatrix(int dim)
		{
			if (dim <= 0)
				throw MatrixDimensionError("Matrix::GetUnitMatrix - dimension must be positive", dim, dim, -1, -1);

			Matrix unitMat(dim, dim);
			unitMat.MakeUnitMatrix();

			return unitMat;
		}
		static Matrix GetDiagonalMatrix(const Vector<Type>& diagValues)
		{
			if (diagValues.size() <= 0)
				throw MatrixDimensionError("Matrix::GetDiagonalMatrix - vector size must be positive", diagValues.size(), diagValues.size(), -1, -1);
			Matrix diagMat(diagValues.size(), diagValues.size());
			for (int i = 0; i < diagValues.size(); ++i)
				diagMat._data[i][i] = diagValues[i];
			return diagMat;
		}

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

		///////////////////////              Creating views                //////////////////////
		// Block view (submatrix)
		MatrixView<Type> block(int startRow, int startCol, int numRows, int numCols) {
			if (startRow < 0 || startCol < 0 || numRows <= 0 || numCols <= 0 ||
				startRow + numRows > _rows || startCol + numCols > _cols)
				throw MatrixDimensionError("Matrix::block - invalid block parameters", _rows, _cols, startRow, startCol);
			return MatrixView<Type>(&_data[startRow][startCol], numRows, numCols, _cols);
		}
		const MatrixView<Type> block(int startRow, int startCol, int numRows, int numCols) const {
			if (startRow < 0 || startCol < 0 || numRows <= 0 || numCols <= 0 ||
				startRow + numRows > _rows || startCol + numCols > _cols)
				throw MatrixDimensionError("Matrix::block - invalid block parameters", _rows, _cols, startRow, startCol);
			return MatrixView<Type>(&_data[startRow][startCol], numRows, numCols, _cols);
		}

		// Row view
		MatrixView<Type> row_view(int rowIdx) {
			if (rowIdx < 0 || rowIdx >= _rows)
				throw MatrixAccessBoundsError("Matrix::row_view", rowIdx, 0, _rows, _cols);
			return MatrixView<Type>(_data[rowIdx], 1, _cols, _cols);
		}
		const MatrixView<Type> row_view(int rowIdx) const {
			if (rowIdx < 0 || rowIdx >= _rows)
				throw MatrixAccessBoundsError("Matrix::row_view", rowIdx, 0, _rows, _cols);
			return MatrixView<Type>(_data[rowIdx], 1, _cols, _cols);
		}

		// Column view
		MatrixView<Type> col_view(int colIdx) {
			if (colIdx < 0 || colIdx >= _cols)
				throw MatrixAccessBoundsError("Matrix::col_view", 0, colIdx, _rows, _cols);
			return MatrixView<Type>(&_data[0][colIdx], _rows, 1, _cols);
		}
		const MatrixView<Type> col_view(int colIdx) const {
			if (colIdx < 0 || colIdx >= _cols)
				throw MatrixAccessBoundsError("Matrix::col_view", 0, colIdx, _rows, _cols);
			return MatrixView<Type>(&_data[0][colIdx], _rows, 1, _cols);
		}
		
		///////////////////////              Standard stuff                //////////////////////
		inline int RowNum() const { return _rows; }
		inline int ColNum() const { return _cols; }

		Vector<Type> GetDiagonal() const
		{
			return VectorFromDiagonal();
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
		Matrix GetSubmatrix(int start_row_ind, int start_col_ind, int row_num, int col_num) const
		{
			if (start_row_ind < 0 || start_row_ind >= RowNum() || start_col_ind < 0 || start_col_ind >= ColNum())
				throw MatrixDimensionError("Matrix::GetSubmatrix - invalid starting row or column index", _rows, _cols, start_row_ind, start_col_ind);
			if (row_num <= 0 || col_num <= 0)
				throw MatrixDimensionError("Matrix::GetSubmatrix - rowNum and colNum must be positive", row_num, col_num, -1, -1);
			if (start_row_ind + row_num > RowNum() || start_col_ind + col_num > ColNum())
				throw MatrixDimensionError("Matrix::GetSubmatrix - submatrix out of bounds", _rows, _cols, start_row_ind, start_col_ind);

			return Matrix(*this, start_row_ind, start_col_ind, row_num, col_num);
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

		///////////////////////               Matrix properties            //////////////////////
		bool IsUnit(double eps = Defaults::IsMatrixUnitTolerance) const
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
		bool IsDiagonal(double eps = Defaults::IsMatrixDiagonalTolerance) const
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
				{
					if (IsMatrixTypeComplex())
					{
						// for Complex types, we check for conjugate symmetry
						if ((*this)[i][j] != -std::conj((*this)[j][i]))
							return false;
					}
					else
					{
						if ((*this)[i][j] != -(*this)[j][i])
							return false;
					}
				}
			return true;
		}

		///////////////////////             Matrix norm calculations       //////////////////////
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

		///////////////////////             Iterator support               //////////////////////
		// Flat iterator for Matrix (row-major order)
		class iterator {
			Type* ptr;
		public:
			using iterator_category = std::random_access_iterator_tag;
			using value_type = Type;
			using difference_type = std::ptrdiff_t;
			using pointer = Type*;
			using reference = Type&;

			iterator(Type* p) : ptr(p) {}
			iterator& operator++() { ++ptr; return *this; }
			iterator operator++(int) { iterator tmp = *this; ++ptr; return tmp; }
			iterator& operator--() { --ptr; return *this; }
			iterator operator--(int) { iterator tmp = *this; --ptr; return tmp; }
			iterator operator+(difference_type n) const { return iterator(ptr + n); }
			iterator operator-(difference_type n) const { return iterator(ptr - n); }
			difference_type operator-(const iterator& other) const { return ptr - other.ptr; }
			bool operator==(const iterator& other) const { return ptr == other.ptr; }
			bool operator!=(const iterator& other) const { return ptr != other.ptr; }
			reference operator*() const { return *ptr; }
			pointer operator->() const { return ptr; }
		};

		class const_iterator {
			const Type* ptr;
		public:
			using iterator_category = std::random_access_iterator_tag;
			using value_type = Type;
			using difference_type = std::ptrdiff_t;
			using pointer = const Type*;
			using reference = const Type&;

			const_iterator(const Type* p) : ptr(p) {}
			const_iterator& operator++() { ++ptr; return *this; }
			const_iterator operator++(int) { const_iterator tmp = *this; ++ptr; return tmp; }
			const_iterator& operator--() { --ptr; return *this; }
			const_iterator operator--(int) { const_iterator tmp = *this; --ptr; return tmp; }
			const_iterator operator+(difference_type n) const { return const_iterator(ptr + n); }
			const_iterator operator-(difference_type n) const { return const_iterator(ptr - n); }
			difference_type operator-(const const_iterator& other) const { return ptr - other.ptr; }
			bool operator==(const const_iterator& other) const { return ptr == other.ptr; }
			bool operator!=(const const_iterator& other) const { return ptr != other.ptr; }
			reference operator*() const { return *ptr; }
			pointer operator->() const { return ptr; }
		};

		iterator begin() { return iterator(_rows && _cols ? &_data[0][0] : nullptr); }
		iterator end() { return iterator(_rows && _cols ? &_data[0][0] + _rows * _cols : nullptr); }
		const_iterator begin() const { return const_iterator(_rows && _cols ? &_data[0][0] : nullptr); }
		const_iterator end()   const { return const_iterator(_rows && _cols ? &_data[0][0] + _rows * _cols : nullptr); }
		const_iterator cbegin() const { return begin(); }
		const_iterator cend()   const { return end(); }

		// Row iterator: iterates over elements in a single row
		class row_iterator {
			Type* ptr;
		public:
			using iterator_category = std::random_access_iterator_tag;
			using value_type = Type;
			using difference_type = std::ptrdiff_t;
			using pointer = Type*;
			using reference = Type&;

			row_iterator(Type* p) : ptr(p) {}
			row_iterator& operator++() { ++ptr; return *this; }
			row_iterator operator++(int) { row_iterator tmp = *this; ++ptr; return tmp; }
			row_iterator& operator--() { --ptr; return *this; }
			row_iterator operator--(int) { row_iterator tmp = *this; --ptr; return tmp; }
			row_iterator operator+(difference_type n) const { return row_iterator(ptr + n); }
			row_iterator operator-(difference_type n) const { return row_iterator(ptr - n); }
			difference_type operator-(const row_iterator& other) const { return ptr - other.ptr; }
			bool operator==(const row_iterator& other) const { return ptr == other.ptr; }
			bool operator!=(const row_iterator& other) const { return ptr != other.ptr; }
			reference operator*() const { return *ptr; }
			pointer operator->() const { return ptr; }
		};

		// Column iterator: iterates over elements in a single column
		class col_iterator {
			Type* ptr;
			int stride; // number of columns in the matrix
		public:
			using iterator_category = std::random_access_iterator_tag;
			using value_type = Type;
			using difference_type = std::ptrdiff_t;
			using pointer = Type*;
			using reference = Type&;

			col_iterator(Type* p, int stride_) : ptr(p), stride(stride_) {}
			col_iterator& operator++() { ptr += stride; return *this; }
			col_iterator operator++(int) { col_iterator tmp = *this; ptr += stride; return tmp; }
			col_iterator& operator--() { ptr -= stride; return *this; }
			col_iterator operator--(int) { col_iterator tmp = *this; ptr -= stride; return tmp; }
			bool operator==(const col_iterator& other) const { return ptr == other.ptr; }
			bool operator!=(const col_iterator& other) const { return ptr != other.ptr; }
			reference operator*() const { return *ptr; }
			pointer operator->() const { return ptr; }
		};

		// Proxy for a row
		class row_range {
			Type* row_ptr;
			int cols;
		public:
			row_range(Type* p, int n) : row_ptr(p), cols(n) {}
			row_iterator begin() { return row_iterator(row_ptr); }
			row_iterator end() { return row_iterator(row_ptr + cols); }
		};

		// Proxy for a column
		class col_range {
			Type* col_ptr;
			int rows, stride;
		public:
			col_range(Type* p, int n, int stride_) : col_ptr(p), rows(n), stride(stride_) {}
			col_iterator begin() { return col_iterator(col_ptr, stride); }
			col_iterator end() { return col_iterator(col_ptr + rows * stride, stride); }
		};

		row_range row(int rowInd) {
			if (rowInd < 0 || rowInd >= _rows) throw MatrixAccessBoundsError("Matrix::row", rowInd, 0, _rows, _cols);
			return row_range(_data[rowInd], _cols);
		}
		col_range col(int colInd) {
			if (colInd < 0 || colInd >= _cols) throw MatrixAccessBoundsError("Matrix::col", 0, colInd, _rows, _cols);
			return col_range(&_data[0][colInd], _rows, _cols);
		}

		// Const row and column accessors
		const row_range row(int rowInd) const {
			if (rowInd < 0 || rowInd >= _rows) throw MatrixAccessBoundsError("Matrix::row", rowInd, 0, _rows, _cols);
			return row_range(_data[rowInd], _cols);
		}
		const col_range col(int colInd) const {
			if (colInd < 0 || colInd >= _cols) throw MatrixAccessBoundsError("Matrix::col", 0, colInd, _rows, _cols);
			return col_range(&_data[0][colInd], _rows, _cols);
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

		bool IsEqualTo(const Matrix<Type>& b, Type eps = Defaults::MatrixIsEqualTolerance) const
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
		static bool AreEqual(const Matrix& a, const Matrix& b, Type eps = Defaults::MatrixIsEqualTolerance)
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

		// Load a matrix from a CSV file (comma-separated values)
		static bool LoadFromCSV(const std::string& filename, Matrix& outMat) {
			std::ifstream file(filename);
			if (!file.is_open()) {
				std::cerr << "Error: could not open file " << filename << " for reading." << std::endl;
				return false;
			}
			std::vector<std::vector<Type>> data;
			std::string line;
			while (std::getline(file, line)) {
				std::vector<Type> row;
				std::stringstream ss(line);
				std::string cell;
				while (std::getline(ss, cell, ',')) {
					std::stringstream cellStream(cell);
					Type value;
					cellStream >> value;
					row.push_back(value);
				}
				if (!row.empty())
					data.push_back(row);
			}
			file.close();
			if (data.empty()) return false;
			int rows = static_cast<int>(data.size());
			int cols = static_cast<int>(data[0].size());
			outMat.Resize(rows, cols);
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					outMat[i][j] = data[i][j];
			return true;
		}

		// Save a matrix to a CSV file (comma-separated values)
		static bool SaveToCSV(const Matrix& mat, const std::string& filename) {
			std::ofstream file(filename);
			if (!file.is_open()) {
				std::cerr << "Error: could not create file " << filename << " for writing." << std::endl;
				return false;
			}
			for (int i = 0; i < mat.RowNum(); ++i) {
				for (int j = 0; j < mat.ColNum(); ++j) {
					file << mat(i, j);
					if (j < mat.ColNum() - 1)
						file << ",";
				}
				file << "\n";
			}
			file.close();
			return true;
		}

		// Save a matrix to a simple JSON array of arrays
		static bool SaveToJSON(const Matrix& mat, const std::string& filename) {
			std::ofstream file(filename);
			if (!file.is_open()) {
				std::cerr << "Error: could not create file " << filename << " for writing." << std::endl;
				return false;
			}
			file << "[";
			for (int i = 0; i < mat.RowNum(); ++i) {
				file << "[";
				for (int j = 0; j < mat.ColNum(); ++j) {
					file << mat(i, j);
					if (j < mat.ColNum() - 1)
						file << ",";
				}
				file << "]";
				if (i < mat.RowNum() - 1)
					file << ",";
			}
			file << "]\n";
			file.close();
			return true;
		}

		// Load a matrix from a simple JSON array of arrays (no error checking, expects valid format)
		static bool LoadFromJSON(const std::string& filename, Matrix& outMat) {
			std::ifstream file(filename);
			if (!file.is_open()) {
				std::cerr << "Error: could not open file " << filename << " for reading." << std::endl;
				return false;
			}
			std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
			file.close();
			std::vector<std::vector<Type>> data;
			size_t pos = 0;
			while ((pos = content.find('[')) != std::string::npos) {
				size_t end = content.find(']', pos);
				if (end == std::string::npos) break;
				std::string rowStr = content.substr(pos + 1, end - pos - 1);
				std::vector<Type> row;
				std::stringstream ss(rowStr);
				std::string cell;
				while (std::getline(ss, cell, ',')) {
					std::stringstream cellStream(cell);
					Type value;
					cellStream >> value;
					row.push_back(value);
				}
				if (!row.empty())
					data.push_back(row);
				content = content.substr(end + 1);
			}
			if (data.empty()) return false;
			int rows = static_cast<int>(data.size());
			int cols = static_cast<int>(data[0].size());
			outMat.Resize(rows, cols);
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					outMat[i][j] = data[i][j];
			return true;
		}

		// Save a matrix to a MATLAB .m file (as a variable assignment)
		static bool SaveToMatlab(const Matrix& mat, const std::string& filename, const std::string& varName = "A") {
			std::ofstream file(filename);
			if (!file.is_open()) {
				std::cerr << "Error: could not create file " << filename << " for writing." << std::endl;
				return false;
			}
			file << varName << " = [";
			for (int i = 0; i < mat.RowNum(); ++i) {
				for (int j = 0; j < mat.ColNum(); ++j) {
					file << mat(i, j);
					if (j < mat.ColNum() - 1)
						file << " ";
				}
				if (i < mat.RowNum() - 1)
					file << "; ";
			}
			file << "];\n";
			file.close();
			return true;
		}

		// Load a matrix from a MATLAB .m file (expects a single variable assignment, e.g., A = [1 2; 3 4];)
		// Only supports numeric types, no error checking for malformed files
		static bool LoadFromMatlab(const std::string& filename, Matrix& outMat) {
			std::ifstream file(filename);
			if (!file.is_open()) {
				std::cerr << "Error: could not open file " << filename << " for reading." << std::endl;
				return false;
			}
			std::string line, content;
			while (std::getline(file, line)) {
				auto pos = line.find('=');
				if (pos != std::string::npos) {
					content = line.substr(pos + 1);
					break;
				}
			}
			file.close();
			if (content.empty()) return false;
			// Remove brackets and semicolons
			content.erase(std::remove(content.begin(), content.end(), '['), content.end());
			content.erase(std::remove(content.begin(), content.end(), ']'), content.end());
			content.erase(std::remove(content.begin(), content.end(), ';'), content.end());
			std::stringstream ss(content);
			std::vector<std::vector<Type>> data;
			std::string rowStr;
			while (std::getline(ss, rowStr, ';')) {
				std::stringstream rowStream(rowStr);
				std::vector<Type> row;
				Type value;
				while (rowStream >> value) {
					row.push_back(value);
				}
				if (!row.empty())
					data.push_back(row);
			}
			if (data.empty()) return false;
			int rows = static_cast<int>(data.size());
			int cols = static_cast<int>(data[0].size());
			outMat.Resize(rows, cols);
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					outMat[i][j] = data[i][j];
			return true;
		}

		// Serialize matrix to a binary file using the portable NPY (NumPy) format
		static bool SaveToNpy(const Matrix& mat, const std::string& filename) {
			std::ofstream file(filename, std::ios::binary);
			if (!file.is_open()) {
				std::cerr << "Error: could not create file " << filename << " for writing." << std::endl;
				return false;
			}
			// NPY header
			std::string magic = "\x93NUMPY";
			uint8_t major = 1, minor = 0;
			std::string dtype;
			if constexpr (std::is_same_v<Type, double>) dtype = "<f8";
			else if constexpr (std::is_same_v<Type, float>) dtype = "<f4";
			else if constexpr (std::is_same_v<Type, int>) dtype = "<i4";
			else {
				std::cerr << "Error: unsupported type for NPY serialization." << std::endl;
				return false;
			}
			std::ostringstream header;
			header << "{'descr': '" << dtype << "', 'fortran_order': False, 'shape': (" << mat.RowNum() << ", " << mat.ColNum() << "), }";
			std::string headerStr = header.str();
			// Pad header to 16-byte alignment
			size_t header_len = 10 + headerStr.size();
			size_t pad = 16 - (header_len % 16);
			headerStr.append(pad, ' ');
			headerStr.back() = '\n';
			uint16_t header_size = static_cast<uint16_t>(headerStr.size());
			// Write header
			file.write(magic.c_str(), 6);
			file.put(major);
			file.put(minor);
			file.write(reinterpret_cast<const char*>(&header_size), 2);
			file.write(headerStr.c_str(), headerStr.size());
			// Write data (row-major)
			for (int i = 0; i < mat.RowNum(); ++i)
				file.write(reinterpret_cast<const char*>(mat._data[i]), sizeof(Type) * mat.ColNum());
			file.close();
			return true;
		}

		// Load matrix from a binary NPY (NumPy) file
		static bool LoadFromNpy(const std::string& filename, Matrix& outMat) {
			std::ifstream file(filename, std::ios::binary);
			if (!file.is_open()) {
				std::cerr << "Error: could not open file " << filename << " for reading." << std::endl;
				return false;
			}
			char magic[6];
			file.read(magic, 6);
			if (std::string(magic, 6) != "\x93NUMPY") return false;
			uint8_t major = file.get();
			uint8_t minor = file.get();
			uint16_t header_size;
			file.read(reinterpret_cast<char*>(&header_size), 2);
			std::string header(header_size, ' ');
			file.read(&header[0], header_size);
			// Parse shape
			size_t pos1 = header.find('(');
			size_t pos2 = header.find(')', pos1);
			if (pos1 == std::string::npos || pos2 == std::string::npos) return false;
			std::string shape = header.substr(pos1 + 1, pos2 - pos1 - 1);
			int rows = 0, cols = 0;
			std::replace(shape.begin(), shape.end(), ',', ' ');
			std::istringstream iss(shape);
			iss >> rows >> cols;
			if (rows <= 0 || cols <= 0) return false;
			outMat.Resize(rows, cols);
			// Read data
			for (int i = 0; i < rows; ++i)
				file.read(reinterpret_cast<char*>(outMat._data[i]), sizeof(Type) * cols);
			file.close();
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