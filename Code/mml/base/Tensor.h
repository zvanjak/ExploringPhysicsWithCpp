#if !defined MML_TENSORS_H
#define MML_TENSORS_H

#include "MMLBase.h"

#include "interfaces/ITensor.h"

#include "base/VectorN.h"
#include "base/MatrixNM.h"

namespace MML
{
	template <int N>
	class Tensor2 : public ITensor2<N>
	{
		MatrixNM<Real, N, N> _coeff;
	public:
		int _numContravar = 0;
		int _numCovar = 0;
		bool _isContravar[2];

		Tensor2(int nCovar, int nContraVar) : _numCovar(nCovar), _numContravar(nContraVar)
		{
			if ( _numContravar < 0 || _numCovar  < 0 || _numContravar + _numCovar != 2 )
				throw TensorCovarContravarNumError("Tensor2 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;
		}
		Tensor2(int nCovar, int nContraVar, std::initializer_list<Real> values) : _numCovar(nCovar), _numContravar(nContraVar)
		{
			if ( _numContravar < 0 || _numCovar  < 0 || _numContravar + _numCovar != 2 )
				throw TensorCovarContravarNumError("Tensor2 ctor, wrong number of covariant and contravariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;

			auto val = values.begin();
			for (size_t i = 0; i < N; ++i)
				for (size_t j = 0; j < N; ++j)
					if (val != values.end())
						_coeff[i][j] = *val++;
					else
						_coeff[i][j] = 0.0;
		}

		int  NumContravar() const override { return _numContravar; }
		int  NumCovar()     const override { return _numCovar; }

		bool IsContravar(int i) const { return _isContravar[i]; }
		bool IsCovar(int i) const			{ return !_isContravar[i]; }

		Real  operator()(int i, int j) const override { return _coeff[i][j]; }
		Real& operator()(int i, int j) override				{ return _coeff[i][j]; }

		Tensor2 operator+(const Tensor2& other) const
		{
			if (_numContravar != other._numContravar || _numCovar != other._numCovar)
				throw TensorCovarContravarArithmeticError("Tensor2 operator+, wrong number of contravariant and covariant indices", _numContravar, _numCovar, other._numContravar, other._numCovar);

			Tensor2 result(_numContravar, _numCovar);

			result._coeff = _coeff + other._coeff;

			return result;
		}
		Tensor2 operator-(const Tensor2& other) const
		{
			if (_numContravar != other._numContravar || _numCovar != other._numCovar)
				throw TensorCovarContravarArithmeticError("Tensor2 operator-, wrong number of contravariant and covariant indices", _numContravar, _numCovar, other._numContravar, other._numCovar);

			Tensor2 result(_numContravar, _numCovar);

			result._coeff = _coeff - other._coeff;

			return result;
		}
		Tensor2 operator*=(Real scalar) const
		{
			Tensor2 result(_numContravar, _numCovar);

			result._coeff = _coeff * scalar;

			return result;
		}
		Tensor2 operator/(Real scalar) const
		{
			Tensor2 result(_numContravar, _numCovar);

			result._coeff = _coeff / scalar;

			return result;
		}

		friend Tensor2 operator*(Real scalar, const Tensor2& b)
		{
			Tensor2 result(b.NumContravar(), b.NumCovar());

			result._coeff = b._coeff * scalar;

			return result;
		}

		Real Contract() const
		{
			if (_numContravar != 1 || _numCovar != 1)
				throw TensorCovarContravarNumError("Tensor2 Contract, wrong number of contravariant and covariant indices", _numContravar, _numCovar);

			Real result = 0.0;
			for (int i = 0; i < N; i++)
				result += _coeff[i][i];

			return result;
		}
		Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2) const
		{
			Real sum = 0.0;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					sum += _coeff[i][j] * v1[i] * v2[j];

			return sum;
		}

		void   Print(std::ostream& stream, int width, int precision) const
		{
			stream << std::fixed << "N = " << N << std::endl;

			for (size_t i = 0; i < N; i++)
			{
				stream << "[ ";
				for (size_t j = 0; j < N; j++)
					stream << std::setw(width) << std::setprecision(precision) << _coeff[i][j] << ", ";
				stream << " ]" << std::endl;
			}
		}
		friend std::ostream& operator<<(std::ostream& stream, const Tensor2& a)
		{
			a.Print(stream, 15, 10);

			return stream;
		}
	};

	template <int N>
	class Tensor3 : public ITensor3<N>
	{
		Real _coeff[N][N][N] = { 0 };
	public:
		int _numContravar;
		int _numCovar;
		bool _isContravar[3];

		Tensor3(int nCovar, int nContraVar) : _numCovar(nCovar), _numContravar(nContraVar)
		{
			if (_numContravar + _numCovar != 3)
				throw TensorCovarContravarNumError("Tensor3 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;
		}

		int   NumContravar() const override { return _numContravar; }
		int   NumCovar()     const override { return _numCovar; }

		Real  operator()(int i, int j, int k) const override { return _coeff[i][j][k]; }
		Real& operator()(int i, int j, int k) override { return _coeff[i][j][k]; }

		Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2, const VectorN<Real, N>& v3) const
		{
			Real sum = 0.0;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						sum += _coeff[i][j][k] * v1[i] * v2[j] * v3[k];

			return sum;
		}

		VectorN<Real, N> Contract(int ind1, int ind2) const
		{
			VectorN<Real, N> result;

			if (ind1 < 0 || ind1 >= N || ind2 < 0 || ind2 >= N)
				throw TensorIndexError("Tensor3 Contract, wrong index");

			if (ind1 == ind2)
				throw TensorIndexError("Tensor3 Contract, indices are the same");

			if (ind1 == 0 && ind2 == 1)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						result[i] += _coeff[j][j][i];
			}
			else if (ind1 == 0 && ind2 == 2)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						result[i] += _coeff[j][i][j];
			}
			else if (ind1 == 1 && ind2 == 2)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						result[i] += _coeff[i][j][j];
			}
			else
				throw TensorIndexError("Tensor3 Contract, wrong indices");

			return result;
		}
	};

	template <int N>
	class Tensor4 : public ITensor4<N>
	{
		Real _coeff[N][N][N][N] = { 0 };
	public:
		int _numContravar;
		int _numCovar;
		bool _isContravar[4];

		Tensor4(int nCovar, int nContraVar) : _numCovar(nCovar), _numContravar(nContraVar) 
		{
			if (_numContravar + _numCovar != 4)
				throw TensorCovarContravarNumError("Tensor4 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;
		}

		int   NumContravar() const override { return _numContravar; }
		int   NumCovar()     const override { return _numCovar; }

		Real  operator()(int i, int j, int k, int l) const override { return _coeff[i][j][k][l]; }
		Real& operator()(int i, int j, int k, int l) override { return _coeff[i][j][k][l]; }

		Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2, const VectorN<Real, N>& v3, const VectorN<Real, N>& v4) const
		{
			Real sum = 0.0;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
							sum += _coeff[i][j][k][l] * v1[i] * v2[j] * v3[k] * v4[l];

			return sum;
		}

		Tensor2<N> Contract(int ind1, int ind2) const
		{
			MatrixNM<Real, N, N> result;
			Tensor2<N> ret;

			// check covar contravar!!!!
			if (ind1 < 0 || ind1 >= N || ind2 < 0 || ind2 >= N)
				throw TensorIndexError("Tensor4 Contract, wrong index");

			if (ind1 == ind2)
				throw TensorIndexError("Tensor4 Contract, indices are the same");

			if (ind1 == 0 && ind2 == 1)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						for(int k = 0; k < N; k++)
							result[i][j] += _coeff[k][k][i][j];
			}
			else if (ind1 == 0 && ind2 == 2)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						for(int k = 0; k < N; k++)
							result[i][j] += _coeff[k][i][k][j];
			}
			else if (ind1 == 0 && ind2 == 3)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						for(int k = 0; k < N; k++)
							result[i][j] += _coeff[k][i][j][k];
			}
			else if (ind1 == 1 && ind2 == 2)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						for(int k = 0; k < N; k++)
							result[i][j] += _coeff[i][k][k][j];
			}
			else if (ind1 == 1 && ind2 == 3)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						for(int k = 0; k < N; k++)
							result[i][j] += _coeff[i][k][j][k];
			}
			else if (ind1 == 2 && ind2 == 3)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						for(int k = 0; k < N; k++)
							result[i][j] += _coeff[i][j][k][k];
			}
			else
				throw TensorIndexError("Tensor4 Contract, wrong indices");

			ret._coeff = result;

			return ret;
		}
	};

	template <int N>
	class Tensor5 : public ITensor5<N>
	{
		Real _coeff[N][N][N][N][N];
	public:
		int _numContravar;
		int _numCovar;
		bool _isContravar[5];

		Tensor5(int nCovar, int nContraVar) : _numCovar(nCovar), _numContravar(nContraVar) 
		{
			if (_numContravar + _numCovar != 5)
				throw TensorCovarContravarNumError("Tensor5 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;

		}

		int   NumContravar() const override { return _numContravar; }
		int   NumCovar()     const override { return _numCovar; }

		Real  operator()(int i, int j, int k, int l, int m) const override { return _coeff[i][j][k][l][m]; }
		Real& operator()(int i, int j, int k, int l, int m) override { return _coeff[i][j][k][l][m]; }
	};
}
#endif