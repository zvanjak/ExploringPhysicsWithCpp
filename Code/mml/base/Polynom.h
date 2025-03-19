#if !defined  MML_POLYNOM_H
#define MML_POLYNOM_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/MatrixNM.h"

namespace MML
{
	template <typename _Field>
	class Polynom
	{
	protected:
		// polynom coefficients - with _vecCoef[0] being the constant term
		std::vector<Real> _vecCoef;
	public:
		Polynom() {}
		Polynom(int n) { _vecCoef.resize(n + 1); }
		Polynom(const std::vector<Real>& vecCoef) : _vecCoef(vecCoef) {}
		Polynom(std::initializer_list<Real> list) : _vecCoef(list) {}
		Polynom(const Polynom& Copy) = default;

		Polynom& operator=(const Polynom& Copy) { _vecCoef = Copy._vecCoef; return *this; }

		int  GetDegree() const { return (int)_vecCoef.size() - 1; }
		void SetDegree(int newDeg) { _vecCoef.resize(newDeg + 1); }
		bool IsNullPolynom() const { return _vecCoef.size() == 0; }

		Real  operator[] (int i) const { return _vecCoef[i]; }
		Real& operator[] (int i) { return _vecCoef[i]; }

		_Field operator() (const _Field& x) const {
			int j = GetDegree();
			_Field p = _vecCoef[j] * _Field(1);

			while (j > 0)
				p = p * x + _vecCoef[--j];
			return p;
		}
		// Given the coefficients of a polynomial of degree nc as an array c[0..nc] of size nc+1 (with
		// c[0] being the constant term), and given a value x, this routine fills an output array pd of size
		// nd+1 with the value of the polynomial evaluated at x in pd[0], and the first nd derivatives at
		// x in pd[1..nd].
		void Derive(const Real x, Vector<Real>& pd)
		{
			int  nnd, j, i;
			int  nc = GetDegree();
			int  nd = pd.size() - 1;
			Real cnst = 1.0;

			pd[0] = (*this)[nc];
			for (j = 1; j < nd + 1; j++)
				pd[j] = 0.0;

			for (i = nc - 1; i >= 0; i--)
			{
				nnd = (nd < (nc - i) ? nd : nc - i);
				for (j = nnd; j > 0; j--)
					pd[j] = pd[j] * x + pd[j - 1];

				pd[0] = pd[0] * x + (*this)[i];
			}
			for (i = 2; i < nd + 1; i++) {
				cnst *= i;
				pd[i] *= cnst;
			}
		}

		bool IsEqual(const Polynom& b) const
		{
			if (_vecCoef.size() != b._vecCoef.size())
				return false;
			for (int i = 0; i < _vecCoef.size(); i++)
				if (_vecCoef[i] != b._vecCoef[i])
					return false;
			return true;
		}

		bool operator==(const Polynom& b) const
		{
			// TODO - nije nuzno - mogu biti isti, a razlicitog reda?
			// ako ne reduciramo red kad su highest coef 0!
			if (_vecCoef.size() != b._vecCoef.size())
				return false;
			for (int i = 0; i < _vecCoef.size(); i++)
				if (_vecCoef[i] != b._vecCoef[i])
					return false;
			return true;
		}

		Polynom operator+(const Polynom& b) const
		{
			Polynom result;
			int n = (int)std::max(_vecCoef.size(), b._vecCoef.size());
			result._vecCoef.resize(n);
			for (int i = 0; i < n; i++)
			{
				if (i < _vecCoef.size())
					result._vecCoef[i] += _vecCoef[i];
				if (i < b._vecCoef.size())
					result._vecCoef[i] += b._vecCoef[i];
			}
			return result;
		}

		Polynom operator-(const Polynom& b) const
		{
			Polynom result;
			int n = (int)std::max(_vecCoef.size(), b._vecCoef.size());
			result._vecCoef.resize(n);
			for (int i = 0; i < n; i++)
			{
				if (i < _vecCoef.size())
					result._vecCoef[i] += _vecCoef[i];
				if (i < b._vecCoef.size())
					result._vecCoef[i] -= b._vecCoef[i];
			}
			return result;
		}

		Polynom operator*(const Polynom& b) const
		{
			Polynom result;

			int n = (int)(_vecCoef.size() + b._vecCoef.size() - 1);
			result._vecCoef.resize(n);
			for (int i = 0; i < _vecCoef.size(); i++)
				for (int j = 0; j < b._vecCoef.size(); j++)
					result._vecCoef[i + j] += _vecCoef[i] * b._vecCoef[j];
			return result;
		}

		static void poldiv(const Polynom& u, const Polynom& v, Polynom& qout, Polynom& rout)
		{
			int k, j, n = u.GetDegree(), nv = v.GetDegree();

			// find real degree of v
			while (nv >= 0 && v._vecCoef[nv] == 0.)
				nv--;

			if (nv < 0)
				throw std::domain_error("poldiv divide by zero polynomial");

			Polynom r = u;
			Polynom q(u.GetDegree());
			for (k = n - nv; k >= 0; k--)
			{
				q[k] = r[nv + k] / v[nv];

				for (j = nv + k - 1; j >= k; j--)
					r[j] -= q[k] * v[j - k];
			}
			for (j = nv; j <= n; j++)
				r[j] = 0.0;


			int nq = q.GetDegree();
			while (nq >= 0 && q[nq] == 0.)
				nq--;

			// setting exact size for quotient
			qout.SetDegree(nq);
			for (j = 0; j <= nq; j++)
				qout[j] = q[j];

			// setting exact size for remainder
			rout.SetDegree(nv - 1);
			for (j = 0; j < nv; j++)
				rout[j] = r[j];
		}

		friend Polynom operator*(const Polynom& a, Real b)
		{
			Polynom ret;
			ret._vecCoef.resize(a.GetDegree() + 1);
			for (int i = 0; i <= a.GetDegree(); i++)
				ret._vecCoef[i] = a._vecCoef[i] * b;
			return ret;
		}

		friend Polynom operator*(Real a, const Polynom& b)
		{
			Polynom ret;
			ret._vecCoef.resize(b.GetDegree() + 1);
			for (int i = 0; i <= b.GetDegree(); i++)
				ret._vecCoef[i] = a * b._vecCoef[i];
			return ret;
		}

		friend Polynom operator/(const Polynom& a, Real b)
		{
			Polynom ret;
			ret._vecCoef.resize(a.GetDegree() + 1);
			for (int i = 0; i <= a.GetDegree(); i++)
				ret._vecCoef[i] = a._vecCoef[i] / b;
			return ret;
		}

		///////////////////////////               I/O                 ///////////////////////////
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}

		std::ostream& Print(std::ostream& stream, int width, int precision) const
		{
			// first, save current stream state
			std::ios_base::fmtflags f(stream.flags());

			// change formatting
			stream << std::fixed << std::setw(width) << std::setprecision(precision);

			// print the polynomial
			Print(stream);

			// restore stream state
			stream.flags(f);

			return stream;
		}

		std::ostream& Print(std::ostream& stream) const
		{
			for (int i = (int)_vecCoef.size() - 1; i >= 0; i--)
			{
				if (std::abs(_vecCoef[i]) == 0.0)
					continue;

				// handling first term
				if (i == _vecCoef.size() - 1)
				{
					if (i == 0) // means we have only constant term
					{
						stream << _vecCoef[i];
						return stream;
					}

					if (_vecCoef[i] < 0.0)
					{
						if (_vecCoef[i] == -1.0)
							stream << "-x";
						else
							stream << _vecCoef[i] << "*x";;
					}
					else
					{
						if (_vecCoef[i] == 1.0)
							stream << "x";
						else
							stream << _vecCoef[i] << "*x";;
					}
					// adding x exponent
					if (i > 1)
						stream << "^" << i;
				}
				else // handling other terms
				{
					if (i == 0)
					{
						if (_vecCoef[i] > 0.0)
							stream << " + " << _vecCoef[i];
						else
							stream << " - " << std::abs(_vecCoef[i]);
						return stream;
					}

					if (_vecCoef[i] > 0.0)
					{
						if (_vecCoef[i] == 1.0)
							stream << " + x";
						else
							stream << " + " << _vecCoef[i] << "*x";;
					}
					else
					{
						if (_vecCoef[i] == -1.0)
							stream << " - x";
						else
							stream << " - " << std::abs(_vecCoef[i]) << "*x";;
					}

					if (i > 1)
						stream << "^" << i;
				}
			}

			return stream;
		}

		friend std::ostream& operator<<(std::ostream& stream, const Polynom<_Field>& a)
		{
			a.Print(stream);

			return stream;
		}
	};

	class RealPolynom : public Polynom<Real>, public IRealFunction
	{
	public:
		RealPolynom() {}
		RealPolynom(int n) { _vecCoef.resize(n + 1); }
		RealPolynom(const std::vector<Real>& vecCoef) : Polynom(vecCoef) {}
		RealPolynom(std::initializer_list<Real> list) : Polynom(list) {}
		RealPolynom(const Polynom& Copy) : Polynom(Copy) {}
		~RealPolynom() {}

		Real operator()(Real x) const { return Polynom::operator()(x); }
	};

	//typedef Polynom<Real>         RealPolynom;
	typedef Polynom<Complex>   ComplexPolynom;

	typedef Polynom<MatrixNM<Real, 2, 2>>       Matrix2Polynom;
	typedef Polynom<MatrixNM<Real, 3, 3>>       Matrix3Polynom;
	typedef Polynom<MatrixNM<Real, 4, 4>>       Matrix4Polynom;
}

#endif