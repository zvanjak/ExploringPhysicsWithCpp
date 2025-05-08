#if !defined MML_INTERPOLATEDFUNCTION_H
#define MML_INTERPOLATEDFUNCTION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "base/Function.h"

namespace MML
{
	class RealFunctionInterpolated : public IRealFunction
	{
	private:
		int		_numPoints, _usedPoints;
		Vector<Real> _x, _y;						// we are storing copies of given values!

	public:
		RealFunctionInterpolated(const Vector<Real> &x, const Vector<Real> &y, 
														 int usedPointsInInterpolation)
							:  _x(x), _y(y), _numPoints(x.size()), _usedPoints(usedPointsInInterpolation)
		{
			// throw if not enough points
			if (_numPoints < 2 || _usedPoints < 2 || _usedPoints > _numPoints)
				throw RealFuncInterpInitError("RealFunctionInterpolated size error");
		}

		virtual ~RealFunctionInterpolated() {}

		Real virtual calcInterpValue(int startInd, Real x) const = 0;

		inline Real MinX() const { return X(0); }
		inline Real MaxX() const { return X(_numPoints-1); }

		inline Real X(int i) const { return _x[i]; }
		inline Real Y(int i) const { return _y[i]; }

		inline int	getNumPoints() const { return _numPoints; }
		inline int  getInterpOrder() const { return _usedPoints; }

		Real operator()(Real x) const
		{
			int startInd = locate(x);
			return calcInterpValue(startInd, x);
		}

		// Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
		// xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
		// increasing or decreasing. The returned value is not less than 0, nor greater than _numPoints-1.
		int locate(const Real x) const
		{
			int  indUpper, indMid, indLower;
			bool isAscending = (_x[_numPoints - 1] >= _x[0]);

			indLower = 0;
			indUpper = _numPoints - 1;

			while (indUpper - indLower > 1) {
				indMid = (indUpper + indLower) / 2;
				if (x >= _x[indMid] == isAscending)
					indLower = indMid;
				else
					indUpper = indMid;
			}
			return std::max(0, std::min(_numPoints - _usedPoints, indLower - ((_usedPoints - 2) / 2)));
		}
	};

	////////////////////////             LINEAR INTERPOLATION						  ///////////////////////
	class LinearInterpRealFunc : public RealFunctionInterpolated
	{
	public:
		LinearInterpRealFunc(const Vector<Real>& xv, Vector<Real>& yv) 
						: RealFunctionInterpolated(xv, yv, 2) {}

		Real calcInterpValue(int j, Real x) const {
			if (X(j) == X(j + 1)) 
				return Y(j);
			else 
				return Y(j) + ((x - X(j)) / (X(j + 1) - X(j)) * (Y(j + 1) - Y(j)));
		}
	};

	////////////////////////           POLYNOMIAL INTERPOLATION						///////////////////////
	// Polynomial interpolation object.Construct with x and y vectors, and the number M of points
	// to be used locally(polynomial order plus one), then call interp for interpolated values.
	class PolynomInterpRealFunc : public RealFunctionInterpolated
	{
	private:
		mutable Real _errorEst;
	public:
		PolynomInterpRealFunc(const Vector<Real>& xv, Vector<Real>& yv, int m)
						: RealFunctionInterpolated(xv, yv, m), _errorEst(0.) { }

		Real getLastErrorEst() const { return _errorEst; }	

		// Given a value x, and using pointers to data xx and yy, this routine returns an interpolated
		// value y, and stores an error estimate _errorEst.The returned value is obtained by mm-point polynomial
		// interpolation on the subrange xx[startInd..startInd + mm - 1].
		Real calcInterpValue(int startInd, Real x) const
		{
			int i, m, ns = 0;
			Real y, den, dif, dift, ho, hp, w, dy;
			Vector<Real> c(getInterpOrder()), d(getInterpOrder());

			dif = std::abs(x - X(startInd));

			// First we find the index ns of the closest table entry
			for (i = 0; i < getInterpOrder(); i++) {
				if ((dift = std::abs(x - X(startInd + i))) < dif) {
					ns = i;
					dif = dift;
				}
				c[i] = Y(startInd + i);
				d[i] = Y(startInd + i);
			}
			y = Y(startInd + ns--);				// This is the initial approximation to y
			
			for (m = 1; m < getInterpOrder(); m++)  	// For each column of the tableau
			{
				for (i = 0; i < getInterpOrder() - m; i++) 
				{
					// we loop over the current c’s and d’s and update
					ho = X(startInd + i) - x;
					hp = X(startInd + i + m) - x;
					w = c[i + 1] - d[i];
					
					if ((den = ho - hp) == 0.0) 
						throw("Poly_interp error");
					
					den = w / den;
					d[i] = hp * den;
					c[i] = ho * den;
				}
				dy = (2 * (ns + 1) < (getInterpOrder() - m) ? c[ns + 1] : d[ns--]);
				// After each column in the tableau is completed, we decide which correction, c or d, we
				//	want to add to our accumulating value of y, i.e., which path to take through the tableau
				//	— forking up or down.We do this in such a way as to take the most “straight line”
				//	route through the tableau to its apex, updating ns accordingly to keep track of where
				//	we are.This route keeps the partial approximations centered(insofar as possible) on
				//	the target x.The last dy added is thus the error indication.
				y += dy;
			}
			_errorEst = dy;
			return y;
		}
	};

	///////////////////////           CUBIC SPLINE INTERPOLATION						//////////////////////
	// Construct with x and y vectors, and (optionally) values of
	// the first derivative at the endpoints, then call interp for interpolated values.
	struct SplineInterpRealFunc : RealFunctionInterpolated
	{
		Vector<Real> _secDerY;

		SplineInterpRealFunc(Vector<Real>& xv, Vector<Real>& yv, Real yp1 = 1.e99, Real ypn = 1.e99)
			: RealFunctionInterpolated(xv, yv, 2), _secDerY(xv.size())
		{
			initSecDerivs(&xv[0], &yv[0], yp1, ypn);
		}

		// This routine stores an array _secDerY[0..numPoints-1] with second derivatives of the interpolating function
		// at the tabulated points pointed to by xv, using function values pointed to by yv. If yp1 and/or
		// ypn are equal to 1  1099 or larger, the routine is signaled to set the corresponding boundary
		// condition for a natural spline, with zero second derivative on that boundary; otherwise, they are
		// the values of the first derivatives at the endpoints.
		void initSecDerivs(const Real* xv, const Real* yv, Real yp1, Real ypn)
		{
			int		i, k;
			Real	p, qn, sig, un;

			int numPoints = (int)_secDerY.size();
			Vector<Real> u(numPoints - 1);

			if (yp1 > 0.99e99)
				_secDerY[0] = u[0] = 0.0;
			else 
			{
				_secDerY[0] = -0.5;
				u[0] = (3.0 / (xv[1] - xv[0])) * ((yv[1] - yv[0]) / (xv[1] - xv[0]) - yp1);
			}
			for (i = 1; i < numPoints - 1; i++) 
			{
				sig = (xv[i] - xv[i - 1]) / (xv[i + 1] - xv[i - 1]);
				p = sig * _secDerY[i - 1] + 2.0;

				_secDerY[i] = (sig - 1.0) / p;

				u[i] = (yv[i + 1] - yv[i]) / (xv[i + 1] - xv[i]) - (yv[i] - yv[i - 1]) / (xv[i] - xv[i - 1]);
				u[i] = (6.0 * u[i] / (xv[i + 1] - xv[i - 1]) - sig * u[i - 1]) / p;
			}
			if (ypn > 0.99e99)
				qn = un = 0.0;
			else 
			{
				qn = 0.5;
				un = (3.0 / (xv[numPoints - 1] - xv[numPoints - 2])) * (ypn - (yv[numPoints - 1] - yv[numPoints - 2]) / (xv[numPoints - 1] - xv[numPoints - 2]));
			}

			_secDerY[numPoints - 1] = (un - qn * u[numPoints - 2]) / (qn * _secDerY[numPoints - 2] + 1.0);

			for (k = numPoints - 2; k >= 0; k--)
				_secDerY[k] = _secDerY[k] * _secDerY[k + 1] + u[k];
		}

		// Given a value x, and using pointers to data xx and yy, and the stored vector of second derivatives
		// _secDerY, this routine returns the cubic spline interpolated value y.        
		Real calcInterpValue(int startInd, Real x) const
		{
			int indLow = startInd, indUpp = startInd + 1;
			Real y, h, b, a;
			h = X(indUpp) - X(indLow);
			
			if (h == 0.0) 
				throw("Bad input to routine splint");
			
			a = (X(indUpp) - x) / h;
			b = (x - X(indLow)) / h;

			y = a * Y(indLow) + b * Y(indUpp) + 
				  ((POW3(a) - a) * _secDerY[indLow] + (POW3(b) - b) * _secDerY[indUpp]) * (h * h) / 6.0;
			
			return y;
		}
	};

	//////////////////          PARAMETRIC CURVE SPLINE INTERPOLATION						/////////////////
	// Object for interpolating a curve specified by _numPoints points in N dimensions.
	template<int N>
	class SplineInterpParametricCurve : public IParametricCurve<N>
	{
		Real _minT, _maxT;
		int  _dim, _numPoints, _bemba;
		bool _isCurveClosed;

		Matrix<Real> _curvePoints;
		Vector<Real> s;
		Vector<Real> ans;

		std::vector<SplineInterpRealFunc*> srp;
	
	public:
		// Constructor. The _numPoints   _dim matrix ptsin inputs the data points. Input close as 0 for
		// an open curve, 1 for a closed curve. (For a closed curve, the last data point should not
		// duplicate the first — the algorithm will connect them.)
		SplineInterpParametricCurve(Real minT, Real maxT, const Matrix<Real>& ptsin, bool close = 0)
			: _numPoints(ptsin.RowNum()), _dim(ptsin.ColNum()), _bemba(close ? 2 * _numPoints : _numPoints),
			_isCurveClosed(close), _curvePoints(_dim, _bemba), s(_bemba), ans(_dim), srp(_dim),
			_minT(minT), _maxT(maxT)
		{
			// check N == dim
			if (N != _dim)
				throw("SplineInterpParametricCurve: N != dim");

			int i, ii, im, j, ofs;
			Real ss, soff, db, de;

			ofs = close ? _numPoints / 2 : 0;
			s[0] = 0.;
			for (i = 0; i < _bemba; i++) 
			{
				ii = (i - ofs + _numPoints) % _numPoints;
				im = (ii - 1 + _numPoints) % _numPoints;

				for (j = 0; j < _dim; j++) 
					_curvePoints[j][i] = ptsin[ii][j];
				
				if (i > 0) 
				{
					s[i] = s[i - 1] + rad(&ptsin[ii][0], &ptsin[im][0]);
					
					if (s[i] == s[i - 1]) 
						throw("error in Curve_interp");
					// Consecutive points may not be identical. For a closed curve, the last data
					// point should not duplicate the first.                    
				}
			}
			ss = close ? s[ofs + _numPoints] - s[ofs] : s[_numPoints - 1] - s[0];
			soff = s[ofs];
			
			for (i = 0; i < _bemba; i++) 
				s[i] = (s[i] - soff) / ss;
			
			for (j = 0; j < _dim; j++) 
			{
				db = _bemba < 4 ? 1.e99 : fprime(&s[0], &_curvePoints[j][0], 1);
				de = _bemba < 4 ? 1.e99 : fprime(&s[_bemba - 1], &_curvePoints[j][_bemba - 1], -1);

				Vector<Real> vec = _curvePoints.VectorFromRow(j);

				srp[j] = new SplineInterpRealFunc(s, vec, db, de);
			}
		}

		SplineInterpParametricCurve(const Matrix<Real>& ptsin, bool close = 0)
			: SplineInterpParametricCurve(0.0, 1.0, ptsin, close) {	}

		~SplineInterpParametricCurve() {
			for (int j = 0; j < _dim; j++) 
				delete srp[j];
		}
		
		Real getMinT() const { return _minT; }
		Real getMaxT() const { return _maxT; }

		// Interpolate a point on the stored curve. The point is parameterized by t, in the range [0,1].
		// For open curves, values of t outside this range will return extrapolations (dangerous!). For
		// closed curves, t is periodic with period 1
		VectorN<Real, N> operator()(Real t) const
		{
			if (t < _minT || t > _maxT)
				throw("SplineInterpParametricCurve: t outside interval");

			VectorN<Real, N> ans;

			if (_isCurveClosed)
				t = t - floor(t);

			// we have to map t from [minT, maxT] to [0, 1]
			t = (t - _minT) / (_maxT - _minT);
			for (int j = 0; j < _dim; j++)
				ans[j] = (*srp[j])(t);

			return ans;
		}

		// Utility for estimating the derivatives at the endpoints. x and y point to the abscissa and
		// ordinate of the endpoint. If pm is C1, points to the right will be used (left endpoint); if it
		// is -1, points to the left will be used (right endpoint). 
		Real fprime(Real* x, Real* y, int pm) {
			Real s1 = x[0] - x[pm * 1], s2 = x[0] - x[pm * 2], s3 = x[0] - x[pm * 3],
				s12 = s1 - s2, s13 = s1 - s3, s23 = s2 - s3;
			return -(s1 * s2 / (s13 * s23 * s3)) * y[pm * 3] + (s1 * s3 / (s12 * s2 * s23)) * y[pm * 2]
				- (s2 * s3 / (s1 * s12 * s13)) * y[pm * 1] + (1. / s1 + 1. / s2 + 1. / s3) * y[0];
		}

		Real rad(const Real* p1, const Real* p2) {
			Real sum = 0.;
			for (int i = 0; i < _dim; i++)
				sum += POW2(p1[i] - p2[i]);
			return sqrt(sum);
		}
	};


	template<int N>
	class InterpolatedSurface : public IParametricSurfaceRect<N>
	{
	public:
		InterpolatedSurface() {}

		VectorN<Real, N> operator()(const VectorN<Real, 2>& x) const { return VectorN<Real, N>{}; }
	};
}

#endif
