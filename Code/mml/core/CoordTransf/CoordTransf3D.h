#if !defined MML_COORD_TRANSF_3D_H
#define MML_COORD_TRANSF_3D_H

#include "MMLBase.h"

#include "core/CoordTransf.h"


namespace MML
{
	class CoordTransfCart3DRotationXAxis : 
		public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

	public:
		CoordTransfCart3DRotationXAxis(Real inAngle) : _angle(inAngle),
			_f1(std::function<Real(const VectorN<Real, 3>&)> 
						{ std::bind(&CoordTransfCart3DRotationXAxis::func1, this, std::placeholders::_1) }
				 ),
			_f2(std::function<Real(const VectorN<Real, 3>&)> 
						{ std::bind(&CoordTransfCart3DRotationXAxis::func2, this, std::placeholders::_1) }
				 ),
			_f3(std::function<Real(const VectorN<Real, 3>&)> 
						{ std::bind(&CoordTransfCart3DRotationXAxis::func3, this, std::placeholders::_1) }
				 ),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> 
						{ std::bind(&CoordTransfCart3DRotationXAxis::funcInverse1, this, std::placeholders::_1) }
				 ),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> 
						{ std::bind(&CoordTransfCart3DRotationXAxis::funcInverse2, this, std::placeholders::_1) }
				 ),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> 
						{ std::bind(&CoordTransfCart3DRotationXAxis::funcInverse3, this, std::placeholders::_1) }
				 )
		{
			_transf[0][0] = 1.0;
			_transf[1][1] = cos(_angle);
			_transf[1][2] = -sin(_angle);
			_transf[2][1] = sin(_angle);
			_transf[2][2] = cos(_angle);

			_inverse[0][0] = 1.0;
			_inverse[1][1] = cos(_angle);
			_inverse[1][2] = sin(_angle);
			_inverse[2][1] = -sin(_angle);
			_inverse[2][2] = cos(_angle);
		}

		// this could be made much, much more efficient :)
		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_inverse * q)[2]; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i)const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};
	class CoordTransfCart3DRotationYAxis : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

	public:
		CoordTransfCart3DRotationYAxis(Real inAngle) : _angle(inAngle),
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::funcInverse3, this, std::placeholders::_1) })
		{
			_transf[0][0] = cos(_angle);
			_transf[0][2] = sin(_angle);
			_transf[1][1] = 1.0;
			_transf[2][0] = -sin(_angle);
			_transf[2][2] = cos(_angle);

			_inverse[0][0] = cos(_angle);
			_inverse[0][2] = -sin(_angle);
			_inverse[1][1] = 1.0;
			_inverse[2][0] = sin(_angle);
			_inverse[2][2] = cos(_angle);
		}

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_inverse * q)[2]; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};
	class CoordTransfCart3DRotationZAxis : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

	public:
		CoordTransfCart3DRotationZAxis(Real inAngle) : _angle(inAngle),
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::funcInverse3, this, std::placeholders::_1) })
		{
			_transf[0][0] = cos(_angle);
			_transf[0][1] = -sin(_angle);
			_transf[1][0] = sin(_angle);
			_transf[1][1] = cos(_angle);
			_transf[2][2] = 1.0;

			_inverse[0][0] = cos(_angle);
			_inverse[0][1] = sin(_angle);
			_inverse[1][0] = -sin(_angle);
			_inverse[1][1] = cos(_angle);
			_inverse[2][2] = 1.0;
		}

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_inverse * q)[2]; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};

	// TODO 0.9 - LOW, TESKO, Euler angles, finish CoordTransfCart3DRotationGeneralAxis
	// class CoordTransfCart3DRotationGeneralAxis  : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	// {

	// };

	// Performs tranformation from original (Cartesian) system to orthogonal system defined by 
	// its base vectors expressed in original system.
	class CoordTransf3DCartOrthogonal : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Vector3Cartesian _base[3];

		MatrixNM<Real, 3, 3> _transf;
		MatrixNM<Real, 3, 3> _transfInverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[2]; }

	public:
		CoordTransf3DCartOrthogonal(const VectorN<Real, 3>& b1, const VectorN<Real, 3>& b2, const VectorN<Real, 3>& b3) :
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::funcInverse3, this, std::placeholders::_1) })
		{
			_base[0] = b1;
			_base[1] = b2;
			_base[2] = b3;

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++)
				{
					_transf(i, j) = _base[i][j];
					_transfInverse(i, j) = _base[j][i];
				}
			}
		}

		MatrixNM<Real, 3, 3> getTransfMatrix() { return _transf; }
		MatrixNM<Real, 3, 3> getInvTransfMatrix() { return _transfInverse; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		bool IsRightHanded()
		{
			Vector3Cartesian cross = VectorProduct(_base[0], _base[1]);
			if (ScalarProduct(cross, _base[2]) > 0.0)
				return true;
			else
				return false;
		}
	};

	// General 3D Cartesian transformation, given by matrix
	class CoordTransf3DCartGeneral : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Vector3Cartesian _base[3];

		MatrixNM<Real, 3, 3> _transf;
		MatrixNM<Real, 3, 3> _transfInverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[2]; }

	public:
		CoordTransf3DCartGeneral(const MatrixNM<Real, 3, 3>& transfMat) :
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::funcInverse3, this, std::placeholders::_1) })
		{
			_transf = transfMat;

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++)
				{
					_base[i][j] = _transf(i, j);
					_transfInverse(i, j) = _transf[j][i];
				}
			}

			// check if it is orthogonal!
			// TODO - implement this

		}

		MatrixNM<Real, 3, 3> getTransfMatrix() { return _transf; }
		MatrixNM<Real, 3, 3> getInvTransfMatrix() { return _transfInverse; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		bool IsRightHanded()
		{
			Vector3Cartesian cross = VectorProduct(_base[0], _base[1]);
			if (ScalarProduct(cross, _base[2]) > 0.0)
				return true;
			else
				return false;
		}
	};
	class CoordTransf3DCartOblique : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Vector3Cartesian _base[3];
		Vector3Cartesian _dual[3];

		MatrixNM<Real, 3, 3> _alpha;
		MatrixNM<Real, 3, 3> _transf;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		// TODO - ovo popraviti i napraviti kako spada
		Real func1(const VectorN<Real, 3>& q) const { return ScalarProduct(q, MML::Vector3Cartesian(_dual[0])); }
		Real func2(const VectorN<Real, 3>& q) const { return ScalarProduct(q, MML::Vector3Cartesian(_dual[1])); }
		Real func3(const VectorN<Real, 3>& q) const { return ScalarProduct(q, MML::Vector3Cartesian(_dual[2])); }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

	public:
		CoordTransf3DCartOblique(VectorN<Real, 3> b1, VectorN<Real, 3> b2, VectorN<Real, 3> b3) :
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::funcInverse3, this, std::placeholders::_1) })
		{
			_base[0] = b1;
			_base[1] = b2;
			_base[2] = b3;

			Real V = ScalarProduct(_base[0], VectorProduct(_base[1], _base[2]));

			_dual[0] = (1 / V) * VectorProduct(_base[1], _base[2]);
			_dual[1] = (1 / V) * VectorProduct(_base[2], _base[0]);
			_dual[2] = (1 / V) * VectorProduct(_base[0], _base[1]);

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					_alpha(i, j) = _base[i][j];
					_transf(i, j) = _base[j][i];     // transponirano
				}
			}
		}

		MatrixNM<Real, 3, 3> getAlpha() { return _alpha; }
		MatrixNM<Real, 3, 3> getTransf() { return _alpha; }

		Vector3Cartesian    Base(int i) { return _base[i]; }
		Vector3Cartesian    Dual(int i) { return _dual[i]; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		bool IsRightHanded()
		{
			Vector3Cartesian cross = VectorProduct(_base[0], _base[1]);
			if (ScalarProduct(cross, _base[2]) > 0.0)
				return true;
			else
				return false;
		}
	};
}

#endif
