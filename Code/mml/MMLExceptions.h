#if !defined MML_EXCEPTIONS_H
#define MML_EXCEPTIONS_H

#include <stdexcept>

namespace MML
{
	//////////             Vector error exceptions            ///////////
	class VectorInitializationError : public std::invalid_argument
	{
	public:
		int _size1;
		VectorInitializationError(std::string inMessage, int size1) : std::invalid_argument(inMessage), _size1(size1)
		{ }
	};
	class VectorDimensionError : public std::invalid_argument
	{
	public:
		int _size1, _size2;
		VectorDimensionError(std::string inMessage, int size1, int size2) : std::invalid_argument(inMessage), _size1(size1), _size2(size2)
		{ }
	};
	class VectorAccessBoundsError : public std::out_of_range
	{
	public:
		int _i, _n;
		VectorAccessBoundsError(std::string inMessage, int i, int n) : std::out_of_range(inMessage), _i(i), _n(n)
		{ }
	};

	//////////             Matrix error exceptions            ///////////
	class MatrixAllocationError : public std::out_of_range
	{
	public:
		int _rows, _cols;
		MatrixAllocationError(std::string inMessage, int rows, int cols) : std::out_of_range(inMessage), _rows(rows), _cols(cols)
		{ }
	};
	class MatrixAccessBoundsError : public std::out_of_range
	{
	public:
		int _i, _j, _rows, _cols;
		MatrixAccessBoundsError(std::string inMessage, int i, int j, int rows, int cols) : std::out_of_range(inMessage), _i(i), _j(j), _rows(rows), _cols(cols)
		{ }
	};
	class MatrixDimensionError : public std::invalid_argument
	{
	public:
		int _rows1, _cols1, _rows2, _cols2;

		MatrixDimensionError(std::string inMessage, int r1, int c1, int r2, int c2) : std::invalid_argument(inMessage), _rows1(r1), _cols1(c1), _rows2(r2), _cols2(c2)
		{ }
	};
	class SingularMatrixError : public std::domain_error
	{
	public:
		SingularMatrixError(std::string inMessage) : std::domain_error(inMessage)
		{ }
	};

	//////////             Integration exceptions            ///////////
	class IntegrationTooManySteps : public std::domain_error
	{
	public:
		IntegrationTooManySteps(std::string inMessage) : std::domain_error(inMessage)
		{ }
	};

	//////////           Interpolation exceptions            ///////////
	class RealFuncInterpInitError: public std::domain_error
	{
	public:
		RealFuncInterpInitError(std::string inMessage) : std::domain_error(inMessage)
		{ }
	};

	////////////             Tensor exceptions             /////////////
	class TensorCovarContravarNumError : public std::invalid_argument
	{
	public:
		int _numContra, _numCo;
		TensorCovarContravarNumError(std::string inMessage, int size1, int size2) : std::invalid_argument(inMessage), _numContra(size1), _numCo(size2)
		{ }
	};
	class TensorCovarContravarArithmeticError : public std::invalid_argument
	{
	public:
		int _numContra, _numCo;
        int _bContra, _bCo;

		TensorCovarContravarArithmeticError(std::string inMessage, int contra, int co, int b_contra, int b_co) : std::invalid_argument(inMessage), _numContra(contra), _numCo(co), _bContra(b_contra), _bCo(b_co)
		{ }
	};
	class TensorIndexError : public std::invalid_argument
	{
	public:
		TensorIndexError(std::string inMessage) : std::invalid_argument(inMessage)
		{ }
	};    

	//////////           Root finding exceptions            ///////////
	class RootFindingError: public std::runtime_error
	{
	public:
		RootFindingError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

	//////////            ODE solver exceptions             ///////////
	class ODESolverError: public std::runtime_error
	{
	public:
		ODESolverError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

	//////////            Statistics exceptions             ///////////
	class StatisticsError : public std::runtime_error
	{
	public:
		StatisticsError(std::string inMessage) : std::runtime_error(inMessage)
		{
		}
	};

}

#endif // MML_EXCEPTIONS_H