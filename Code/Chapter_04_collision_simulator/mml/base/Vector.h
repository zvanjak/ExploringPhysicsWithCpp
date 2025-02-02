#if !defined  MML_Vector_H
#define MML_Vector_H

#include "MMLBase.h"

namespace MML
{
	template<class Type>
	class	Vector
	{
	private:
		std::vector<Type> _elems;

	public:
		typedef Type value_type;      // make T available externally

		/////////////////////                  Constructors                ////////////////////
		Vector() {}
		explicit Vector(int n) {
			if(n < 0)
				throw VectorInitializationError("Vector::Vector - negative size", n);

			_elems.resize(n);
		}
		explicit Vector(int n, const Type &val) {
			if (n < 0)
				throw VectorInitializationError("Vector::Vector - negative size", n);

			_elems.resize(n, val);
		}
		explicit Vector(int n, Type* vals) 
		{
			if (n < 0)
				throw VectorInitializationError("Vector::Vector - negative size", n);

			_elems.resize(n);
			for (int i = 0; i < n; ++i)
				_elems[i] = vals[i];
		}
		explicit Vector(std::vector<Type> values) : _elems(values) {}
		explicit Vector(std::initializer_list<Type> list) : _elems(list) {}

		static Vector GetUnitVector(int dimVec, int indUnit)
		{
			if (indUnit < 0 || indUnit >= dimVec)
				throw VectorDimensionError("Vector::GetUnitVector - wrong unit index", dimVec, indUnit);

			Vector ret(dimVec);
			ret[indUnit] = Type{ 1.0 };
			return ret;
		}

		// not really needed, but let's be explicit
		Vector(const Vector& b) = default; 
		Vector(Vector&& b) = default;
		Vector& operator=(const Vector& b) = default; 
		Vector& operator=(Vector&& b) = default;

		////////////////            std::vector forwarding                 ////////////////////
		int  size()  const { return (int)_elems.size(); }
		bool empty() const { return _elems.empty(); }

		void clear()						{ _elems.clear(); }
		void resize(int newLen)	{ _elems.resize(newLen); }

		/////////////////////            Accessing elements             ///////////////////////
		inline Type&       operator[](int n)       { return _elems[n]; }
		inline const Type& operator[](int n) const { return _elems[n]; }

		// checked access
		Type& at(int n)	{
			if(n < 0 || n >= size())
				throw VectorDimensionError("Vector::at - index out of bounds", size(), n);
			else
				return _elems[n];
		}
		Type  at(int n) const { 
			if(n < 0 || n >= size())
				throw VectorDimensionError("Vector::at - index out of bounds", size(), n);
			else
				return _elems[n];
		}

		///////////////////////          Arithmetic operators         ////////////////////////
		Vector  operator-() const         // unary minus
		{
			Vector ret(size());
			for (int i = 0; i < size(); i++)
				ret._elems[i] = Type{ -1 } * (*this)[i];
			return ret;
		}
		
		Vector  operator+(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator+() - vectors must be equal size", size(), b.size());

			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = (*this)[i] + b._elems[i];
			return ret;
		}
		Vector& operator+=(const Vector& b)
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator+=() - vectors must be equal size", size(), b.size());
			
			for (int i = 0; i < b.size(); i++)
				_elems[i] += b._elems[i];
			return *this;
		}
		Vector  operator-(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator-() - vectors must be equal size", size(), b.size());

			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = (*this)[i] - b._elems[i];
			return ret;
		}
		Vector& operator-=(const Vector& b)
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator-=() - vectors must be equal size", size(), b.size());
			
			for (int i = 0; i < b.size(); i++)
				_elems[i] -= b._elems[i];
			return *this;
		}

		Vector  operator*(Type b)
		{
			Vector ret(size());;
			for (int i = 0; i < size(); i++)
				ret._elems[i] = b * _elems[i];
			return ret;
		}
		Vector& operator*=(Type b)
		{
			for (int i = 0; i < size(); i++)
				_elems[i] *= b;
			return *this;
		}
		Vector  operator/(Type b)
		{
			Vector ret(size());
			for (int i = 0; i < size(); i++)
				ret._elems[i] = _elems[i] / b;
			return ret;
		}
		Vector& operator/=(Type b)
		{
			for (int i = 0; i < size(); i++)
				_elems[i] /= b;
			return *this;
		}
		
		friend Vector operator*(Type a, const Vector& b)
		{
			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = a * b._elems[i];
			return ret;
		}

		///////////////////////             Testing equality             ////////////////////////
		bool operator==(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator==() - vectors must be equal size", size(), b.size());

			for (int i = 0; i < size(); i++)
				if ((*this)[i] != b[i])
					return false;

			return true;
		}
		bool operator!=(const Vector& b) const
		{
			return !(*this == b);
		}
		bool IsEqualTo(const Vector& b, Real eps = Defaults::VectorEqualityPrecision) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::IsEqual - vectors must be equal size", size(), b.size());

			for (int i = 0; i < size(); i++)
			{
				if (Abs((*this)[i] - b[i]) > eps)
					return false;
			}
			return true;
		}	
		bool IsNullVec() const
		{
			for (int i = 0; i < size(); i++)
				if (Abs((*this)[i]) != 0.0 )
					return false;
			return true;
		}
		
		//////////////////////                 Operations                 ///////////////////////
		Real NormL1() const
		{
			Real norm{ 0.0 };
			for (int i = 0; i < size(); i++)
				norm += Abs((*this)[i]);
			return norm;
		}
		Real NormL2() const
		{
			Real norm{ 0.0 };
			for (int i = 0; i < size(); i++)
				norm += (*this)[i] * (*this)[i];
			return std::sqrt(norm);
		}
		Real NormLInf() const
		{
			Real norm{ 0.0 };
			for (int i = 0; i < size(); i++)
				norm = std::max(norm, Abs((*this)[i]));
			return norm;
		}

		///////////////////////////               I/O                 ///////////////////////////
		std::ostream& Print(std::ostream& stream, int width, int precision) const
		{
			stream << "[";
			bool first = true;
			for (const Type& x : _elems)
			{
				if (!first)
					stream << ", ";
				else
					first = false;

				stream << std::setw(width) << std::setprecision(precision) << x;
			}
			stream << "]";

      return stream;
		}
		std::ostream& Print(std::ostream& stream, int width, int precision, Real zeroThreshold) const
		{
			stream << "[";
			bool first = true;
			for (const Type& x : _elems)
			{
				if (!first)
					stream << ", ";
				else
					first = false;

        if( Abs(x) > zeroThreshold )
				  stream << std::setw(width) << std::setprecision(precision) << x;
        else
          stream << std::setw(width) << std::setprecision(precision) << 0.0;
			}
			stream << "]";

			return stream;
		}    
    std::ostream& PrintLine(std::ostream& stream, const std::string &msg, int width, int precision) const
		{
			stream << msg;
			Print(stream, width, precision);
			stream << std::endl;

			return stream;
		}
		
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}    
		friend std::ostream& operator<<(std::ostream &stream, const Vector &a)
		{
			a.Print(stream, Defaults::VectorPrintWidth, Defaults::VectorPrintPrecision);

			return stream;
		}  
	};

	typedef Vector<int>     VectorInt;
	typedef Vector<float>   VectorFlt;
	typedef Vector<double>  VectorDbl;
	typedef Vector<Complex> VectorComplex;

	typedef Vector<int>     VecI;
	typedef Vector<float>   VecF;
	typedef Vector<double>  VecD;
	typedef Vector<Complex> VecC;
}

#endif // MML_Vector_H