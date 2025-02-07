#if !defined MML_ITENSOR_H
#define MML_ITENSOR_H

#include "MMLBase.h"

namespace MML
{
	enum TensorIndexType { CONTRAVARIANT, COVARIANT };
    
	template<int N>
	class ITensor2
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  operator()(int i, int j) const = 0;
		virtual Real& operator()(int i, int j) = 0;
	};

	template<int N>
	class ITensor3
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  operator()(int i, int j, int k) const = 0;
		virtual Real& operator()(int i, int j, int k) = 0;
	};

	template<int N>
	class ITensor4
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  operator()(int i, int j, int k, int l) const = 0;
		virtual Real& operator()(int i, int j, int k, int l) = 0;
	};

	template<int N>
	class ITensor5
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  operator()(int i, int j, int k, int l, int m) const = 0;
		virtual Real& operator()(int i, int j, int k, int l, int m) = 0;
	};
}
#endif


