#if !defined MML_ICoordSystemTransf_H
#define MML_ICoordSystemTransf_H

#include "MMLBase.h"

#include "IFunction.h"

#include "base/VectorN.h"

namespace MML
{
	template<typename VectorFrom, typename VectorTo, int N>
	class ICoordTransf : public IVectorFunction<N>
	{
	public:
		virtual       VectorTo            transf(const VectorFrom& in) const = 0;
		virtual const IScalarFunction<N>& coordTransfFunc(int i) const = 0;

		virtual ~ICoordTransf() {}
	};

	template<typename VectorFrom, typename VectorTo, int N>
	class ICoordTransfWithInverse : public virtual ICoordTransf<VectorFrom, VectorTo, N>
	{
	public:
		virtual       VectorFrom          transfInverse(const VectorTo& in) const = 0;
		virtual const IScalarFunction<N>& inverseCoordTransfFunc(int i) const = 0;

		virtual ~ICoordTransfWithInverse() {}
	};
}
#endif