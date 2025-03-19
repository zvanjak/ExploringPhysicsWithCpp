#if !defined MML_IINTERVAL_H
#define MML_IINTERVAL_H

#include "MMLBase.h"

namespace MML
{
	// group
	class IInterval
	{
	public:
		virtual Real getLowerBound() const = 0;
		virtual Real getUpperBound() const = 0;
		virtual Real getLength() const = 0;

		virtual bool isContinuous() const = 0;
		virtual bool contains(Real x) const = 0;

		virtual void GetEquidistantCovering(int numPoints) = 0;

		virtual ~IInterval() {}
	};
}
#endif