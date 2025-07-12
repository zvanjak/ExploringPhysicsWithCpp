#if !defined MML_VOLUME_INTEGRATION_H
#define MML_VOLUME_INTEGRATION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/Geometry3D.h"

#include "core/Integration.h"

namespace MML
{
	class VolumeIntegration
	{

		// moze se definirati od tijela, koe ima funkciju isWithin(), i onda komponirati

	public:
		static Real VolumeIntegral(const IScalarFunction<3>& scalarField, const BodyWithRectSurfaces& solid)
		{
			return 0.0;
		}
	};
} // end namespace
#endif