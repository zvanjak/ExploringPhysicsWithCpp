#if !defined MPL_BASE_H
#define MPL_BASE_H


namespace MPL
{
	enum CollisionSimulatorRunType
	{
		RunTypeExact = 0,					// exact simulation, 
		RunTypeFast,							// fast simulation, with subdivision of the box into subcontainers,
		RunTypeFastMultithread,		// fast simulation, with subdivision and using multithreading
		RunTypeFastMultithreadTP	// fast simulation, with subdivision and using multithreading with thread pool
	};
}

#endif // MPLBase.h