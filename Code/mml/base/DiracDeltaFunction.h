#if !defined MML_DELTA_FUNCTION_H
#define MML_DELTA_FUNCTION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

namespace MML
{
	class DiracFunction : public IRealFunction
	{
	protected:
		int _N;
	public:
		DiracFunction(int N) : _N(N) {}
	};

	class DiracStep : public DiracFunction
	{
	public:
		DiracStep(int N) : DiracFunction(N) {}

		Real operator()(const Real x) const
		{
			if (x < -1.0 / (2 * _N) || x > 1.0 / (2 * _N))
				return 0.0;
			else
				return _N;
		}
	};
	class DiracExp : public DiracFunction
	{
	public:
		DiracExp(int N) : DiracFunction(N) {}

		Real operator()(const Real x) const { return _N / sqrt(2 * Constants::PI) * exp(-x * x * _N * _N); }
	};
	class DiracSqr : public DiracFunction
	{
	public:
		DiracSqr(int N) : DiracFunction(N) {}

		Real operator()(const Real x) const { return _N / Constants::PI / (1 + _N * _N * x * x); }
	};
	class DiracSin : public DiracFunction
	{
	public:
		DiracSin(int N) : DiracFunction(N) {}

		Real operator()(const Real x) const { return sin(_N * x) / (Constants::PI * x); }
	};
}

#endif
