#if !defined MML_IODESYSTEM_H
#define MML_IODESYSTEM_H

#include "MMLBase.h"

#include "base/Vector.h"

namespace MML
{
	class IODESystem
	{
	public:
		virtual int   getDim() const = 0;
		virtual void  derivs(const Real t, const Vector<Real> &x, Vector<Real> &dxdt) const = 0;

		// overridable function for providing variable names
		virtual std::string getVarName(int ind) const
		{
			if (ind < 0 || ind >= getDim())
				return "var" + std::to_string(ind);

			return "var" + std::to_string(ind + 1); // 1-based index 
		}
	};

	class IODESystemParametrized : public IODESystem
	{
	public:
		virtual int		getNumParam() const = 0;
		virtual Real	getParam(int i) const = 0;
		virtual void	setParam(int i, Real val) = 0;

		virtual Vector<Real>	getParams() const = 0;
		virtual void					setParams(const Vector<Real>&) = 0;

	};

}
#endif