#if !defined MML_STATISTICS_H
#define MML_STATISTICS_H

#include "MMLBase.h"

#include "base/Vector.h"

namespace MML
{
	class Statistics
	{
	public:
		static Real Avg(const Vector<Real>& data)
		{
			Real outAvg = 0.0;
			int n = data.size();
			
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in Avg");

			for (int j = 0; j < n; j++)
				outAvg += data[j];
			outAvg /= n;

			return outAvg;
		}

		static void AvgVar(const Vector<Real>& data, Real& outAvg, Real& outVar)
		{
			Real s, ep;
			int j, n = data.size();

			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in AvgVar");

			outAvg = Avg(data);

			outVar = ep = 0.0;
			for (j = 0; j < n; j++) {
				s = data[j] - outAvg;
				ep += s;
				outVar += s * s;
			}
			outVar = (outVar - ep * ep / n) / (n - 1);
		}
		static void AvgStdDev(const Vector<Real>& data, Real& outAvg, Real& outStdDev)
		{
			Real var;
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in AvgStdDev");
			AvgVar(data, outAvg, var);
			outStdDev = sqrt(var);
		}


		static void Moments(const Vector<Real>& data, Real& ave, Real& adev, Real& sdev, Real& var, Real& skew, Real& curt)
		{
			int j, n = data.size();
			Real ep = 0.0, s, p;

			if (n <= 1)
				throw StatisticsError("Vector size must be greater than 1 in Moments");

			s = 0.0;
			for (j = 0; j < n; j++)
				s += data[j];
			ave = s / n;

			adev = var = skew = curt = 0.0;
			for (j = 0; j < n; j++) {
				adev += std::abs(s = data[j] - ave);
				ep += s;
				var += (p = s * s);
				skew += (p *= s);
				curt += (p *= s);
			}
			adev /= n;
			var = (var - ep * ep / n) / (n - 1);
			sdev = sqrt(var);

			if (var != 0.0) {
				skew /= (n * var * sdev);
				curt = curt / (n * var * var) - 3.0;
			}
			else
				throw("No skew/kurtosis when variance = 0 (in Moments)");
		}
	};
}
#endif