#if !defined MML_INTERVALS_H
#define MML_INTERVALS_H

#include "MMLBase.h"

#include "interfaces/IInterval.h"

namespace MML
{
	// TODO - finalize intervals properly (and use them in test beds)
	///////////////////////////////////////////////   Interfaces    ///////////////////////////////////////////
	enum class EndpointType
	{
		OPEN,
		CLOSED,
		NEG_INF,
		POS_INF
	};

	class BaseInterval : public IInterval
	{
	protected:
		Real _lower, _upper;
		EndpointType _lowerType, _upperType;

		BaseInterval(Real lower, EndpointType lowerType, Real upper, EndpointType upperType) : _lower(lower), _lowerType(lowerType), _upper(upper), _upperType(upperType) { }
	public:
		virtual ~BaseInterval() {}

		Real getLowerBound() const { return _lower; }
		Real getUpperBound() const { return _upper; }
		Real getLength()     const { return _upper - _lower; }

		virtual bool isContinuous()  const { return true; }     // we suppose continuous intervals by default

		void GetEquidistantCovering(int numPoints) 
		{ 
			// TODO 1.1 - implement equidistant covering
		}
	};

	class CompleteRInterval : public BaseInterval
	{
	public:
		CompleteRInterval() : BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, std::numeric_limits<double>::max(), EndpointType::POS_INF) { }
		bool contains(Real x) const { return true; }
	};
	class CompleteRWithReccuringPointHoles : public BaseInterval
	{
		Real _hole0, _holeDelta;
	public:
		CompleteRWithReccuringPointHoles(Real hole0, Real holeDelta) : BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, std::numeric_limits<double>::max(), EndpointType::POS_INF)
		{
			_hole0 = hole0;
			_holeDelta = holeDelta;
		}

		bool contains(Real x) const {


			if (x == _hole0)
				return false;

			Real diff = (x - _hole0) / _holeDelta;
			if (diff == (int)diff)
				return false;

			return true;
		}
		bool isContinuous() const { return false; }
	};
	class OpenInterval : public BaseInterval
	{
		// for equidistant covering
		double _lowerRealDif = 0.0000001;
	public:
		OpenInterval(Real lower, Real upper) : BaseInterval(lower, EndpointType::OPEN, upper, EndpointType::OPEN) { }
		bool contains(Real x) const { return (x > _lower) && (x < _upper); }
	};
	class OpenClosedInterval : public BaseInterval
	{
	public:
		OpenClosedInterval(Real lower, Real upper) : BaseInterval(lower, EndpointType::OPEN, upper, EndpointType::CLOSED) { }
		bool contains(Real x) const { return (x > _lower) && (x <= _upper); }
	};
	class ClosedInterval : public BaseInterval
	{
	public:
		ClosedInterval(Real lower, Real upper) : BaseInterval(lower, EndpointType::CLOSED, upper, EndpointType::CLOSED) { }

		bool contains(Real x) const {
			return (x >= _lower) && (x <= _upper);
		}
	};
	class ClosedIntervalWithReccuringPointHoles : public BaseInterval
	{
		Real _hole0, _holeDelta;
	public:
		ClosedIntervalWithReccuringPointHoles(Real lower, Real upper, Real hole0, Real holeDelta)
			: BaseInterval(lower, EndpointType::CLOSED, upper, EndpointType::CLOSED)
		{
			_hole0 = hole0;
			_holeDelta = holeDelta;
		}

		bool contains(Real x) const {
			// check for hole!
			if (x == _hole0)
				return false;

			Real diff = (x - _hole0) / _holeDelta;
			if (diff == (int)diff)
				return false;

			return (x >= _lower) && (x <= _upper);
		}
		bool isContinuous() const { return false; }
	};
	class ClosedOpenInterval : public BaseInterval
	{
	public:
		ClosedOpenInterval(Real lower, Real upper) : BaseInterval(lower, EndpointType::CLOSED, upper, EndpointType::OPEN) { }

		bool contains(Real x) const {
			return (x >= _lower) && (x < _upper);
		}
	};
	class NegInfToOpenInterval : public BaseInterval
	{
	public:
		NegInfToOpenInterval(Real upper) : BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, upper, EndpointType::OPEN) { }

		bool contains(Real x) const {
			return x < _upper;
		}
	};
	class NegInfToClosedInterval : public BaseInterval
	{
	public:
		NegInfToClosedInterval(Real upper) : BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, upper, EndpointType::CLOSED) { }

		bool contains(Real x) const {
			return x <= _upper;
		}
	};
	class OpenToInfInterval : public BaseInterval
	{
	public:
		OpenToInfInterval(Real lower) : BaseInterval(lower, EndpointType::OPEN, std::numeric_limits<double>::max(), EndpointType::POS_INF) { }

		bool contains(Real x) const {
			return x > _lower;
		}
	};
	class ClosedToInfInterval : public BaseInterval
	{
	public:
		ClosedToInfInterval(Real lower) : BaseInterval(lower, EndpointType::CLOSED, std::numeric_limits<double>::max(), EndpointType::POS_INF) { }

		bool contains(Real x) const {
			return x >= _lower;
		}
	};

	class Interval : public IInterval
	{
		Real _lower, _upper;
		std::vector<std::shared_ptr<BaseInterval>> _intervals;
	public:
		Interval() {}
		Interval(std::initializer_list<BaseInterval*> intervals)
		{
			for (BaseInterval* interval : intervals)
			{
				_intervals.emplace_back(std::shared_ptr<BaseInterval>(interval));
			}

			// sortirati po lower bound
			// i redom provjeriti presijecanja
		}

		template<class _IntervalType>
		Interval& AddInterval(const _IntervalType& interval)
		{
			_intervals.emplace_back(std::make_shared<_IntervalType>(interval));
			return *this;
		}

		static Interval Intersection(const IInterval& a, const IInterval& b)
		{
			Interval ret;
			// TODO - implement intersection
			return ret;
		}
		static Interval Difference(const IInterval& a, const IInterval& b)
		{
			Interval ret;
			// TODO - implement diff
			return ret;
		}
		static Interval Complement(const IInterval& a)
		{
			Interval ret;
			// TODO - implement complement
			return ret;
		}
		Real getLowerBound() const { return 0; }
		Real getUpperBound() const { return 0; }
		Real getLength()     const { return 0; }

		bool isContinuous()  const { return false; }
		bool contains(Real x) const
		{
			// check for each interval if it contains x
			for (auto& interval : _intervals)
			{
				if (interval->contains(x))
					return true;
			}
			return false;
		}
		// TODO - implement this
		// bool contains(const IInterval &other) const = 0;
		// bool intersects(const IInterval &other) const = 0;

		void GetEquidistantCovering(int numPoints) { }
	};
}

#endif