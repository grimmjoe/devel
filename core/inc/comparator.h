#ifndef _COMPARATOR_H__
#define _COMPARATOR_H__

#include <cmath>

namespace core
{
	template <class T>
	class comparator
	{
		T m_threshold;
	public:
		comparator(T t)
			: m_threshold(t)
		{}

		bool is_equal(T a, T b) const
		{
			return std::fabs(a-b) <= m_threshold;
		}
	};
} // namespace core

#endif // _COMPARATOR_H__
