#ifndef _COMPARATOR_H__
#define _COMPARATOR_H__

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
			// TODO - This should be implemented the right way
			return a == b;
		}
	};
} // namespace core

#endif // _COMPARATOR_H__
