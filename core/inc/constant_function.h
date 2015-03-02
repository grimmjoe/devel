#ifndef __CONSTANT_FUNCTION_H_
#define __CONSTANT_FUNCTION_H_

#include "function.h"
#include <sstream>

namespace core
{
	template <class T>
	class const_function : public function<T>
	{
		T m_value;
	public:
		const_function(const T& b)
			: m_value(b)
		{}
		virtual T operator()(const T&)
		{
			return m_value;
		}
		virtual std::string toString() const
		{
			std::stringstream ss;
			ss << m_value;
			return ss.str();
		}
		virtual std::shared_ptr<function<T> > derivative(int num)
		{
			if (num > 0)
				return std::shared_ptr<functn<T> >(new const_function<T>(0));
			else if (num == 0)
				return this->clone();
			else
				return std::shared_ptr<function<T> >(nullptr);
		}

		virtual void optimize()
		{
			// TODO - Nothing is here for now
		}

		virtual std::shared_ptr<function<T> > clone()
		{
			return std::shared_ptr<function<T> >(new const_function<T>(m_value));
		}
	};
}

#endif // __CONSTANT_FUNCTION_H_
