#ifndef __EXP_FUNCTION_H_
#define __EXP_FUNCTION_H_

#include "function.h"
#include <sstream>

namespace core
{
	template <class T>
	class exp : public function<T>
	{
		T m_base;
		std::shared_ptr<function<T> > m_power;
	public:
		exp(T b, std::shared_ptr<function<T> > p)
			: m_base(b)
			, m_power(p)
		{}
		virtual T operator()(const T& t)
		{
			return std::pow(m_base, (*m_power)(t));
		}
		virtual std::string toString() const
		{
			if (m_power == 0)
				return "1";
			else if (m_power == 1)
				return m_base.toString();

			std::stringstream ss;
			ss << m_base << "^" << m_power.toString();
			return ss.str();
		}
		virtual std::shared_ptr<function<T> > derivative(int num)
		{
			std::shared_ptr<function<T> > thisFunc = this->clone();
			while (num > 0)
			{
				--num;
			}
			return thisFunc;
		}

		virtual void optimize()
		{
			// TODO - Nothing is here for now
		}

		virtual std::shared_ptr<function<T> > clone()
		{
			return std::shared_ptr<function<T> >(new exp<T>(m_base, m_power));
		}
	};
}

#endif // __EXP_FUNCTION_H_
