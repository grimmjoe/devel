#ifndef __POWER_FUNCTION_H_
#define __POWER_FUNCTION_H_

#include "function.h"
#include <sstream>

namespace core
{
	template <class T>
	class power : public function<T>
	{
		std::shared_ptr<function<T> > m_base;
		double m_power;
	public:
		explicit power(std::shared_ptr<function<T> > b, double p)
			: m_base(b)
			, m_power(p)
		{}
		virtual T operator()(const T& t)
		{
			return std::pow((*m_base)(t), m_power);
		}
		virtual std::string toString() const
		{
			if (m_power == 0)
				return "1";
			else if (m_power == 1)
				return m_base.toString();

			std::stringstream ss;
			ss << m_base.toString() << "^" << m_power;
			return ss.str();
		}
		virtual std::shared_ptr<function<T> > derivative(int num)
		{
			std::shared_ptr<function<T> > thisFunc = this->clone();
			while (num > 0)
			{
				--num;
			}
		}

		virtual void optimize()
		{
			// TODO - Nothing is here for now
		}

		virtual std::shared_ptr<function<T> > clone()
		{
			return std::shared_ptr<function<T> >(new power<T>(m_base, m_power));
		}
	};
}

#endif // __POWER_FUNCTION_H_
