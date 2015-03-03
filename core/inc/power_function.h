#ifndef __POWER_FUNCTION_H_
#define __POWER_FUNCTION_H_

#include "function.h"
#include <sstream>
#include <cmath>
#include <memory>

namespace core
{
	template <class T>
	class power : public function<T>
	{
	public:
		typedef function<T> tBase;
		using typename tBase::tFunctionPtr;
	private:
		tFunctionPtr m_base;
		T m_power;
	public:
		power(tFunctionPtr b, T p)
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
		virtual tFunctionPtr derivative(int num)
		{
			tFunctionPtr thisFunc = this->clone();
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

		virtual tFunctionPtr clone()
		{
			return tFunctionPtr(new power<T>(m_base, m_power));
		}
	};
}

#endif // __POWER_FUNCTION_H_
