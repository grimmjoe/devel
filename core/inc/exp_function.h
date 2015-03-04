#ifndef __EXP_FUNCTION_H_
#define __EXP_FUNCTION_H_

#include "function.h"
#include "multiply.h"
#include <sstream>
#include <cmath>

namespace core
{
	template <class T>
	class exp : public function<T>
	{
	public:
		typedef function<T> tBase;
		using typename tBase::tFunctionPtr;
	private:
		T m_base;
		tFunctionPtr m_power;
	public:
		exp(T b, tFunctionPtr p)
			: m_base(b)
			, m_power(p)
		{}
		virtual T operator()(const T& t)
		{
			return std::pow(m_base, (*m_power)(t));
		}
		virtual std::string toString() const
		{
			// TODO
			//if (m_power == 0)
			//	return "1";
			//else if (m_power == 1)
			//	return m_base.toString();

			std::stringstream ss;
			ss << m_base << "^" << "(" << m_power->toString() << ")";
			return ss.str();
		}
		virtual tFunctionPtr derivative(int num)
		{
			if (num == 0)
				return this->clone();
			if (num == 1)
				return tFunctionPtr(new multiply<T>(
						tFunctionPtr(new multiply<T, eMultNum>(std::log(m_base), this->clone())),
						m_power->derivative(1)));
			tFunctionPtr thisFunc = this->clone();
			while (num > 0)
			{
				--num;
				thisFunc = thisFunc->derivative(1);
			}
			return thisFunc;
		}

		virtual void optimize()
		{
			// TODO - Nothing is here for now
		}

		virtual tFunctionPtr clone()
		{
			return tFunctionPtr(new exp<T>(m_base, m_power));
		}
	};
}

#endif // __EXP_FUNCTION_H_
