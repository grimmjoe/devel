#ifndef __ADD_H_
#define __ADD_H_

#include "function.h"
#include "strings.h"
#include <sstream>

namespace core
{
	template <class T>
	class add : public function<T>
	{
	public:
		typedef function<T> tBase;
		using typename tBase::tFunctionPtr;
	private:
		tFunctionPtr m_op1;
		tFunctionPtr m_op2;
	public:
		add(tFunctionPtr op1, tFunctionPtr op2)
			: m_op1(op1)
			, m_op2(op2)
		{}
		virtual T operator()(const T& t)
		{
			return (*m_op1)(t)+(*m_op2)(t);
		}
		virtual std::string toString() const
		{
			std::stringstream ss;
			ss << m_op1->toString() << strings::sPlus << m_op2->toString();
			return ss.str();
		}
		virtual tFunctionPtr derivative(int num)
		{
			if (num == 0)
				return this->clone();
			return tFunctionPtr(new add<T>(m_op1->derivative(num), m_op2->derivative(num)));
		}

		virtual void optimize()
		{
			// TODO - Nothing is here for now
		}

		virtual tFunctionPtr clone()
		{
			return tFunctionPtr(new add<T>(m_op1, m_op2));
		}
	}; // class add

}

#endif // __ADD_H_
