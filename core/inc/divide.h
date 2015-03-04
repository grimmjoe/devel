#ifndef __DIVIDE_H_
#define __DIVIDE_H_

#include "function.h"
#include "strings.h"
#include "multiply.h"
#include "add.h"
#include "subtract.h"
#include <sstream>

namespace core
{
	template <class T>
	class divide : public function<T>
	{
	public:
		typedef function<T> tBase;
		using typename tBase::tFunctionPtr;
	private:
		tFunctionPtr m_op1;
		tFunctionPtr m_op2;
	public:
		divide(tFunctionPtr op1, tFunctionPtr op2)
			: m_op1(op1)
			, m_op2(op2)
		{}
		virtual T operator()(const T& t)
		{
			return (*m_op1)(t)/(*m_op2)(t);
		}
		virtual std::string toString() const
		{
			std::stringstream ss;
			ss << "(" << m_op1->toString() << ")" << strings::sDiv << "(" << m_op2->toString() << ")";
			return ss.str();
		}
		virtual tFunctionPtr derivative(int num)
		{
			if (num == 0)
				return this->clone();
			tFunctionPtr ptr(new divide<T>(
					tFunctionPtr(new subtract<T>(
						tFunctionPtr(new multiply<T>(m_op1->derivative(1), m_op2)),
						tFunctionPtr(new multiply<T>(m_op1, m_op2->derivative(1)))
						)),
					tFunctionPtr(new power<T>(m_op2, 2))
					));
			return ptr->derivative(num-1);
		}

		virtual void optimize()
		{
			// TODO - Nothing is here for now
		}

		virtual tFunctionPtr clone()
		{
			return tFunctionPtr(new divide<T>(m_op1, m_op2));
		}
	}; // class divide
}

#endif // __DIVIDE_H_
