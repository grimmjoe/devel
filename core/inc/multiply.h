#ifndef __MULTIPLY__H_
#define __MULTIPLY__H_

#include "function.h"
#include "add.h"
#include "strings.h"
#include <sstream>

namespace core
{
	enum eMultType
	{
		eMultFunc = 0,
		eMultNum = 1
	};
	template <class T, eMultType MT = eMultFunc>
	class multiply : public function<T>
	{
	public:
		typedef function<T> tBase;
		using typename tBase::tFunctionPtr;
	private:
		tFunctionPtr m_op1;
		tFunctionPtr m_op2;
	public:
		multiply(tFunctionPtr op1, tFunctionPtr op2)
			: m_op1(op1)
			, m_op2(op2)
		{}
		virtual T operator()(const T& t)
		{
			return (*m_op1)(t)*(*m_op2)(t);
		}
		virtual std::string toString() const
		{
			std::stringstream ss;
			ss << "(" << m_op1->toString() << ")" << strings::sMult << "(" << m_op2->toString() << ")";
			return ss.str();
		}
		virtual tFunctionPtr derivative(int num)
		{
			return tFunctionPtr(new add<T>(
						tFunctionPtr(new multiply<T>(m_op1->derivative(num), m_op2)),
						tFunctionPtr(new multiply<T>(m_op1, m_op2->derivative(num)))));
		}

		virtual void optimize()
		{
			// TODO - Nothing is here for now
		}

		virtual tFunctionPtr clone()
		{
			return tFunctionPtr(new multiply<T>(m_op1, m_op2));
		}
	}; // class multiply

	template <class T>
	class multiply<T, eMultNum> : public function<T>
	{
	public:
		typedef function<T> tBase;
		using typename tBase::tFunctionPtr;
	private:
		T m_base;
		tFunctionPtr m_op2;
	public:
		multiply(T b, tFunctionPtr op2)
			: m_base(b)
			, m_op2(op2)
		{}
		virtual T operator()(const T& t)
		{
			return m_base*(*m_op2)(t);
		}
		virtual std::string toString() const
		{
			std::stringstream ss;
			ss << m_base << strings::sMult << "(" << m_op2->toString() << ")";
			return ss.str();
		}
		virtual tFunctionPtr derivative(int num)
		{
			return tFunctionPtr(new multiply<T, eMultNum>(m_base, m_op2->derivative(num)));
		}

		virtual void optimize()
		{
			// TODO - Nothing is here for now
		}

		virtual tFunctionPtr clone()
		{
			return tFunctionPtr(new multiply<T, eMultNum>(m_base, m_op2));
		}
	}; // class multiply

}

#endif // __MULTIPLY__H_
