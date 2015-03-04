#ifndef __NEGATE_FUNCTION_H_
#define __NEGATE_FUNCTION_H_

#include "function.h"
#include "strings.h"
#include "constant_function.h"
#include <sstream>

namespace core
{
	template <class T>
	class negate : public function<T>
	{
	public:
		typedef function<T> tBase;
		using typename tBase::tFunctionPtr;
	private:
		tFunctionPtr m_base;
	public:
		negate(tFunctionPtr b)
			: m_base(b)
		{}
		virtual T operator()(const T& t)
		{
			return -(*m_base)(t);
		}
		virtual std::string toString() const
		{
			std::stringstream ss;
			ss << "-" << "(" << m_base->toString() << ")";
			return ss.str();
		}
		virtual tFunctionPtr derivative(int num)
		{
			return tFunctionPtr(new negate<T>(m_base->derivative(num)));
		}

		virtual void optimize()
		{
			// TODO - Nothing is here for now
		}

		virtual tFunctionPtr clone()
		{
			return tFunctionPtr(new negate<T>(m_base));
		}
	};
}

#endif // __NEGATE_FUNCTION_H_
