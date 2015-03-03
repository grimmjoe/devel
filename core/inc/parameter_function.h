#ifndef __PARAMETER_FUNCTION_H_
#define __PARAMETER_FUNCTION_H_

#include "function.h"
#include "strings.h"
#include "constant_function.h"
#include <sstream>

namespace core
{
	template <class T>
	class parameter : public function<T>
	{
	public:
		typedef function<T> tBase;
		using typename tBase::tFunctionPtr;
	private:
	public:
		parameter()
		{}
		virtual T operator()(const T& t)
		{
			return t;
		}
		virtual std::string toString() const
		{
			std::stringstream ss;
			ss << strings::sParam;
			return ss.str();
		}
		virtual tFunctionPtr derivative(int num)
		{
			if (num > 0)
				return tFunctionPtr(new const_function<T>(0));
			else if (num == 0)
				return this->clone();
			else
				return tFunctionPtr(nullptr);
		}

		virtual void optimize()
		{
			// TODO - Nothing is here for now
		}

		virtual tFunctionPtr clone()
		{
			return tFunctionPtr(new parameter<T>());
		}
	};
}

#endif // __PARAMETER_FUNCTION_H_
