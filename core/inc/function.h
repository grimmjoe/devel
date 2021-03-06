#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#include <string>
#include <memory>

namespace core
{
	template <class T>
	class function
	{
	public:
		typedef function<T> tFunction;
		typedef std::shared_ptr<tFunction> tFunctionPtr;

		virtual T operator()(const T& t) = 0;
		virtual std::string toString() const = 0;
		virtual tFunctionPtr derivative(int num) = 0;
		virtual void optimize() = 0;
		virtual tFunctionPtr clone() = 0;

		virtual ~function()
		{}
	};
}

#endif // __FUNCTION_H__
