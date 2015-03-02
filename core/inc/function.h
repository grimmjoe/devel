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
		virtual T operator()(const T& t) = 0;
		virtual std::string toString() const = 0;
		virtual std::shared_ptr<function<T> > derivative(int num) = 0;
		virtual void optimize() = 0;
		virtual std::shared_ptr<function<T> > clone() = 0;
	};
}

#endif // __FUNCTION_H__
