#ifndef __PARSER_H_
#define __PARSER_H_

#include "function.h"
#include "power_function.h"
#include "exp_function.h"
#include "constant_function.h"
#include "parameter_function.h"
#include "negate_function.h"
#include "arithmetic_functions.h"
#include "strings.h"
#include <memory>
#include <string>
#include <stack>
#include <vector>
#include <iterator>
#include <iostream>
#include <cassert>

namespace core
{
	class parserException : public std::exception
	{
		std::string err;
	public:
		parserException(const std::string& e)
			: err(e)
		{}
		virtual const char* what() const throw()
		{
			return err.c_str();
		}
	}; // parserException

	template <class T>
	class parser
	{
	public:
		parser()
		{}
		std::shared_ptr<function<T> > parse(const std::string& expression)
		{
			typedef typename function<T>::tFunctionPtr tFunctionPtr;
			//std::cout << __func__ << std::endl;
			//std::cout << "expression = " << expression << std::endl;
			std::vector<std::string> post;
			post.reserve(expression.size());
			get_postfix(expression, std::back_inserter(post));
			//std::copy(post.begin(), post.end(), std::ostream_iterator<std::string>(std::cout, ""));
			//std::cout << std::endl;
			std::stack<std::shared_ptr<function<T> > > theStack;
			for (auto it = post.begin(); it != post.end(); ++it)
			{
				std::string val = *it;
				assert (!val.empty());
				//std::cout << "Processing " << val << std::endl;

				if (is_operator(val[0]))
				{
					if (val.size() > 1 && val[0] == strings::sMinus)
					{
						//std::cout << "Need to negate\n";
						tFunctionPtr f = theStack.top();
						theStack.pop();
						theStack.push(tFunctionPtr(new negate<T>(f)));
						//std::cout << "Pushed the new function\n";
						continue;
					}
					switch (val[0])
					{
					case strings::sPlus:
					{
						tFunctionPtr f1 = theStack.top();
						theStack.pop();
						tFunctionPtr f2 = theStack.top();
						theStack.pop();
						theStack.push(tFunctionPtr(new add<T>(f2, f1)));
					}
					break;
					case strings::sMinus:
					{
						tFunctionPtr f1 = theStack.top();
						theStack.pop();
						tFunctionPtr f2 = theStack.top();
						theStack.pop();
						theStack.push(tFunctionPtr(new subtract<T>(f2, f1)));
					}
					break;
					case strings::sMult:
					{
						tFunctionPtr f1 = theStack.top();
						theStack.pop();
						tFunctionPtr f2 = theStack.top();
						theStack.pop();
						theStack.push(tFunctionPtr(new multiply<T>(f2, f1)));
					}
					break;
					case strings::sDiv:
					{
						tFunctionPtr f1 = theStack.top();
						theStack.pop();
						tFunctionPtr f2 = theStack.top();
						theStack.pop();
						theStack.push(tFunctionPtr(new divide<T>(f2, f1)));
					}
					break;
					case strings::sPower:
					{
						tFunctionPtr f1 = theStack.top();
						theStack.pop();
						tFunctionPtr f2 = theStack.top();
						theStack.pop();
						const std::string f1_str = f1->toString();
						std::stringstream ss;
						ss << f1_str;
						T p;
						ss >> p;
						//std::cout << "power = " << f1_str << std::endl;
						//std::cout << "power in double = " << p << std::endl;
						theStack.push(tFunctionPtr(new power<T>(f2, p)));
					}
					break;
					}
					continue;
				}
				if (is_parameter(val))
				{
					theStack.push(tFunctionPtr(new parameter<T>()));
					continue;
				}
				if (check_token(val))
				{
					if (val == EXPONENT)
					{
						tFunctionPtr f1 = theStack.top();
						theStack.pop();
						theStack.push(tFunctionPtr(new exp<T>(std::exp(1), f1)));
					}
					else
					{
						assert (! "Unknown function");
					}
					continue;
				}
				std::stringstream ss;
				T v;
				ss << val;
				ss >> v;
				//std::cout << "Const = " << v << std::endl;
				theStack.push(tFunctionPtr(new const_function<T>(v)));
			}
			return (theStack.empty() ? tFunctionPtr(nullptr) : theStack.top());
		}
	protected:
		bool is_operator(char op) const
		{
			switch (op)
			{
			case strings::sPlus:
			case strings::sMinus:
			case strings::sDiv:
			case strings::sMult:
			case strings::sPower:
				return true;
			default:
				return false;
			}
			return false;
		}

		bool is_operator(const std::string& op) const
		{
			if (op.size() > 1 && op[0] == strings::sMinus)
				return true;
			return (op.size() != 1 ? false : is_operator(op[0]));
		}
		
		bool is_parethesis(char c) const
		{
			return (c == strings::sLeftPar) || (c == strings::sRightPar);
		}

		bool is_parameter(const std::string& op)
		{
			return (op.size() != 1 ? false : is_parameter(op[0]));
		}

		bool is_parameter(char c) const
		{
			return c == strings::sParam;
		}

		bool is_numeric(char c) const
		{
			return (c >= '0' && c <= '9') || (c == '.');
		}

		int precedence(const std::string& s) const
		{
			if (s.size() > 1 && s[0] == strings::sMinus)
				return precedence(strings::sDiv);
			return precedence(s[0]);
		}

		int precedence(char c) const
		{
			switch (c)
			{
			case strings::sPower:
				return 4;
			case strings::sMult:
			case strings::sDiv:
				return 3;
			case strings::sPlus:
			case strings::sMinus:
				return 2;
			default:
				return 0;
			}
			return 0;
		}

		bool check_token(const std::string& t) const
		{
			return t == EXPONENT;
		}

		template <class IT>
		void get_postfix(const std::string& expression, IT output) const
		{
			std::stack<std::string> ops;
			std::string token = "";
			std::string number = "";
			std::string param = "";
			bool is_start = true;
			for (auto it = expression.begin(); it != expression.end(); ++it)
			{
				//std::cout << "Processing " << *it << std::endl;
				if (is_operator(*it))
				{
					//std::cout << "Was operator\n";
					if (!number.empty())
					{
						//std::cout << "Number was not empty, adding " << number << std::endl;
						*output++ = number;
						number.clear();
					}
					param.clear();
					std::string current(1, *it);
					if (is_start && (current[0] == strings::sMinus))
					{
						current.push_back(strings::sMinus);
						is_start = false;
					}
					is_start = false;
					int pit = precedence(current[0]);
					if (current.size() > 1)
						pit = precedence(strings::sDiv);
					while (!ops.empty())
					{
						std::string op = ops.top();
						//std::cout << "In the stack checking " << op << std::endl;
						if (is_operator(op) && (precedence(op) >= pit))
						{
							*output++ = op;
							ops.pop();
							//std::cout << op << " added, now poping" << std::endl;
						}
						else
							break;
					}
					ops.push(current);
				}
				else if (is_parameter(*it))
				{
					is_start = false;
					//std::cout << "Was parameter\n";
					if (!number.empty() || !param.empty() || !token.empty())
						throw parserException("Invalid number or parameter");
					param = *it;
					*output++ = param;
					param.clear();
				}
				else if (is_numeric(*it))
				{
					is_start = false;
					//std::cout << "Was numeric\n";
					if (!param.empty() || !token.empty())
						throw parserException("Invalid number or parameter");
					number += *it;
				}
				else if (is_parethesis(*it))
				{
					is_start = true;
					//std::cout << "Was parenthesis\n";
					if (!number.empty())
					{
						//std::cout << "Number was not empty, adding " << number << std::endl;
						*output++ = number;
						number.clear();
					}
					param.clear();
					if (*it == strings::sLeftPar)
					{
						//std::cout << "Was left\n";
						if (!token.empty() && check_token(token))
							ops.push(token);
						token.clear();
						//std::cout << "Added to the stack\n";
						ops.push(std::string(1, *it));
						continue;
					}
					bool found = false;
					while (!ops.empty())
					{
						std::string op = ops.top();
						//std::cout << "In the stack checking " << op << std::endl;
						if (op[0] == strings::sLeftPar)
						{
							//std::cout << "Found the left par\n";
							found = true;
							ops.pop();
							break;
						}
						//std::cout << "Not found the left par yet, so popping the ops\n";
						*output++ = op;
						ops.pop();
					}
					if (!found)
						throw parserException("Closing parenthesis without opening one");
					if (!ops.empty() && !is_operator(ops.top()) && !is_parethesis(ops.top()[0]))
					{
						//std::cout << "Found a function - " << ops.top() << std::endl;
						*output++ = ops.top();
						ops.pop();
					}
					is_start = false;
				}
				else
				{
					if (!number.empty() || !param.empty())
						throw parserException("Invalid number or parameter");
					//std::cout << "This is a token, accumulate " << *it << std::endl;
					token += *it;
				}
			}
			if (!number.empty())
				*output++ = number;
			if (!param.empty())
				*output++ = param;
			if (!token.empty())
				throw parserException("Function without arguments detected!");
			while (!ops.empty())
			{
				std::string op = ops.top();
				if (is_parethesis(op[0]))
					throw parserException("Parenthesis mismatch");
				*output++ = op;
				ops.pop();
			}
		}
	}; // class parser
} // namespace core

#endif // __PARSER_H_
