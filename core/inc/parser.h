#ifndef __PARSER_H_
#define __PARSER_H_

#include "function.h"
#include <memory>
#include <string>
#include <stack>
#include <vector>
#include <iterator>
#include <iostream>

namespace core
{
	//extern const char sLeftPar;
	//extern const char sRightPar;
	//extern const char sPlus;
	//extern const char sMinus;
	//extern const char sMult;
	//extern const char sDiv;
	//extern const char sPower;
	//extern const char sParam;
	//extern const std::string sExp;
	const char sLeftPar 	= '(';
	const char sRightPar 	= ')';
	const char sPlus 		= '+';
	const char sMinus		= '-';
	const char sMult		= '*';
	const char sDiv			= '/';
	const char sPower		= '^';
	const char sParam		= 't';
	const std::string sExp	= "exp";

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
			std::cout << __func__ << std::endl;
			std::cout << "expression = " << expression << std::endl;
			std::shared_ptr<function<T> > ret(nullptr);
			std::vector<std::string> post;
			post.reserve(expression.size());
			get_postfix(expression, std::back_inserter(post));
			std::copy(post.begin(), post.end(), std::ostream_iterator<std::string>(std::cout, ""));
			std::cout << std::endl;
			return ret;
		}
	protected:
		bool is_operator(char op) const
		{
			switch (op)
			{
			case sPlus:
			case sMinus:
			case sDiv:
			case sMult:
			case sPower:
				return true;
			default:
				return false;
			}
			return false;
		}

		bool is_operator(const std::string& op) const
		{
			return (op.size() != 1 ? false : is_operator(op[0]));
		}
		
		bool is_parethesis(char c) const
		{
			return (c == sLeftPar) || (c == sRightPar);
		}

		bool is_parameter(char c) const
		{
			return c == sParam;
		}

		bool is_numeric(char c) const
		{
			return (c >= '0' && c <= '9') || (c == '.');
		}

		int precedence(char c) const
		{
			switch (c)
			{
			case sPower:
				return 4;
			case sMult:
			case sDiv:
				return 3;
			case sPlus:
			case sMinus:
				return 2;
			default:
				return 0;
			}
			return 0;
		}

		bool check_token(const std::string& t) const
		{
			return t == sExp;
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
				std::cout << "Processing " << *it << std::endl;
				if (is_operator(*it))
				{
					std::cout << "Was operator\n";
					if (!number.empty())
					{
						std::cout << "Number was not empty, adding " << number << std::endl;
						*output++ = number;
						number.clear();
					}
					param.clear();
					if (is_start)
					{
						ops.push(std::string(2, *it));
						is_start = false;
						continue;
					}
					is_start = false;
					int pit = precedence(*it);
					while (!ops.empty())
					{
						std::string op = ops.top();
						std::cout << "In the stack checking " << op << std::endl;
						if (is_operator(op) && (precedence(op[0]) >= pit))
						{
							*output++ = op;
							ops.pop();
							std::cout << op << " added, now poping" << std::endl;
						}
						else
							break;
					}
					ops.push(std::string(1, *it));
				}
				else if (is_parameter(*it))
				{
					is_start = false;
					std::cout << "Was parameter\n";
					if (!number.empty() || !param.empty() || !token.empty())
						throw parserException("Invalid number or parameter");
					param = *it;
					*output++ = param;
					param.clear();
				}
				else if (is_numeric(*it))
				{
					is_start = false;
					std::cout << "Was numeric\n";
					if (!param.empty() || !token.empty())
						throw parserException("Invalid number or parameter");
					number += *it;
				}
				else if (is_parethesis(*it))
				{
					is_start = true;
					std::cout << "Was parenthesis\n";
					if (!number.empty())
					{
						std::cout << "Number was not empty, adding " << number << std::endl;
						*output++ = number;
						number.clear();
					}
					param.clear();
					if (*it == sLeftPar)
					{
						std::cout << "Was left\n";
						if (!token.empty() && check_token(token))
							ops.push(token);
						token.clear();
						std::cout << "Added to the stack\n";
						ops.push(std::string(1, *it));
						continue;
					}
					bool found = false;
					while (!ops.empty())
					{
						std::string op = ops.top();
						std::cout << "In the stack checking " << op << std::endl;
						if (op[0] == sLeftPar)
						{
							std::cout << "Found the left par\n";
							found = true;
							ops.pop();
							break;
						}
						std::cout << "Not found the left par yet, so popping the ops\n";
						*output++ = op;
						ops.pop();
					}
					if (!found)
						throw parserException("Closing parenthesis without opening one");
					if (!ops.empty() && !is_operator(ops.top()) && !is_parethesis(ops.top()[0]))
					{
						std::cout << "Found a function - " << ops.top() << std::endl;
						*output++ = ops.top();
						ops.pop();
					}
					is_start = false;
				}
				else
				{
					if (!number.empty() || !param.empty())
						throw parserException("Invalid number or parameter");
					std::cout << "This is a token, accumulate " << *it << std::endl;
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
