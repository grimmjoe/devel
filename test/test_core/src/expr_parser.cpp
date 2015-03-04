#include "parser.h"
#include "function.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>

using namespace core;

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cerr << "Wrong num of args, please specify the input expression\n";
		return 1;
	}
	parser<double> p;
	std::shared_ptr<function<double> > expression = p.parse(argv[1]);
	if (expression == nullptr)
		std::cerr << "Parse error\n";
	else
		std::cout << "The expression is: " << expression->toString() << std::endl;
	if (argc >= 3)
	{
		std::stringstream ss;
		ss << argv[2];
		int n;
		ss >> n;
		for (int i = 0; i <= n; ++i)
		{
			typename core::function<double>::tFunctionPtr ptr = expression->derivative(i);
			std::cout << "i: " << i << ", derivative: " << (ptr == nullptr ? "NULL" : ptr->toString()) << std::endl;
		}
	}
	return 0;
}
