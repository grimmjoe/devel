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
	if (argc != 2)
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
	return 0;
}
