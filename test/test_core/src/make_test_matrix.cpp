
#include "arithmetic.h"
#include "function.h"
#include "matrix.h"
#include "parser.h"
#include "apparatus.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <iterator>

#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cstdio>

#include <chrono>

typedef core::apparatus<double> tApp;
typedef typename tApp::tMatrixDiscrete tDiscrete;
typedef typename tApp::tMatrixDiscretes tDiscretes;
typedef typename tApp::tFunction tFunction;
typedef typename tApp::tFunctionPtr tFunctionPtr;
typedef typename tApp::tFuncMatrix tMatrix;
typedef typename tApp::tDiffInfo tDiffInfo;

struct testData
{
	int number;
	std::string path;
	std::string prefix;

	int m;
	int n;
	int r;
	int deg;
	double sp1;
	double sp2;
	int k;
	double tv;
	double epsilon;
};

template <class T>
std::string makeStr(const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}

void makeTestMatrices(testData& data)
{
	std::cout << "Making:\n";
	std::cout << "number = " << data.number << std::endl;
	std::cout << "m = " << data.m << std::endl;
	std::cout << "n = " << data.n << std::endl;
	std::cout << "r = " << data.r << std::endl;
	std::cout << "deg = " << data.deg << std::endl;
	std::cout << "sp1 = " << data.sp1 << std::endl;
	std::cout << "sp2 = " << data.sp2 << std::endl;
	std::cout << "k = " << data.k << std::endl;
	std::cout << "tv = " << data.tv << std::endl;
	std::cout << "epsilon = " << data.epsilon << std::endl;

	int numOfZeros = data.m*data.n*(1-data.sp1);
	int numOfNonZeroes = data.m*data.n-numOfZeros;
	int numOfNZCoeff = data.sp2*data.deg;
	std::cout << "numOfZeros = " << numOfZeros << std::endl;
	std::cout << "numOfNZCoeff = " << numOfNZCoeff << std::endl;
	srand(time(NULL));
	for (int n = 1; n <= data.number; ++n)
	{
		std::stringstream ss;
		ss << data.path << "/" << data.prefix;
		ss << std::setfill('0') << std::setw(3) << n;
		std::ofstream out(ss.str());
		assert(out.good());
		out << data.m << " " << data.n << " " << data.k << " " << data.tv << " " << data.epsilon << std::endl;
		int nz = 0;
		std::vector<std::vector<std::string> > theMatrix (data.m, std::vector<std::string>(data.n, " "));
		int r = data.r;
		int m = data.m;
		for (int j = 1; j <= r; ++j)
		{
			for (int i = 1; i <= r; ++i)
			{
				if (nz > numOfNonZeroes)
				{
					theMatrix[i-1][j-1] = "0";
					continue;
				}
				++nz;
				std::string val = "";
				for (int p = 1; p <= numOfNZCoeff; ++p)
				{
					int d = rand() % (data.deg+1);
					if (!val.empty())
						val += "+";
					if (d == 0)
					{
						int rnum = rand()%1000+1;
						val += makeStr(rnum);
						continue;
					}
					int coeff = rand()%1000+1;
					if (coeff != 1)
					{
						val += makeStr(coeff);
						val += "*";
					}
					val += "t";
					if (d != 1)
					{
						val += "^";
						val += makeStr(d);
					}
				}
				theMatrix[i-1][j-1] = (val.empty() ? "0" : val);
			}
		}

		for (int j = r+1; j <= data.n; ++j)
		{
			for (int i = 1; i <= r; ++i)
			{
				theMatrix[i-1][j-1] = theMatrix[i-1][r-1];
			}
		}
		for (int i = r+1; i <= m; ++i)
		{
			for (int j = 1; j <= data.n; ++j)
				theMatrix[i-1][j-1] = theMatrix[r-1][j-1];
		}
		for (int i = 1; i <= m; ++i)
		{
			for (int j = 1; j < data.n; ++j)
			{
				out << theMatrix[i-1][j-1] << " ";
			}
			out << theMatrix[i-1][n-1] << std::endl;
		}
	}
}


int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		std::cerr << "Please specify the file\n";
		return 1;
	}
	std::ifstream fin(argv[1]);
	assert (fin.good());
	testData data;
	fin >> data.number;
	fin >> data.m >> data.n >> data.r;
	fin >> data.deg >> data.sp1 >> data.sp2;
	fin >> data.k >> data.tv >> data.epsilon;
	fin >> data.path;
	fin >> data.prefix;

	makeTestMatrices(data);

	return 0;
}
