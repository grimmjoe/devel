
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

typedef core::apparatus<double> tApp;
typedef typename tApp::tMatrixDiscrete tDiscrete;
typedef typename tApp::tMatrixDiscretes tDiscretes;
typedef typename tApp::tFunction tFunction;
typedef typename tApp::tFunctionPtr tFunctionPtr;
typedef typename tApp::tFuncMatrix tMatrix;
typedef typename tApp::tDiffInfo tDiffInfo;

void printFuncMatrix(const tMatrix& mat)
{
	for (int i = 1; i <= mat.getNumRows(); ++i)
	{
		for (int j = 1; j <= mat.getNumCols(); ++j)
		{
			std::cout << mat[i][j]->toString() << " ";
		}
		std::cout << std::endl;
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
	int m, n;
	fin >> m >> n;
	int k;
	int tv;
	fin >> k >> tv;
	tMatrix theMatrix(m, n, tFunctionPtr(nullptr));
	core::parser<double> p;
	for (int i = 1; i <= m; ++i)
	{
		for (int j = 1; j <= n; ++j)
		{
			std::string s;
			fin >> s;
			theMatrix[i][j] = p.parse(s);
		}
	}
	std::cout << "The matrix:\n";
	printFuncMatrix(theMatrix);
	tMatrix B(n, m, tFunctionPtr(nullptr));
	for (int i = 1; i <= n; ++i)
	{
		for (int j = 1; j <= m; ++j)
		{
			std::string s;
			fin >> s;
			B[i][j] = p.parse(s);
		}
	}
	std::cout << "B:\n";
	printFuncMatrix(B);
	tDiffInfo di(tv, 1, k);
	tApp app(di);
	tDiscretes theDiscretes;
	tDiscretes theMatrixDiscretes;
	app.applyDiffTrans(theMatrix, theMatrixDiscretes);
	tDiscretes theB;
	app.applyDiffTrans(B, theB);
	app.getBInverse(theMatrixDiscretes, m, theB, theDiscretes);
	std::cout << "\nDiscretes:\n";
	std::copy(theDiscretes.begin(), theDiscretes.end(), std::ostream_iterator<tDiscrete>(std::cout, "\n"));
	std::cout << "\nRestoring the original\n";
	tMatrix origMatrix(n, m, tFunctionPtr(nullptr));
	app.restoreTaylorSingle(theDiscretes, origMatrix);
	std::cout << "The original:\n";
	printFuncMatrix(origMatrix);

	if (app.checkB_Q_BQ_Inverse(theMatrixDiscretes, theDiscretes))
		std::cout << "The calculation is correct" << std::endl;
	else
		std::cout << "The calculation is NOT CORRECT" << std::endl;

	return 0;
}
