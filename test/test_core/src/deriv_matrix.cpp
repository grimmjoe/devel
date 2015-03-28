
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
	double tv;
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
	//std::cout << "Calculating the A^3...\n";
	//tMatrix theM = theMatrix * theMatrix;
	////std::cout << "A^2:\n";
	////theMatrix = theM;
	////printFuncMatrix(theMatrix);
	//theMatrix = theM * theMatrix;
	//std::cout << "A^3:\n";
	//printFuncMatrix(theMatrix);
	//typedef core::matrix<double> tDiscrete;
	//std::vector<tDiscrete> theDiscretes(k+1, tDiscrete(m, n, 0));
	tDiffInfo di(tv, 1, k);
	tDiscretes theDiscretes;
	tApp app(di);
	app.applyDiffTrans(theMatrix, theDiscretes);
	//for (int K = 0; K <= k; ++K)
	//{
	//	for (int i = 1; i <= m; ++i)
	//	{
	//		for (int j = 1; j <= n; ++j)
	//		{
	//			theMatrix[i][j] = theMatrix[i][j]->derivative(K > 0 ? 1 : 0);
	//			theDiscretes[K][i][j] = (*theMatrix[i][j])(tv);
	//		}
	//	}
	//	std::cout << "K: = " << K << std::endl;
	//	printFuncMatrix(theMatrix);
	//}
	std::cout << "\nDiscretes:\n";
	std::copy(theDiscretes.begin(), theDiscretes.end(), std::ostream_iterator<tDiscrete>(std::cout, "\n"));
	std::cout << "\nRestoring the original\n";
	tMatrix origMatrix(m, n, tFunctionPtr(nullptr));
	app.restoreTaylorSingle(theDiscretes, origMatrix);
	std::cout << "The original:\n";
	printFuncMatrix(origMatrix);

	return 0;
}
