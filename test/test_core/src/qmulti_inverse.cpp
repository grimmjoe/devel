
#include "arithmetic.h"
#include "function.h"
#include "add.h"
#include "comparator.h"
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

template <class T>
class garsevan : public core::function<T>
{
	T begin;
	T end;
	tFunctionPtr m_base;
public:
	garsevan(T b, T e, tFunctionPtr mb)
		: begin(b)
		, end(e)
		, m_base(mb)
	{
	}

	virtual T operator()(const T& t)
	{
		if ((t > end) || (t < begin))
			return 0;
		if ((t == begin) || (t == end))
			return 0.5*(*m_base)(t);

		return (*m_base)(t);
	}
	virtual std::string toString() const
	{
		return "";
	}
	virtual tFunctionPtr derivative(int num)
	{
		return tFunctionPtr(nullptr);
	}
	virtual void optimize()
	{
	}
	virtual tFunctionPtr clone()
	{
		return tFunctionPtr(new garsevan<T>(begin, end, m_base));
	}
};

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
	double a;
	double b;
	double h;
	double epsilon;
	fin >> k >> a >> b >> h  >> epsilon;
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
	tDiffInfo di(1, 1, k);
	tApp app(di);
	if (epsilon != 0)
	{
		app.setComparator(std::shared_ptr<core::comparator<double> >(new core::comparator<double>(epsilon)));
	}
	tMatrix realOrig(n, m, tFunctionPtr(nullptr));
	for (double begin = a; begin <= b; begin += h)
	{
		std::cout << "********************************\n";
		std::cout << "[" << begin << ", " << begin+h << "]\n";
		std::cout << "********************************\n";
		tDiffInfo di(begin, h, k);
		app.setDiffInfo(di);
		tDiscretes theDiscretes;
		app.getQInverse(theMatrix, theDiscretes);
		tMatrix origMatrix(n, m, tFunctionPtr(nullptr));
		app.restoreTaylorSingle(theDiscretes, origMatrix, k-1);
		if (realOrig[1][1] == nullptr)
		{
			for (int i = 1; i <= n; ++i)
				for (int j = 1; j <= m; ++j)
					realOrig[i][j] = tFunctionPtr(new garsevan<double>(begin-h, begin+h, origMatrix[i][j]));
		}
		else
		{
			for (int i = 1; i <= n; ++i)
				for (int j = 1; j <= m; ++j)
					realOrig[i][j] = tFunctionPtr(new core::add<double>(realOrig[i][j],
											tFunctionPtr(new garsevan<double>(begin, begin+h, origMatrix[i][j]))));
		}
	}
	std::cout << "Multi-point restored:\n";
	std::ofstream qmulti("qmulti_out.txt");
	assert (qmulti.good());
	qmulti << "t";
	for (int i = 1; i <= n; ++i)
		for (int j = 1; j <= m; ++j)
			qmulti << " qinv[" << i << ", " << j << "]";
	qmulti << std::endl;
	double delta = 0.1;
	for (double begin = a-1; begin <= b+1; begin += delta)
	{
		qmulti << begin;
		for (int i = 1; i <= n; ++i)
			for (int j = 1; j <= m; ++j)
				qmulti << " " << (*realOrig[i][j])(begin);
		qmulti << std::endl;
	}
	qmulti.close();
	//tDiscretes theDiscretes;
	//tDiscretes theMatrixDiscretes;
	//app.applyDiffTrans(theMatrix, theMatrixDiscretes);
	////app.getQInverse(theMatrixDiscretes, theDiscretes);
	//app.getQInverse(theMatrix, theDiscretes);
	//std::cout << "\nDiscretes:\n";
	//std::copy(theDiscretes.begin(), theDiscretes.end(), std::ostream_iterator<tDiscrete>(std::cout, "\n"));
	//std::cout << "\nRestoring the original\n";
	//tMatrix origMatrix(n, m, tFunctionPtr(nullptr));
	//app.restoreTaylorSingle(theDiscretes, origMatrix, k-1);
	//std::cout << "The original:\n";
	//printFuncMatrix(origMatrix);

	//if (app.checkB_Q_BQ_Inverse(theMatrixDiscretes, theDiscretes, di.K>0? di.K-1:0))
	//	std::cout << "The calculation is correct" << std::endl;
	//else
	//	std::cout << "The calculation is NOT CORRECT" << std::endl;

	return 0;
}
