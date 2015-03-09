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

void inverse(tApp& app, const tDiscretes& c, int m, int n, int r, int K, tDiscretes& c1)
{
	tDiscretes L;
	tDiscretes U;

	//app.getLU(c, r, L, U);
	app.getInverse(c, r, c1);

	//tDiscrete L(K, tMatrix(r, n, 0));
	//tDiscrete U(K, tMatrix(m, r, 0));
	//// First
	//for (int k = 0; k < K; ++k)
	//	for (int i = 1; i <= r; ++i)
	//		L[k][i][1] = c[k][i][1];
	//
	//for (int k = 0; k < K; ++k)
	//{
	//	for (int j = 1; j <= r ; ++j)
	//	{
	//		double sum = 0;
	//		for (int l = 1; l <= k; ++l)
	//		{
	//			sum += U[k-l][1][j]*L[l][1][1];
	//		}
	//		U[k][1][j] = (c[k][1][j]-sum)/L[0][1][1];
	//	}
	//}

	//// Next
	//for (int k = 0; k < K; ++k)
	//{
	//	for (int i = 2; i <= r; ++i)
	//	{
	//		for (int p = 2; p <= i; ++p)
	//		{
	//			double sum = 0;
	//			for (int j = 1; j <= p-1; ++j)
	//			{
	//				for (int l = 0; l <= k; ++l)
	//				{
	//					sum += L[l][i][j]*U[k-l][j][p];
	//				}
	//			}
	//			L[k][i][p] = c[k][i][p] - sum;
	//		}
	//	}
	//}

	//for (int k = 0; k < K; ++k)
	//{
	//	for (int i = 2; i <= r; ++i)
	//	{
	//		for (int p = i; p <= r; ++p)
	//		{
	//			double sum1 = 0;
	//			for (int j = 1; j <= i-1; ++j)
	//				for (int g = 0; g <= k; ++g)
	//					sum1 += L[g][i][j]*U[k-g][j][p];

	//			double sum2 = 0;
	//			for (int l = 1; l <= k; ++l)
	//			{
	//				sum2 += U[k-l][i][p]*L[l][i][i];
	//			}
	//			U[k][i][p] = (c[k][i][p] - sum1 - sum2)/L[0][i][i];
	//		}
	//	}
	//}
	//std::cout << "L:\n";
	//std::copy(L.begin(), L.end(), std::ostream_iterator<tMatrix>(std::cout, "\n"));
	//std::cout << "U:\n";
	//std::copy(U.begin(), U.end(), std::ostream_iterator<tMatrix>(std::cout, "\n"));

	// Now the inverses
	//tDiscretes L1(K, tDiscrete(n, r, 0));
	//tDiscretes U1(K, tDiscrete(r, m, 0));
	//tDiscretes e(K, tDiscrete(m, n, 0));
	//for (int i = 1; i <= m; ++i)
	//	e[0][i][i] = 1;
	//
	//// Now L1
	//for (int k = 0; k < K; ++k)
	//{
	//	for (int i = 1; i <= r; ++i)
	//	{
	//		for (int j = 1; j <= r; ++j)
	//		{
	//			double sum1 = 0;
	//			for (int p = 1; p <= i-1; ++p)
	//				for (int l = 0; l <= k; ++l)
	//					sum1 += L[l][i][p]*L1[k-l][p][j];

	//			double sum2 = 0;
	//			for (int l = 1; l <= k; ++l)
	//				sum2 += L1[k-l][i][j]*L[l][i][i];
	//			L1[k][i][j] = (e[k][i][j] - sum1 - sum2)/L[0][i][i];
	//		}
	//	}
	//}

	//// Now U1
	//for (int k = 0; k < K; ++k)
	//{
	//	for (int i = r; i >= 1; --i)
	//	{
	//		for (int j = 1; j <= r; ++j)
	//		{
	//			double sum = 0;
	//			for (int p = i+1; p <= m; ++p)
	//				for (int g = 0; g <= k; ++g)
	//					sum += U[g][i][p]*U1[k-g][p][j];
	//			U1[k][i][j] = e[k][i][j] - sum;
	//		}
	//	}
	//}
	////std::cout << "L1:\n";
	////std::copy(L1.begin(), L1.end(), std::ostream_iterator<tDiscrete>(std::cout, "\n"));
	////std::cout << "U1:\n";
	////std::copy(U1.begin(), U1.end(), std::ostream_iterator<tDiscrete>(std::cout, "\n"));
	//for (int k = 0; k < K; ++k)
	//{
	//	for (int l = 0; l <= k; ++l)
	//		c1[k] += U1[l]*L1[k-l];
	//}
}

int main(int argc, char* argv[])
{
	int m = 3;
	int n = 3;
	int r = 3;
	int K = 5;
	if (argc != 2)
	{
		std::cerr << "Wrong num of args, please specify the input file\n";
		return 1;
	}
	std::ifstream fstr(argv[1]);
	fstr >> m >> n >> r >> K;
	tDiscretes c(K, tDiscrete(m, n, 0));
	for (int k = 0; k < K; ++k)
	{
		for (int i = 1; i <= m; ++i)
			for (int j = 1; j <= n; ++j)
				fstr >> c[k][i][j];
	}
	std::cout << "Original:\n";
	std::copy(c.begin(), c.end(), std::ostream_iterator<tDiscrete>(std::cout, "\n"));
	tDiscretes c1(K, tDiscrete(m, n, 0));
	tDiffInfo di(0, 1, K);
	tApp app(di);
	inverse(app, c, m, n, r, K, c1);
	std::cout << "Inverse:\n";
	std::copy(c1.begin(), c1.end(), std::ostream_iterator<tDiscrete>(std::cout, "\n"));
	return 0;
}
