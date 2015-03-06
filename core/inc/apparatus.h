#ifndef __APPARATUS_H__
#define __APPARATUS_H__

#include "function.h"
#include "matrix.h"
#include <vector>
#include <cmath>

namespace core
{
	template <class T>
	class apparatus
	{
	public:
		typedef T tArgType;
		typedef typename function<T>::tFunction tFunction;
		typedef typename function<T>::tFunctionPtr tFunctionPtr;
		typedef core::matrix<tFunctionPtr> tFuncMatrix;
		typedef core::matrix<T> tMatrixDiscrete;
		typedef std::vector<tMatrixDiscrete> tMatrixDiscretes;
	
		struct diffInfo
		{
			T tv;
			T H;
			int K;
			diffInfo(T t, T h, T k)
				: tv(t)
				, H(h)
				, K(k)
			{}
		};
	
		typedef diffInfo tDiffInfo;
	
		//
		/// Apply the Taylor-diff transformation @a di on the @a sourceMatrix into the discretes @a discs
		//
		bool applyDiffTrans(const tFuncMatrix& sourceMatrix, const tDiffInfo& di, tMatrixDiscretes& discs);
	
		//
		/// Get LU decomposition of the matrix
		//
		bool getLU(const tFuncMatrix& sourceMatrix, const diffInfo& di, tMatrixDiscretes& l, tMatrixDiscretes& u);
	};
} // namespace core

template <class T>
bool core::apparatus<T>::getLU(const tFuncMatrix& sourceMatrix, const diffInfo& di, tMatrixDiscretes& l, tMatrixDiscretes& u)
{
	const int m = sourceMatrix.getNumRows();
	const int n = sourceMatrix.getNumCols();
	return true;
}

template <class T>
bool core::apparatus<T>::applyDiffTrans(const tFuncMatrix& sourceMatrix, const tDiffInfo& di, tMatrixDiscretes& discs)
{
	const int m = sourceMatrix.getNumRows();
	const int n = sourceMatrix.getNumCols();
	discs.resize(di.K+1, tMatrixDiscrete(m, n, 0));
	tFuncMatrix theMatrix = sourceMatrix;
	long long kfac = 1;
	for (int k = 0; k <= di.K; ++k)
	{
		kfac *= (k > 0 ? k : 1);
		for (int i = 1; i <= m; ++i)
		{
			for (int j = 1; j <= n; ++j)
			{
				theMatrix[i][j] = theMatrix[i][j]->derivative(k > 0 ? 1 : 0);
				discs[k][i][j] = std::pow(di.H, k)*(*theMatrix[i][j])(di.tv)/(T)kfac;
			}
		}
	}
	return true;
}


#endif // __APPARATUS_H__
