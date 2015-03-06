#ifndef __APPARATUS_H__
#define __APPARATUS_H__

#include "function.h"
#include "constant_function.h"
#include "power_function.h"
#include "add.h"
#include "subtract.h"
#include "multiply.h"
#include "matrix.h"
#include <vector>
#include <cmath>

#define EPSILON_FOR_ZERO 0.00000001

namespace core
{
	template <class T, bool isParallel = false>
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

		bool is_equal(T a, T b)
		{
			// TODO - This should be implemented the right way
			return a==b;
		}

		apparatus(const tDiffInfo& di)
			: m_di(di)
		{}

		const tDiffInfo& getDiffInfo() const
		{
			return m_di;
		}

		void setDiffInfo(const tDiffInfo& di)
		{
			m_di = di;
		}
	
		//
		/// Apply the Taylor-diff transformation @a m_di on the @a sourceMatrix into the discretes @a discs
		//
		bool applyDiffTrans(const tFuncMatrix& sourceMatrix, tMatrixDiscretes& discs);

		//
		/// Add 2 discretes
		//
		bool addDiscretes(const tMatrixDiscretes& x, const tMatrixDiscretes& y, tMatrixDiscretes& out);

		//
		/// Subtract 2 discretes
		//
		bool subtractDiscretes(const tMatrixDiscretes& x, const tMatrixDiscretes& y, tMatrixDiscretes& out);

		//
		/// Multiply 2 discretes
		//
		bool multDiscretes(const tMatrixDiscretes& x, const tMatrixDiscretes& y, tMatrixDiscretes& out);

		//
		/// Restore the original @a out from the discretes @a theDiscretes with reverse single-point
		/// Taylor transformations
		//
		bool restoreTaylorSingle(const tMatrixDiscretes& theDiscretes, tFuncMatrix& out);

		//
		/// Get LU decomposition of the matrix
		//
		bool getLU(const tFuncMatrix& sourceMatrix, tMatrixDiscretes& l, tMatrixDiscretes& u);
	private:
		tDiffInfo m_di;
	};
} // namespace core

template <class T, bool isParallel>
bool core::apparatus<T, isParallel>::getLU(const tFuncMatrix& sourceMatrix, tMatrixDiscretes& l, tMatrixDiscretes& u)
{
	const int m = sourceMatrix.getNumRows();
	const int n = sourceMatrix.getNumCols();
	return true;
}

template <class T, bool isParallel>
bool core::apparatus<T, isParallel>::applyDiffTrans(const tFuncMatrix& sourceMatrix, tMatrixDiscretes& discs)
{
	const int m = sourceMatrix.getNumRows();
	const int n = sourceMatrix.getNumCols();
	discs.resize(m_di.K+1, tMatrixDiscrete(m, n, 0));
	tFuncMatrix theMatrix = sourceMatrix;
	long long kfac = 1;
	for (int k = 0; k <= m_di.K; ++k)
	{
		kfac *= (k > 0 ? k : 1);
		for (int i = 1; i <= m; ++i)
		{
			for (int j = 1; j <= n; ++j)
			{
				theMatrix[i][j] = theMatrix[i][j]->derivative(k > 0 ? 1 : 0);
				discs[k][i][j] = std::pow(m_di.H, k)*(*theMatrix[i][j])(m_di.tv)/(T)kfac;
			}
		}
	}
	return true;
}

template <class T, bool isParallel>
bool core::apparatus<T, isParallel>::addDiscretes(const tMatrixDiscretes& x, const tMatrixDiscretes& y, tMatrixDiscretes& out)
{
	assert (m_di.K == x.size()+1);
	assert (m_di.K == y.size()+1);
	const int xm = x.getNumRows();
	const int xn = x.getNumCols();
	const int ym = y.getNumRows();
	const int yn = y.getNumCols();
	if (xm != ym || xn != yn)
		throw matrixException("Matrix addition not possible, boundary error");
	out.resize(m_di.K+1, tMatrixDiscrete(xm, xn, 0));
	for (int k = 0; k <= m_di.K; ++k)
	{
		out[k] = x[k]+y[k];
	}
	return true;
}

template <class T, bool isParallel>
bool core::apparatus<T, isParallel>::subtractDiscretes(const tMatrixDiscretes& x, const tMatrixDiscretes& y, tMatrixDiscretes& out)
{
	assert (m_di.K == x.size()+1);
	assert (m_di.K == y.size()+1);
	const int xm = x.getNumRows();
	const int xn = x.getNumCols();
	const int ym = y.getNumRows();
	const int yn = y.getNumCols();
	if (xm != ym || xn != yn)
		throw matrixException("Matrix addition not possible, boundary error");
	out.resize(m_di.K+1, tMatrixDiscrete(xm, xn, 0));
	for (int k = 0; k <= m_di.K; ++k)
	{
		out[k] = x[k]-y[k];
	}
	return true;
}

template <class T, bool isParallel>
bool core::apparatus<T, isParallel>::multDiscretes(const tMatrixDiscretes& x, const tMatrixDiscretes& y, tMatrixDiscretes& out)
{
	assert (m_di.K == x.size()+1);
	assert (m_di.K == y.size()+1);
	const int xm = x.getNumRows();
	const int xn = x.getNumCols();
	const int ym = y.getNumRows();
	const int yn = y.getNumCols();
	if (xn != ym)
		throw matrixException("Matrix multiplication not possible, boundary error");
	out.resize(m_di.K+1, tMatrixDiscrete(xm, yn, 0));
	for (int k = 0; k <= m_di.K; ++k)
	{
		for (int l = 0; l <= k; ++l)
		{
			out[k] += x[k]*y[k-l];
		}
	}
	return true;
}

template <class T, bool isParallel>
bool core::apparatus<T, isParallel>::restoreTaylorSingle(const tMatrixDiscretes& theDiscretes, tFuncMatrix& out)
{
	assert (m_di.K+1 == theDiscretes.size());
	assert (m_di.K >= 0);
	const int m = out.getNumRows();
	const int n = out.getNumCols();
	for (int i = 1; i <= m; ++i)
	{
		for (int j = 1; j <= n; ++j)
		{
			out[i][j] = tFunctionPtr(new const_function<tArgType>(theDiscretes[0][i][j]));
			for (int k = 1; k <= m_di.K; ++k)
			{
				if (is_equal(theDiscretes[k][i][j], 0))
					continue;
				tFunctionPtr pow_base(
						new multiply<tArgType, eMultNum>(
								1./m_di.H,
								tFunctionPtr(
									new subtract<tArgType>(
										tFunctionPtr(new parameter<tArgType>()),
										tFunctionPtr(new const_function<tArgType>(m_di.tv))
										)
										)));
				out[i][j] = tFunctionPtr(
						new add<tArgType>(
								out[i][j],
								tFunctionPtr(
									new multiply<tArgType, eMultNum>(
										theDiscretes[k][i][j],
										(k <= 1)
											? pow_base
											: tFunctionPtr(new power<tArgType>(pow_base, k))))));
			}
		}
	}

	return true;
}




#endif // __APPARATUS_H__
