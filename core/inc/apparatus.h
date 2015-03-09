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

//#include <iostream>
//#include <algorithm>
//#include <iterator>

#define EPSILON_FOR_ZERO 0.00000001

namespace core
{
	enum eAlgoType
	{
		eAlgoBasic = 0,
		eAlgoParallel = 1
	};


	class algoException : public std::exception
	{
		std::string mErr;
	public:
		algoException(const std::string& err)
			: mErr(err)
		{}
		const char* what() const noexcept
		{
			return mErr.c_str();
		}
	};

	template <class T, int algo = eAlgoBasic>
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

		//
		/// Get LU decomposition of the matrix
		//
		bool getLU(const tFuncMatrix& sourceMatrix, int r, tMatrixDiscretes& l, tMatrixDiscretes& u);

		//
		/// Get LU decomposition of the matrix discretes
		//
		bool getLU(const tMatrixDiscretes& discretes, int r, tMatrixDiscretes& L, tMatrixDiscretes& U);

		//
		/// Get the rank of the matrix
		//
		int getRank(const tFuncMatrix& theMatrix) const;

		//
		/// Get inverse of the matrix
		//
		bool getInverse(const tFuncMatrix& theMatrix, tMatrixDiscretes& inv);

		//
		/// Get inverse of the matrix
		//
		bool getInverse(const tFuncMatrix& theMatrix, int r, tMatrixDiscretes& inv);

		//
		/// Get inverse of the matrix
		//
		bool getInverse(const tMatrixDiscretes& theMatrix, int r, tMatrixDiscretes& inv);

		//
		/// Get inverse of the matrix
		//
		bool getInverseLowerTriangular(const tFuncMatrix& theMatrix, tMatrixDiscretes& inv)
		{
			return this->getInverseLowerTriangular(theMatrix, getRank(theMatrix), inv);
		}

		//
		/// Get inverse of the matrix
		//
		bool getInverseLowerTriangular(const tFuncMatrix& theMatrix, int r, tMatrixDiscretes& inv);

		//
		/// Get inverse of the matrix
		//
		bool getInverseLowerTriangular(const tMatrixDiscretes& theMatrix, int r, tMatrixDiscretes& inv);

		//
		/// Get inverse of the matrix
		//
		bool getInverseUpperTriangular(const tFuncMatrix& theMatrix, tMatrixDiscretes& inv)
		{
			return this->getInverseUpperTriangular(theMatrix, getRank(theMatrix), inv);
		}

		//
		/// Get inverse of the matrix
		//
		bool getInverseUpperTriangular(const tFuncMatrix& theMatrix, int r, tMatrixDiscretes& inv);

		//
		/// Get inverse of the matrix
		//
		bool getInverseUpperTriangular(const tMatrixDiscretes& theMatrix, int r, tMatrixDiscretes& inv);

		//
		/// Check if the matrix is upper triangular or not
		//
		bool isUpperTriangular(const tFuncMatrix& theMatrix) const;

		//
		/// Check if the matrix is lower triangular or not
		//
		bool isLowerTriangular(const tFuncMatrix& theMatrix) const;

		//
		/// Check if the matrix is triangular or not
		//
		bool isTriangular(const tFuncMatrix& theMatrix) const
		{
			return this->isLowerTriangular(theMatrix) || this->isUpperTriangular(theMatrix);
		}

		//
		/// Check if the matrix is invertible or not
		//
		bool isInvertible(const tFuncMatrix& theMatrix) const
		{
			return this->isInvertible(theMatrix, this->getRank(theMatrix));
		}

		//
		/// Check if the matrix is invertible or not
		//
		bool isInvertible(const tFuncMatrix& theMatrix, int r) const;

		//
		/// Check if the matrix is invertible or not
		//
		bool isInvertible(const tMatrixDiscretes& theMatrix, int r) const;

		//
		/// Check if the matrix is square or not
		//
		bool isSquare(const tFuncMatrix& theMatrix) const
		{
			return theMatrix.getNumRows() == theMatrix.getNumCols();
		}

		//
		/// Check if the matrix is square or not
		//
		bool isSquare(const tMatrixDiscretes& discretes) const
		{
			assert (m_di.K >= 0);
			return discretes[0].getNumRows() == discretes[0].getNumCols();
		}

		//
		///
		//
	protected:
		enum eTriangularType
		{
			eLowerTriangular = 0,
			eUpperTriangular = 1,
			eNoneTriangular = 2
		};
		bool impl_getInverse(const tMatrixDiscretes& theMatrix, int r, tMatrixDiscretes& inv, eTriangularType et);
		T impl_identity(int k, int i, int j) const
		{
			return (k == 0 && i == j) ? 1 : 0;
		}
	private:
		tDiffInfo m_di;
	};
} // namespace core

template <class T, int algo>
int core::apparatus<T, algo>::getRank(const tFuncMatrix& theMatrix) const
{
	// TODO - We need to calculate the rank correctly!!!
	return theMatrix.getNumRows();
}

template <class T, int algo>
bool core::apparatus<T, algo>::getInverseLowerTriangular(const tFuncMatrix& theMatrix, int r, tMatrixDiscretes& inv)
{
	if (!this->isInvertible(theMatrix, r))
		throw algoException("Inversion: the matrix is not invertible");
	tMatrixDiscretes discretes;
	this->applyDiffTrans(theMatrix, discretes);
	return this->getInverseLowerTriangular(discretes, r, inv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::isInvertible(const tFuncMatrix& theMatrix, int r) const
{
	if (!isSquare(theMatrix))
		return false;
	if (r != theMatrix.getNumRows())
		return false;
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::isInvertible(const tMatrixDiscretes& theMatrix, int r) const
{
	if (!isSquare(theMatrix))
		return false;
	assert (m_di.K >= 0);
	if (r != theMatrix[0].getNumRows())
		return false;
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::getInverseLowerTriangular(const tMatrixDiscretes& L, int r, tMatrixDiscretes& L1)
{
	if (!this->isInvertible(L, r))
		throw algoException("Inversion: the matrix is not invertible");
	assert (m_di.K >= 0);
	const int n = L[0].getNumRows();
	L1.resize(m_di.K+1, tMatrixDiscrete(n, r, 0));
	
	int K = m_di.K;
	// Now L1
	for (int k = 0; k <= K; ++k)
	{
		for (int i = 1; i <= r; ++i)
		{
			for (int j = 1; j <= r; ++j)
			{
				double sum1 = 0;
				for (int p = 1; p <= i-1; ++p)
					for (int l = 0; l <= k; ++l)
						sum1 += L[l][i][p]*L1[k-l][p][j];

				double sum2 = 0;
				for (int l = 1; l <= k; ++l)
					sum2 += L1[k-l][i][j]*L[l][i][i];
				L1[k][i][j] = (impl_identity(k, i, j) - sum1 - sum2)/L[0][i][i];
			}
		}
	}
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::getInverseUpperTriangular(const tFuncMatrix& theMatrix, int r, tMatrixDiscretes& inv)
{
	if (!this->isInvertible(theMatrix, r))
		throw algoException("Inversion: the matrix is not invertible");
	tMatrixDiscretes discretes;
	this->applyDiffTrans(theMatrix, discretes);
	return this->getInverseUpperTriangular(discretes, r, inv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getInverseUpperTriangular(const tMatrixDiscretes& U, int r, tMatrixDiscretes& U1)
{
	if (!this->isInvertible(U, r))
		throw algoException("Inversion: the matrix is not invertible");
	assert (m_di.K >= 0);
	const int m = U[0].getNumCols();
	int K = m_di.K;
	U1.resize(K+1, tMatrixDiscrete(r, m, 0));
	// Now U1
	for (int k = 0; k <= K; ++k)
	{
		for (int i = r; i >= 1; --i)
		{
			for (int j = 1; j <= r; ++j)
			{
				double sum = 0;
				for (int p = i+1; p <= m; ++p)
					for (int g = 0; g <= k; ++g)
						sum += U[g][i][p]*U1[k-g][p][j];
				U1[k][i][j] = impl_identity(k, i, j) - sum;
			}
		}
	}
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::impl_getInverse(const tMatrixDiscretes& theMatrix, int r, tMatrixDiscretes& inv, eTriangularType et)
{
	switch (et)
	{
	case eLowerTriangular:
		return this->getInverseLowerTriangular(theMatrix, r, inv);
	break;
	case eUpperTriangular:
		return this->getInverseUpperTriangular(theMatrix, r, inv);
	break;
	default:
	break;
	}
	tMatrixDiscretes L;
	tMatrixDiscretes U;
	this->getLU(theMatrix, r, L, U);
	tMatrixDiscretes L1;
	this->getInverseLowerTriangular(L, r, L1);
	tMatrixDiscretes U1;
	this->getInverseUpperTriangular(U, r, U1);
	//std::cout << "L1:\n";
	//std::copy(L1.begin(), L1.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//std::cout << "U1:\n";
	//std::copy(U1.begin(), U1.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	return this->multDiscretes(U1, L1, inv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getInverse(const tMatrixDiscretes& theMatrix, int r, tMatrixDiscretes& inv)
{
	return this->impl_getInverse(theMatrix, r, inv, eNoneTriangular);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getInverse(const tFuncMatrix& theMatrix, tMatrixDiscretes& inv)
{
	return this->getInverse(theMatrix, getRank(theMatrix), inv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getInverse(const tFuncMatrix& theMatrix, int r, tMatrixDiscretes& inv)
{
	if (!this->isInvertible(theMatrix, r))
		throw algoException("Inversion: the matrix is not invertible");
	tMatrixDiscretes discretes;
	this->applyDiffTrans(theMatrix, discretes);
	return this->getInverse(discretes, r, inv);

}

template <class T, int algo>
bool core::apparatus<T, algo>::getLU(const tMatrixDiscretes& discretes, int r, tMatrixDiscretes& L, tMatrixDiscretes& U)
{
	int K = m_di.K;
	assert (K >= 0);
	const int m = discretes[0].getNumRows();
	const int n = discretes[0].getNumCols();
	if (m != n)
		throw matrixException("Matrix is not square!"); 
	
	if (is_equal(discretes[0][1][1], 0))
		throw algoException("LU: The element 1,1 is 0!");
	L.resize(m_di.K+1, tMatrixDiscrete(m, n, 0));
	U.resize(m_di.K+1, tMatrixDiscrete(m, n, 0));


	// First
	for (int k = 0; k < K; ++k)
		for (int i = 1; i <= r; ++i)
			L[k][i][1] = discretes[k][i][1];
	
	for (int k = 0; k < K; ++k)
	{
		for (int j = 1; j <= r ; ++j)
		{
			double sum = 0;
			for (int l = 1; l <= k; ++l)
			{
				sum += U[k-l][1][j]*L[l][1][1];
			}
			U[k][1][j] = (discretes[k][1][j]-sum)/L[0][1][1];
		}
	}

	// Next
	for (int k = 0; k < K; ++k)
	{
		for (int i = 2; i <= r; ++i)
		{
			for (int p = 2; p <= i; ++p)
			{
				double sum = 0;
				for (int j = 1; j <= p-1; ++j)
				{
					for (int l = 0; l <= k; ++l)
					{
						sum += L[l][i][j]*U[k-l][j][p];
					}
				}
				L[k][i][p] = discretes[k][i][p] - sum;
			}
		}
	}

	for (int k = 0; k < K; ++k)
	{
		for (int i = 2; i <= r; ++i)
		{
			for (int p = i; p <= r; ++p)
			{
				double sum1 = 0;
				for (int j = 1; j <= i-1; ++j)
					for (int g = 0; g <= k; ++g)
						sum1 += L[g][i][j]*U[k-g][j][p];

				double sum2 = 0;
				for (int l = 1; l <= k; ++l)
				{
					sum2 += U[k-l][i][p]*L[l][i][i];
				}
				U[k][i][p] = (discretes[k][i][p] - sum1 - sum2)/L[0][i][i];
			}
		}
	}
	
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::getLU(const tFuncMatrix& sourceMatrix, int r, tMatrixDiscretes& L, tMatrixDiscretes& U)
{
	if (!this->isSquare(sourceMatrix))
		throw matrixException("Matrix is not square!"); 
	tMatrixDiscretes discretes;
	this->applyDiffTrans(sourceMatrix, discretes);
	return this->getLU(discretes, r, L, U);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getLU(const tFuncMatrix& sourceMatrix, tMatrixDiscretes& L, tMatrixDiscretes& U)
{
	if (!this->isSquare(sourceMatrix))
		throw matrixException("Matrix is not square!"); 
	return this->getLU(sourceMatrix, this->getRank(sourceMatrix), L, U);
}

template <class T, int algo>
bool core::apparatus<T, algo>::applyDiffTrans(const tFuncMatrix& sourceMatrix, tMatrixDiscretes& discs)
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

template <class T, int algo>
bool core::apparatus<T, algo>::addDiscretes(const tMatrixDiscretes& x, const tMatrixDiscretes& y, tMatrixDiscretes& out)
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

template <class T, int algo>
bool core::apparatus<T, algo>::subtractDiscretes(const tMatrixDiscretes& x, const tMatrixDiscretes& y, tMatrixDiscretes& out)
{
	assert (m_di.K == x.size()+1);
	assert (m_di.K == y.size()+1);
	const int xm = x.getNumRows();
	const int xn = x.getNumCols();
	const int ym = y.getNumRows();
	const int yn = y.getNumCols();
	if (xm != ym || xn != yn)
		throw matrixException("Matrix subtraction not possible, boundary error");
	out.resize(m_di.K+1, tMatrixDiscrete(xm, xn, 0));
	for (int k = 0; k <= m_di.K; ++k)
	{
		out[k] = x[k]-y[k];
	}
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::multDiscretes(const tMatrixDiscretes& x, const tMatrixDiscretes& y, tMatrixDiscretes& out)
{
	assert (m_di.K+1 == x.size());
	assert (m_di.K+1 == y.size());
	assert (m_di.K >= 0);
	const int xm = x[0].getNumRows();
	const int xn = x[0].getNumCols();
	const int ym = y[0].getNumRows();
	const int yn = y[0].getNumCols();
	if (xn != ym)
		throw matrixException("Matrix multiplication not possible, boundary error");
	out.resize(m_di.K+1, tMatrixDiscrete(xm, yn, 0));
	for (int k = 0; k <= m_di.K; ++k)
	{
		for (int l = 0; l <= k; ++l)
		{
			out[k] += x[l]*y[k-l];
		}
	}
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::restoreTaylorSingle(const tMatrixDiscretes& theDiscretes, tFuncMatrix& out)
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
