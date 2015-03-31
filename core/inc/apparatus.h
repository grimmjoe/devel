#ifndef __APPARATUS_H__
#define __APPARATUS_H__

#include "function.h"
#include "constant_function.h"
#include "power_function.h"
#include "parameter_function.h"
#include "add.h"
#include "subtract.h"
#include "multiply.h"
#include "divide.h"
#include "matrix.h"
#include "comparator.h"
#include <vector>
#include <set>
#include <cmath>
#include <memory>
#include <iostream>

#include <chrono>

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
		typedef typename tFuncMatrix::tMatrixRow tFuncMatrixRow;
		typedef typename tFuncMatrix::tMatrixColumn tFuncMatrixColumn;
		typedef core::matrix<T> tMatrixDiscrete;
		typedef typename tMatrixDiscrete::tMatrixRow tMatrixDiscreteRow;
		typedef typename tMatrixDiscrete::tMatrixColumn tMatrixDiscreteColumn;
		typedef std::vector<tMatrixDiscrete> tMatrixDiscretes;
		typedef std::vector<T> tScalarDiscretes;

	
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

		bool is_equal(T a, T b) const
		{
			return m_comparator != nullptr 
						? m_comparator->is_equal(a, b)
						: a == b;
		}

		bool is_equal(const tMatrixDiscrete& a, T b) const
		{
			const int m = a.getNumRows();
			const int n = a.getNumCols();
			for (int i = 1; i <= m; ++i)
				for (int j = 1; j <= n; ++j)
					if (!this->is_equal(a[i][j], b))
						return false;
			return true;
		}

		bool is_equal(const tMatrixDiscrete& a, const tMatrixDiscrete& b) const
		{
			const int m = a.getNumRows();
			const int n = a.getNumCols();
			for (int i = 1; i <= m; ++i)
				for (int j = 1; j <= n; ++j)
					if (!this->is_equal(a[i][j], b[i][j]))
						return false;
			return true;
		}

		apparatus(const tDiffInfo& di)
			: m_di(di)
			, m_comparator(nullptr)
		{}

		const tDiffInfo& getDiffInfo() const
		{
			return m_di;
		}

		void setDiffInfo(const tDiffInfo& di)
		{
			m_di = di;
		}

		void setComparator(std::shared_ptr<comparator<T> > cm)
		{
			m_comparator = cm;
		}

		const std::shared_ptr<comparator<T> >& getComparator() const
		{
			return m_comparator;
		}
	
		//
		/// Apply the Taylor-diff transformation @a m_di on the @a A into the discretes @a discs
		//
		bool applyDiffTrans(const tFuncMatrix& A, tMatrixDiscretes& discs) const;

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
		/// Power up the discretes
		//
		bool powerDiscretes(const tMatrixDiscretes& x, int n, tMatrixDiscretes& out);

		//
		/// Multiply 2 scalar discretes
		//
		bool multDiscretes(const tScalarDiscretes& x, const tScalarDiscretes& y, tScalarDiscretes& out);

		//
		/// Multiply scalar and matrix discretes
		//
		bool multDiscretes(const tScalarDiscretes& x, const tMatrixDiscretes& y, tMatrixDiscretes& out);

		//
		/// Power up the scalar discretes
		//
		bool powerDiscretes(const tScalarDiscretes& x, int n, tScalarDiscretes& out);

		//
		/// Inverse discrete
		//
		bool inverseDiscrete(const tScalarDiscretes& x, tScalarDiscretes& out);

		//
		/// Restore the original @a out from the discretes @a theDiscretes with reverse single-point
		/// Taylor transformations
		//
		bool restoreTaylorSingle(const tMatrixDiscretes& theDiscretes, tFuncMatrix& out)
		{
			return this->restoreTaylorSingle(theDiscretes, out, m_di.K);
		}

		//
		/// Restore the original @a out from the discretes @a theDiscretes with reverse single-point
		/// Taylor transformations
		//
		bool restoreTaylorSingle(const tMatrixDiscretes& theDiscretes, tFuncMatrix& out, int K);

		//
		/// Restore the original @a out from the discretes @a theDiscretes with reverse Pade
		/// transformations
		//
		bool restorePade(const tMatrixDiscretes& theDiscretes, int m, int n, tFuncMatrix& out)
		{
			return this->restorePade(theDiscretes, m, n, out, m_di.K);
		}

		//
		/// Restore the original @a out from the discretes @a theDiscretes with reverse Pade
		/// transformations
		//
		bool restorePade(const tMatrixDiscretes& theDiscretes, int m, int n, tFuncMatrix& out, int K);

		//
		/// Get LU decomposition of the matrix
		//
		bool getLU(const tFuncMatrix& A, tMatrixDiscretes& l, tMatrixDiscretes& u);

		//
		/// Get LU decomposition of the matrix
		//
		bool getLU(const tFuncMatrix& A, int r, tMatrixDiscretes& l, tMatrixDiscretes& u);

		//
		/// Get LU decomposition of the matrix discretes
		//
		bool getLU(const tMatrixDiscretes& discretes, int r, tMatrixDiscretes& L, tMatrixDiscretes& U);

		//
		/// Get the rank of the matrix
		//
		int getRank(const tFuncMatrix& A) const;

		//
		/// Get the rank of the matrix
		//
		int getRank(const tMatrixDiscretes& A) const;

		bool getRankDecomposition(const tMatrixDiscrete& A, int& r, tMatrixDiscrete& C, tMatrixDiscrete& F);
	
	//// Inverses
	public:

		//
		/// Get inverse of the matrix
		//
		bool getInverse(const tFuncMatrix& A, tMatrixDiscretes& inv);

		//
		/// Get inverse of the matrix
		//
		bool getInverse(const tFuncMatrix& A, int r, tMatrixDiscretes& inv);

		//
		/// Get inverse of the matrix
		//
		bool getInverse(const tMatrixDiscretes& A, int r, tMatrixDiscretes& inv);

		//
		/// Get inverse of the matrix
		//
		bool getInverseLowerTriangular(const tFuncMatrix& A, tMatrixDiscretes& inv)
		{
			return this->getInverseLowerTriangular(A, getRank(A), inv);
		}

		//
		/// Get inverse of the matrix
		//
		bool getInverseLowerTriangular(const tFuncMatrix& A, int r, tMatrixDiscretes& inv);

		//
		/// Get inverse of the matrix
		//
		bool getInverseLowerTriangular(const tMatrixDiscretes& A, int r, tMatrixDiscretes& inv);

		//
		/// Get inverse of the matrix
		//
		bool getInverseUpperTriangular(const tFuncMatrix& A, tMatrixDiscretes& inv)
		{
			return this->getInverseUpperTriangular(A, getRank(A), inv);
		}

		//
		/// Get inverse of the matrix
		//
		bool getInverseUpperTriangular(const tFuncMatrix& A, int r, tMatrixDiscretes& inv);

		//
		/// Get inverse of the matrix
		//
		bool getInverseUpperTriangular(const tMatrixDiscretes& A, int r, tMatrixDiscretes& inv);

		//
		/// Check if the matrix is invertible or not
		//
		bool isInvertible(const tFuncMatrix& A) const
		{
			return this->isInvertible(A, this->getRank(A));
		}

		//
		/// Check if the matrix is invertible or not
		//
		bool isInvertible(const tFuncMatrix& A, int r) const;

		//
		/// Check if the matrix is invertible or not
		//
		bool isInvertible(const tMatrixDiscretes& A, int r) const;

	
	//// (B)-, (Q)- and (B, Q)-inverses
	public:
		//
		/// Get the (B)-inverse
		//
		bool getBInverse(const tFuncMatrix& A, tMatrixDiscretes& binv);

		//
		/// Get the (B)-inverse
		//
		bool getBInverse(const tFuncMatrix& A, int r, tMatrixDiscretes& binv);

		//
		/// Get the (B)-inverse
		//
		bool getBInverse(const tMatrixDiscretes& A, tMatrixDiscretes& binv);

		//
		/// Get the (B)-inverse
		//
		bool getBInverse(const tMatrixDiscretes& A, int r, tMatrixDiscretes& binv);

		//
		/// Get the (B)-inverse
		//
		bool getBInverse(const tFuncMatrix& A, const tFuncMatrix& B, tMatrixDiscretes& binv);

		//
		/// Get the (B)-inverse
		//
		bool getBInverse(const tFuncMatrix& A, int r, const tFuncMatrix& B, tMatrixDiscretes& binv);

		//
		/// Get the (B)-inverse
		//
		bool getBInverse(const tMatrixDiscretes& A, int r, const tMatrixDiscretes& B, tMatrixDiscretes& binv);

		//
		/// Check if the matrx is (B)-invertible
		//
		bool isBInvertible(const tFuncMatrix& A, int r) const;

		//
		/// Check if the matrx is (B)-invertible
		//
		bool isBInvertible(const tMatrixDiscretes& A, int r) const;

		//
		/// Check if the matrx is (B)-invertible
		//
		bool isBInvertible(const tFuncMatrix& A, const tFuncMatrix& B, int r) const;

		//
		/// Check if the matrx is (B)-invertible
		//
		bool isBInvertible(const tMatrixDiscretes& A, const tMatrixDiscretes& B, int r) const;

		//
		/// Choose the matrix B so that AB is invertible
		//
		bool chooseB(const tMatrixDiscretes& A, int r, tMatrixDiscretes& B);

		//
		/// Get the (Q)-inverse
		//
		bool getQInverse(const tFuncMatrix& A, tMatrixDiscretes& binv);

		//
		/// Get the (Q)-inverse
		//
		bool getQInverse(const tFuncMatrix& A, int r, tMatrixDiscretes& binv);

		//
		/// Get the (Q)-inverse
		//
		bool getQInverse(const tMatrixDiscretes& A, tMatrixDiscretes& binv);

		//
		/// Get the (Q)-inverse
		//
		bool getQInverse(const tMatrixDiscretes& A, int r, tMatrixDiscretes& binv);

		//
		/// Get the (Q)-inverse
		//
		bool getQInverse(const tFuncMatrix& A, const tFuncMatrix& Q, tMatrixDiscretes& binv);

		//
		/// Get the (Q)-inverse
		//
		bool getQInverse(const tFuncMatrix& A, int r, const tFuncMatrix& Q, tMatrixDiscretes& binv);

		//
		/// Get the (Q)-inverse
		//
		bool getQInverse(const tMatrixDiscretes& A, int r, const tMatrixDiscretes& Q, tMatrixDiscretes& binv);

		//
		/// Check if the matrx is (Q)-invertible
		//
		bool isQInvertible(const tFuncMatrix& A, const tFuncMatrix& Q, int r) const;

		//
		/// Check if the matrx is (Q)-invertible
		//
		bool isQInvertible(const tMatrixDiscretes& A, const tMatrixDiscretes& Q, int r) const;

		//
		/// Check if the matrx is (Q)-invertible
		//
		bool isQInvertible(const tFuncMatrix& A, int r) const;

		//
		/// Check if the matrx is (Q)-invertible
		//
		bool isQInvertible(const tMatrixDiscretes& A, int r) const;

		//
		/// Choose the matrix Q so that QA is invertible
		//
		bool chooseQ(const tMatrixDiscretes& A, int r, tMatrixDiscretes& Q);

		//
		/// Get the (BQ)-inverse
		//
		bool getBQInverse(const tFuncMatrix& A, tMatrixDiscretes& binv);

		//
		/// Get the (BQ)-inverse
		//
		bool getBQInverse(const tMatrixDiscretes& A, tMatrixDiscretes& binv);

		//
		/// Get the (BQ)-inverse
		//
		bool getBQInverse(const tFuncMatrix& A, const tFuncMatrix& S, tMatrixDiscretes& binv);

		//
		/// Get the (BQ)-inverse
		//
		bool getBQInverse(const tFuncMatrix& A, int r, const tFuncMatrix& S, tMatrixDiscretes& binv);

		//
		/// Get the (BQ)-inverse
		//
		bool getBQInverse(const tMatrixDiscretes& A, int r, const tMatrixDiscretes& S, tMatrixDiscretes& binv);

		//
		/// Check if the matrx is (BQ)-invertible
		//
		bool isBQInvertible(const tFuncMatrix& A) const;

		//
		/// Check if the matrx is (BQ)-invertible
		//
		bool isBQInvertible(const tFuncMatrix& A, int r) const;

		//
		/// Check if the matrx is (BQ)-invertible
		//
		bool isBQInvertible(const tMatrixDiscretes& A, int r) const;
		//
		/// Check if the A(t)Ainv(t)A(t)=A(t) condition applies or not for the first @a K discretes
		//
		bool checkB_Q_BQ_Inverse(const tMatrixDiscretes& A, const tMatrixDiscretes& inv)
		{
			return this->checkB_Q_BQ_Inverse(A, inv, m_di.K);
		}
		//
		/// Check if the A(t)Ainv(t)A(t)=A(t) condition applies or not for the first @a K discretes
		//
		bool checkB_Q_BQ_Inverse(const tMatrixDiscretes& A, const tMatrixDiscretes& inv, int K);
	
	/// Drazin inverse
	public:
		//
		/// Has the matrix a Drazin inverse or not
		//
		bool isDrazinInvertible(const tFuncMatrix& A) const;

		//
		/// Check if the Ad(t)A(t)Ad(t)=Ad(t0 and A(t)Ad(t)=Ad(t)A(t)
		//
		bool checkDrazinInverse(const tMatrixDiscretes& A, const tMatrixDiscretes& inv)
		{
			return this->checkDrazinInverse(A, inv, m_di.K);
		}
		//
		/// Check if the Ad(t)A(t)Ad(t)=Ad(t0 and A(t)Ad(t)=Ad(t)A(t)
		//
		bool checkDrazinInverse(const tMatrixDiscretes& A, const tMatrixDiscretes& inv, int K);

		//
		/// Has the matrix a Drazin inverse or not
		//
		bool isDrazinInvertible(const tMatrixDiscretes& A) const;

		//
		/// Get the Drazin inverse recursive
		//
		bool getDrazinInverseRecursive(const tFuncMatrix& A, tMatrixDiscretes& dinv);

		//
		/// Get the Drazin inverse recursive
		//
		bool getDrazinInverseRecursive(const tMatrixDiscretes& A, tMatrixDiscretes& dinv);

		//
		/// Get the Drazin inverse - skeleton decomposition
		//
		bool getDrazinInverseSkeleton(const tFuncMatrix& A, tMatrixDiscretes& dinv);

		//
		/// Get the Drazin inverse - skeleton decomposition
		//
		bool getDrazinInverseSkeleton(const tMatrixDiscretes& A, tMatrixDiscretes& dinv);

		//
		/// Get the Drazin inverse - canonical 
		//
		bool getDrazinInverseCanonical(const tFuncMatrix& A, tMatrixDiscretes& dinv);

		//
		/// Get the Drazin inverse - canonical
		//
		bool getDrazinInverseCanonical(const tMatrixDiscretes& A, tMatrixDiscretes& dinv);
	
	/// Utilities
	public:
		//
		/// Get the skeleton rep of the matrix
		//
		bool getSkeleton(const tMatrixDiscretes& A, tMatrixDiscretes& S, tMatrixDiscretes& R);

	
	//// Checks
	public:

		//
		/// Check if the matrix is upper triangular or not
		//
		bool isUpperTriangular(const tFuncMatrix& A) const;

		//
		/// Check if the matrix is lower triangular or not
		//
		bool isLowerTriangular(const tFuncMatrix& A) const;

		//
		/// Check if the matrix is triangular or not
		//
		bool isTriangular(const tFuncMatrix& A) const
		{
			return this->isLowerTriangular(A) || this->isUpperTriangular(A);
		}

		//
		/// Check if the matrix is square or not
		//
		bool isSquare(const tFuncMatrix& A) const
		{
			return A.getNumRows() == A.getNumCols();
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
		bool impl_getInverse(const tMatrixDiscretes& A, int r, tMatrixDiscretes& inv, eTriangularType et);
		T impl_identity(int k, int i, int j) const
		{
			return (k == 0 && i == j) ? 1 : 0;
		}

		int impl_getRank(const tMatrixDiscretes& A, tMatrixDiscretes& permMatrix) const;
		int impl_getRank(const tMatrixDiscrete& A, tMatrixDiscretes& permMatrix) const;

		T getScalarProduct(const tMatrixDiscreteColumn& c1, const tMatrixDiscreteColumn& c2) const
		{
			assert (c1.getNumRows() == c2.getNumRows());
			T sum = 0;
			for (int i = 1; i <= c1.getNumRows(); ++i)
			{
				sum += c1[i][1]*c2[i][1];
			}
			return sum;
		}
		T getNormNonSqrt(const tMatrixDiscreteColumn& c) const
		{
			T sum = 0;
			for (int i = 1; i <= c.getNumRows(); ++i)
			{
				sum += c[i][1]*c[i][1];
			}
			return sum;
		}
		T getNorm(const tMatrixDiscreteColumn& c) const
		{
			std::sqrt(this->getNormNonSqrt(c));
		}

		void makeIdentity(tMatrixDiscrete& d)
		{
			int m = std::min(d.getNumRows(), d.getNumCols());
			for (int i = 1; i <= m; ++i)
				d[i][i] = 1;
		}

		T getIdentity(int k, int i, int j)
		{
			return (k == 0) && (i == j) ? 1 : 0;
		}

		void checkVal(tMatrixDiscretes& discs, T val = 0)
		{
			std::for_each(discs.begin(), discs.end(),
				[&](tMatrixDiscrete& d)
				{
					this->checkVal(d, val);
				}
				);
		}

		void checkVal(tMatrixDiscrete& d, T val = 0)
		{
			const int m = d.getNumRows();
			const int n = d.getNumCols();
			for (int i = 1; i <= m; ++i)
				for (int j = 1; j <= n; ++j)
					if (this->is_equal(d[i][j], val))
						d[i][j] = val;
		}

	protected:
		tDiffInfo m_di;
		std::shared_ptr<comparator<T> > m_comparator;
	};
} // namespace core

template <class T, int algo>
int core::apparatus<T, algo>::getRank(const tFuncMatrix& A) const
{
	// TODO - We need to calculate the rank correctly!!!
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	return this->getRank(discretes);
}

template <class T, int algo>
int core::apparatus<T, algo>::getRank(const tMatrixDiscretes& A) const
{
	tMatrixDiscretes pm;
	return this->impl_getRank(A, pm);
}

template <class T, int algo>
int core::apparatus<T, algo>::impl_getRank(const tMatrixDiscretes& A, tMatrixDiscretes& permMatrix) const
{
	assert (m_di.K >= 0);
	return this->impl_getRank(A[0], permMatrix);
}

template <class T, int algo>
int core::apparatus<T, algo>::impl_getRank(const tMatrixDiscrete& A, tMatrixDiscretes& permMatrix) const
{
	const int m = A.getNumRows();
	const int n = A.getNumCols();
	typedef std::pair<tMatrixDiscreteColumn, T> tVecNormPair;
	typename tMatrixDiscrete::template numVector<tVecNormPair> d;
	d.reserve(n);
	d.push_back(std::make_pair(A.getColumn(1), 0));
	d[1].second = this->getNormNonSqrt(d[1].first);
	for (int j = 2; j <= n; ++j)
	{
		tMatrixDiscreteColumn aj = A.getColumn(j);
		d.push_back(std::make_pair(aj, 0));
		for (int i = 1; i < j; ++ i)
		{
			if (this->is_equal(d[i].second, 0))
				continue;
			T scalar = this->getScalarProduct(aj, d[i].first);
			if (scalar == 0)
				continue;
			//d[j].first -= (scalar*d[i].first)/d[i].second;
			d[j].first -= (scalar/d[i].second)*d[i].first;
		}
		d[j].second = this->getNormNonSqrt(d[j].first);
		//if (this->is_equal(d[j].second, 0))
		//{
		//	d[j].second = 0;
		//	d[j].first.set(0);
		//}
	}
	//std::cout << "Orthogonalization was:\n";
	//for (int i = 1; i <= d.size(); ++i)
	//{
	//	std::cout << i << ":\n" << d[i].first << std::endl;
	//	std::cout << d[i].second << std::endl;
	//}
	std::vector<int> theIndices(n, 0);
	int r = 0;
	int nn = n;
	for (int i = 1; i <= d.size(); ++i)
	{
		//if (!is_equal(d[i].second, 0))
		//	theIndices[r++] = i;
		//else
		//	theIndices[--nn] = i;
		int j = 1;
		for (; j <= m; ++j)
		{
			if (!is_equal(d[i].first[j], 0))
			{
				theIndices[r++] = i;
				break;
			}
		}
		if (j > m)
			theIndices[--nn] = i;
	}
	//std::cout << "The rank is " << r << std::endl;
	//assert (r <= std::min(m, n));
	//permMatrix.resize(n, n, 0);
	permMatrix.resize(1, tMatrixDiscrete(n, r, 0));
	for (int i = 1; i <= r; ++i)
	{
		permMatrix[0][theIndices[i-1]][i]=1;
	}
	//int i = 1;
	//std::for_each(theIndices.begin(), theIndices.end(), 
	//	[&](int ind)
	//	{
	//		std::cout << "ind = " << ind << std::endl;
	//		permMatrix[ind][i++]=1;
	//	}
	//	);
	
	//std::cout << "The indices:\n";
	//std::copy(theIndices.begin(), theIndices.end(), std::ostream_iterator<int>(std::cout, " "));
	//std::cout << "The permutation matrix is:\n";
	//std::cout << permMatrix[0] << std::endl;

	//return r;
	return std::min(r, std::min(m, n));
}

template <class T, int algo>
bool core::apparatus<T, algo>::isDrazinInvertible(const tFuncMatrix& A) const
{
	return this->isSquare(A);
}

template <class T, int algo>
bool core::apparatus<T, algo>::isDrazinInvertible(const tMatrixDiscretes& A) const
{
	return this->isSquare(A);
}


template <class T, int algo>
bool core::apparatus<T, algo>::getDrazinInverseRecursive(const tFuncMatrix& A, tMatrixDiscretes& dinv)
{
	if (!this->isDrazinInvertible(A))
		throw algoException("Drazin Inversion: the matrix is not Drazin-invertible");
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	return this->getDrazinInverseRecursive(discretes, dinv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getDrazinInverseRecursive(const tMatrixDiscretes& A, tMatrixDiscretes& dinv)
{
	if (!this->isDrazinInvertible(A))
		throw algoException("Drazin Inversion: the matrix is not Drazin-invertible");
	//std::cout << "Here\n";
	assert (m_di.K >= 0);
	const int n = A[0].getNumRows();
	std::vector<tMatrixDiscretes> S(n+1, tMatrixDiscretes());
	S[0].resize(m_di.K+1, tMatrixDiscrete(n, n, 0));
	this->makeIdentity(S[0][0]);
	int j = 1;
	int p = 0;
	bool toStop = j > n;
	std::vector<tScalarDiscretes> betta(n, tScalarDiscretes(m_di.K+1, 0));
	tMatrixDiscrete I(n, n, 0);
	this->makeIdentity(I);
	//std::cout << "Entering the loop\n";
	while (!toStop)
	{
		//std::cout << "j = " << j << std::endl;
		//std::cout << "S[" << j-1 << "]:\n";
		//std::copy(S[j-1].begin(), S[j-1].end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
		p = 0;
		tMatrixDiscretes AS;
		this->multDiscretes(A, S[j-1], AS);
		//std::cout << "AS:\n";
		//std::copy(AS.begin(), AS.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, " "));
		for (int k = 0; k <= m_di.K; ++k)
		{
			betta[n-j][k] = (-1./j) * trace(AS[k]);
			//std::cout << "Trace(AS[" << k << "]=" << trace(AS[k]) << std::endl;
			if (this->is_equal(betta[n-j][k], 0))
				betta[n-j][k] = 0;
		}
		//std::cout << "Calculated bettas, now calculating S\n";
		//std::cout << "betta[" << n-j << "]:\n";
		//std::copy(betta[n-j].begin(), betta[n-j].end(), std::ostream_iterator<T>(std::cout, " "));
		//std::cout << std::endl;
		this->multDiscretes(A, S[j-1], S[j]);
		for (int k = 0; k <= m_di.K; ++k)
		{
			if (!this->is_equal(betta[n-j][k], 0))
				S[j][k] += betta[n-j][k]*I;
			if (!this->is_equal(S[j][k], 0))
				p = -1;
		}
		//std::cout << "Calculated S\n";
		//std::cout << "p = " << p << std::endl;
		if (p != -1)
		{
			p = j;
			break;
		}
		++j;
		toStop = j > n;
	}
	//std::cout << "p = " << p << std::endl;
	//std::cout << "betta:\n";
	//for (int i = 0; i < n; ++i)
	//{
	//	std::cout << "betta[" << i << "]:\n";
	//	std::copy(betta[i].begin(), betta[i].end(), std::ostream_iterator<T>(std::cout, " "));
	//	std::cout << std::endl;
	//	std::cout << "S[" << i << "]:\n";
	//	std::copy(S[i].begin(), S[i].end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, " "));
	//	std::cout << std::endl;
	//}
	//	std::cout << "S[" << n << "]:\n";
	//	std::copy(S[n].begin(), S[n].end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, " "));
	//	std::cout << std::endl;
	int u = n;
	for (; u >= 1; --u)
	{
		int k = 0;
		for (; k <= m_di.K; ++k)
		{
			if (!this->is_equal(betta[n-u][k], 0))
				break;
		}
		if (k <= m_di.K)
			break;
	}
	assert (u >= 1);
	//std::cout << "u = " << u << std::endl;
	int l = n-u;
	tMatrixDiscretes al;
	this->powerDiscretes(A, l, al);
	tMatrixDiscretes su;
	this->powerDiscretes(S[u-1], l+1, su);
	//std::cout << "AL:\n";
	//std::copy(al.begin(), al.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, " "));
	//std::cout << std::endl;

	//std::cout << "S[u-1]^l+1:\n";
	//std::copy(su.begin(), su.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, " "));
	//std::cout << std::endl;

	tMatrixDiscretes alsu;
	this->multDiscretes(al, su, alsu);
	tScalarDiscretes bettal;
	this->powerDiscretes(betta[l], l+1, bettal);
	//std::cout << "bettal:\n";
	//std::copy(bettal.begin(), bettal.end(), std::ostream_iterator<T>(std::cout, " "));
	//std::cout << std::endl;
	// TODO - This might not be true but it seams it is
	if ( n % 2 == 0)
	{
		for (int k = 0; k <= m_di.K; ++k)
			bettal[k] *= -1;
	}
	tScalarDiscretes bettalin;
	this->inverseDiscrete(bettal, bettalin);
	//std::cout << "bettal-inverse:\n";
	//std::copy(bettalin.begin(), bettalin.end(), std::ostream_iterator<T>(std::cout, " "));
	this->multDiscretes(bettalin, alsu, dinv);

	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::getDrazinInverseSkeleton(const tFuncMatrix& A, tMatrixDiscretes& dinv)
{
	if (!this->isDrazinInvertible(A))
		throw algoException("Drazin Inversion: the matrix is not Drazin-invertible");
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	return this->getDrazinInverseRecursive(discretes, dinv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getSkeleton(const tMatrixDiscretes& A, tMatrixDiscretes& S, tMatrixDiscretes& R)
{
	tMatrixDiscretes permMatrix;
	int r = this->impl_getRank(A, permMatrix);
	S.reserve(A.size());
	std::for_each(A.begin(), A.end(), 
		[&](const tMatrixDiscrete& d)
		{
			S.push_back(d*permMatrix[0]);
		}
		);
	//std::cout << "S:\n";
	//std::copy(S.begin(), S.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes Q;
	this->chooseQ(S, r, Q);
	//std::cout << "Chose Q as:\n";
	//std::copy(Q.begin(), Q.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes qinv;
	this->getQInverse(S, r, Q, qinv);
	//std::cout << "Q-Inverse is:\n";
	//std::copy(qinv.begin(), qinv.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	this->multDiscretes(qinv, A, R);
	//std::cout << "R is:\n";
	//std::copy(R.begin(), R.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::getDrazinInverseSkeleton(const tMatrixDiscretes& A, tMatrixDiscretes& dinv)
{
	assert (m_di.K >= 0);
	const int n = A[0].getNumRows();
	tMatrixDiscretes an1;
	this->powerDiscretes(A, n-1, an1);
	tMatrixDiscretes an;
	this->multDiscretes(an1, A, an);
	//std::cout << "An:\n";
	//std::copy(an.begin(), an.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes B;
	tMatrixDiscretes C;
	this->getSkeleton(an, B, C);
	int r = B[0].getNumCols();
	tMatrixDiscretes CB;
	this->multDiscretes(C, B, CB);
	//std::cout << "CB:\n";
	//std::copy(CB.begin(), CB.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//std::cout << std::endl;
	tMatrixDiscretes X1;
	this->getInverse(CB, r, X1);
	//std::cout << "X1:\n";
	//std::copy(X1.begin(), X1.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes x12;
	this->multDiscretes(X1, X1, x12);
	//std::cout << "Got x12\n";
	tMatrixDiscretes an1B;
	this->multDiscretes(an1, B, an1B);
	//std::cout << "Got an1B\n";
	tMatrixDiscretes an1Bx12;
	this->multDiscretes(an1B, x12, an1Bx12);
	//std::cout << "Got an1Bx12\n";
	return this->multDiscretes(an1Bx12, C, dinv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getDrazinInverseCanonical(const tFuncMatrix& A, tMatrixDiscretes& dinv)
{
	if (!this->isDrazinInvertible(A))
		throw algoException("Drazin Inversion: the matrix is not Drazin-invertible");
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	return this->getDrazinInverseCanonical(discretes, dinv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getDrazinInverseCanonical(const tMatrixDiscretes& A, tMatrixDiscretes& dinv)
{
	assert (m_di.K >= 0);
	const int n = A[0].getNumRows();
	tMatrixDiscretes an;
	this->powerDiscretes(A, n, an);
	//std::cout <<"An:\n";
	//std::copy(an.begin(), an.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes R = an;
	for (int i0 = 1; i0 <= n; ++i0)
	{
		int j0 = i0;
		//std::cout << "i0 = " << i0 << std::endl;
		////int imax = i0;
		//for (int i = i0+1; i <= n; ++i)
		//{
		//	if (R[0][i][j0] > R[0][imax][j0])
		//		imax = i;
		//}
		//R[0].swap(imax, i0);
		if (this->is_equal(R[0][i0][j0], 0))
		{
			//std::cout << "i0=" << i0 << ",j0=" << j0 << " is NULL\n";
			continue;
		}
		tScalarDiscretes rk(m_di.K+1, 0);
		for (int k = 0; k <= m_di.K; ++k)
			rk[k] = R[k][i0][j0];
		for (int k = 0; k <= m_di.K; ++k)
		{
			for (int j = j0; j <= n; ++j)
			{
				T sum = 0;
				for (int l = 1; l <= k; ++l)
				{
					sum += R[k-l][i0][j]*rk[l];
				}
				R[k][i0][j] = (R[k][i0][j]-sum)/rk[0];
				if (this->is_equal(R[k][i0][j], 0))
					R[k][i0][j]=0;
			}
		}

		//std::cout << "R[0] = " << R[0] << std::endl;
		//std::cout << "rk[0] = " << rk[0] << std::endl;

		for (int i = 1; i <= n; ++i)
		{
			if (i == i0)
				continue;
			//std::cout << "i = " << i << std::endl;
			std::vector<tScalarDiscretes> b(m_di.K+1, tScalarDiscretes(n-j0+1, 0));
			for (int k = 0; k <= m_di.K; ++k)
			{
				int bi = 0;
				for (int j = j0; j <= n; ++j)
				{
					T sum = 0;
					for (int l = 0; l <= k; ++l)
					{
						sum += R[l][i0][j]*R[k-l][i][j0];
					}
					b[k][bi++] = sum;
				}
			}
			for (int k = 0; k <= m_di.K; ++k)
			{
				int bi = 0;
				for (int j = j0; j <= n; ++j)
				{
					R[k][i][j] = R[k][i][j]-b[k][bi++];
					if (this->is_equal(R[k][i][j], 0))
						R[k][i][j] = 0;
				}
			}
		}

		//std::cout << "R[0] = " << R[0] << std::endl;
	}
	//std::cout << "R:\n";
	//std::copy(R.begin(), R.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//std::cout << std::endl;
	tMatrixDiscretes P(m_di.K+1, tMatrixDiscrete(n, n, 0));
	int pj = 1;
	std::set<int> theCols;
	for (int j = 1; j <= n; ++j)
	{
		int k = 0;
		for (; k <= m_di.K; ++k)
			if (!this->is_equal(R[k][j][j], 0))
				break;

		if (k <= m_di.K)
		{
			for (int k = 0; k <= m_di.K; ++k)
				for (int i = 1; i <= n; ++i)
					P[k][i][pj] = an[k][i][j];
			++pj;
			theCols.insert(j);
		}
	}
	for (int j = 1; j <= n; ++j)
	{
		if (theCols.find(j) != theCols.end())
			continue;
		for (int k = 0; k <= m_di.K; ++k)
			for (int i = 1; i <= n; ++i)
				P[k][i][pj]=this->getIdentity(k, i, j)-R[k][i][j];
		++pj;
	}
	//std::cout << "P:\n";
	//std::copy(P.begin(), P.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes P1;
	this->getInverse(P, n, P1);
	//std::cout << "P1:\n";
	//std::copy(P1.begin(), P1.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes P1A;
	this->multDiscretes(P1, A, P1A);
	tMatrixDiscretes TM;
	this->multDiscretes(P1A, P, TM);
	//std::cout << "TM:\n";
	//std::copy(TM.begin(), TM.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	int cn = n;
	T sum = 0;
	for (; cn >= 1; --cn)
	{
		if (this->is_equal(TM[0][1][cn], 0))
			TM[0][1][cn] = 0;
		if (this->is_equal(TM[0][cn][1], 0))
			TM[0][cn][1] = 0;
		if (!this->is_equal(TM[0][1][cn], 0))
			break;
		else
		{
			assert (this->is_equal(TM[0][1][cn], TM[0][cn][1]));
		}
		if (this->is_equal(TM[0][cn][cn], 0))
			TM[0][cn][cn] = 0;
		sum += TM[0][cn][cn];
		if (sum == 0)
		{
			--cn;
			break;
		}
	}
	//std::cout << "cn = " << cn << std::endl;
	//std::cout << "sum = " << sum << std::endl;
	assert (cn >= 1);
	if (!this->is_equal(sum, 0))
	{
		++cn;
		while (!this->is_equal(sum, 0) && cn <= n)
		{
			sum -= TM[0][cn][cn];
			++cn;
		}
		--cn;
	}
	//std::cout << "cn = " << cn << std::endl;
	assert (cn >= 1);
	tMatrixDiscretes C(m_di.K+1, tMatrixDiscrete(cn, cn, 0));
	for (int k = 0; k <= m_di.K; ++k)
	{
		for (int i = 1; i <= cn; ++i)
			for (int j = 1; j <= cn; ++j)
				C[k][i][j] = TM[k][i][j];
	}
	//std::cout << "C:\n";
	//std::copy(C.begin(), C.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes C1;
	this->getInverse(C, cn, C1);
	//std::cout << "C1:\n";
	//std::copy(C1.begin(), C1.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes C1En(m_di.K+1, tMatrixDiscrete(n, n, 0));
	for (int k = 0; k <= m_di.K; ++k)
		for (int i = 1; i <= cn; ++i)
			for (int j = 1; j <= cn; ++j)
				C1En[k][i][j] = C1[k][i][j];
	tMatrixDiscretes PC1;
	this->multDiscretes(P, C1En, PC1);
	this->multDiscretes(PC1, P1, dinv);

	//if (pj != n)
	//{
	//	//tMatrixDiscretes I(m_di.K+1, tMatrixDiscrete(n, n, 0));
	//	//this->makeIdentity(I[0]);
	//	//tMatrixDiscretes IH;
	//	//this->subtractDiscretes(I, R, IH);
	//	for (int k = 0; k <= m_di.K; ++k)
	//	{
	//		for (int i
	//	}
	//}

	//for (int i0 = 1; i0 <= n; ++i0)
	//{
	//	int j0 = i0;
	//	int imax = i0;
	//	for (int i = i0+1; i <= n; ++i)
	//	{
	//		if (R[0][i][j0] > R[0][imax][j0])
	//			imax = i;
	//	}
	//	for (int k = 0; k <= m_di.K; ++k)
	//	{
	//		R[k].swap(i0, imax);
	//		if (this->is_equal(R[k][i0][j0], 0))
	//			
	//		T sum = 0;
	//		for (int l = 1; l <= k; ++l)
	//		{
	//			sum += R[K-l][i0][j]*R[l][i0][j0];
	//		}
	//		R[k][i0][j] = (A[k][i0][j]-sum)/R[k][i0][j0];
	//	}
	//}

	//tMatrixDiscretes an1;
	//this->powerDiscretes(A, n-1, an1);
	//tMatrixDiscretes an;
	//this->multDiscretes(an1, A, an);
	//std::cout << "An:\n";
	//std::copy(an.begin(), an.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//tMatrixDiscretes B;
	//tMatrixDiscretes C;
	//this->getSkeleton(an, B, C);
	//int r = B[0].getNumCols();
	//tMatrixDiscretes CB;
	//this->multDiscretes(C, B, CB);
	//std::cout << "CB:\n";
	//std::copy(CB.begin(), CB.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//std::cout << std::endl;
	//tMatrixDiscretes X1;
	//this->getInverse(CB, r, X1);
	//std::cout << "X1:\n";
	//std::copy(X1.begin(), X1.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//tMatrixDiscretes x12;
	//this->multDiscretes(X1, X1, x12);
	//std::cout << "Got x12\n";
	//tMatrixDiscretes an1B;
	//this->multDiscretes(an1, B, an1B);
	//std::cout << "Got an1B\n";
	//tMatrixDiscretes an1Bx12;
	//this->multDiscretes(an1B, x12, an1Bx12);
	//std::cout << "Got an1Bx12\n";
	//return this->multDiscretes(an1Bx12, C, dinv);
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::isBQInvertible(const tFuncMatrix& A) const
{
	int r = this->getRank(A);
	return (r < std::min(A.getNumRows(), A.getNumCols()))
			|| this->isInvertible(A, r);
}

template <class T, int algo>
bool core::apparatus<T, algo>::isBQInvertible(const tFuncMatrix& A, int r) const
{
	return (r < std::min(A.getNumRows(), A.getNumCols()))
			|| this->isInvertible(A, r);
}


template <class T, int algo>
bool core::apparatus<T, algo>::isBQInvertible(const tMatrixDiscretes& A, int r) const
{
	assert (m_di.K >= 0);
	return (r < std::min(A[0].getNumRows(), A[0].getNumCols()))
			|| this->isInvertible(A, r);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getBQInverse(const tFuncMatrix& A, tMatrixDiscretes& binv)
{
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	return this->getBQInverse(discretes, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getBQInverse(const tMatrixDiscretes& A, tMatrixDiscretes& binv)
{
	tMatrixDiscretes permMatrix;
	int r = this->impl_getRank(A, permMatrix);
	if (!this->isBQInvertible(A, r))
		throw algoException("(B, Q)-Inversion: The matrix is not (B, Q)-invertible");
	tMatrixDiscretes sdiscs;
	sdiscs.reserve(A.size());
	std::for_each(A.begin(), A.end(), 
		[&](const tMatrixDiscrete& d)
		{
			sdiscs.push_back(d*permMatrix[0]);
		}
		);
	return this->getBQInverse(A, r, sdiscs, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getBQInverse(const tFuncMatrix& A, int r, const tFuncMatrix& S, tMatrixDiscretes& binv)
{
	if (!this->isBQInvertible(A, r))
		throw algoException("(Q)-Inversion: The matrix is not (B, Q)-invertible");
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	tMatrixDiscretes sdiscs;
	this->applyDiffTrans(S, sdiscs);
	return this->getBQInverse(discretes, r, sdiscs, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getBQInverse(const tMatrixDiscretes& A, int r, const tMatrixDiscretes& S, tMatrixDiscretes& binv)
{
	if (!this->isBQInvertible(A, r))
		throw algoException("(Q)-Inversion: The matrix is not (B, Q)-invertible");

	//std::cout << "A:\n" << std::endl;
	//std::copy(A.begin(), A.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes Q;
	this->chooseQ(S, r, Q);
	//std::cout << "Chose Q as:\n";
	//std::copy(Q.begin(), Q.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes qinv;
	this->getQInverse(S, r, Q, qinv);
	//std::cout << "Q-Inverse is:\n";
	//std::copy(qinv.begin(), qinv.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes R;
	this->multDiscretes(qinv, A, R);
	//std::cout << "R is:\n";
	//std::copy(R.begin(), R.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes B;
	this->chooseB(R, r, B);
	//std::cout << "Chose B as:\n";
	//std::copy(B.begin(), B.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes rinv;
	this->getBInverse(R, r, B, rinv);
	//std::cout << "B-Inverse is:\n";
	//std::copy(rinv.begin(), rinv.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	return this->multDiscretes(rinv, qinv, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::checkB_Q_BQ_Inverse(const tMatrixDiscretes& A, const tMatrixDiscretes& inv, int K)
{
	tMatrixDiscretes a_ainv;
	this->multDiscretes(A, inv, a_ainv);
	tMatrixDiscretes a_mult;
	this->multDiscretes(a_ainv, A, a_mult);
	for (int k = 0; k <= K; ++k)
	{
		//std::cout << "k = " << k << std::endl;
		//std::cout << A[k] << std::endl;
		//std::cout << a_mult[k] << std::endl;
		if (!this->is_equal(A[k],a_mult[k]))
		{
			std::cerr << "Discretes not equal at " << k << std::endl;
			return false;
		}
	}
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::checkDrazinInverse(const tMatrixDiscretes& A, const tMatrixDiscretes& inv, int K)
{
	// First
	tMatrixDiscretes ainv_a;
	this->multDiscretes(inv, A, ainv_a);
	tMatrixDiscretes a_mult;
	this->multDiscretes(ainv_a, inv, a_mult);
	for (int k = 0; k <= K; ++k)
	{
		//std::cout << "k = " << k << std::endl;
		//std::cout << inv[k] << std::endl;
		//std::cout << a_mult[k] << std::endl;
		if (!this->is_equal(inv[k], a_mult[k]))
		{
			std::cerr << "Discretes not equal at " << k << std::endl;
			std::cout << inv[k] << std::endl;
			std::cout << a_mult[k] << std::endl;
			return false;
		}
	}
	// Second
	tMatrixDiscretes a_ainv;
	this->multDiscretes(A, inv, a_ainv);
	for (int k = 0; k <= K; ++k)
	{
		//std::cout << "k = " << k << std::endl;
		//std::cout << a_ainv[k] << std::endl;
		//std::cout << ainv_a[k] << std::endl;
		if (!this->is_equal(a_ainv[k], ainv_a[k]))
		{
			std::cerr << "Discretes not equal at " << k << std::endl;
			std::cout << a_ainv[k] << std::endl;
			std::cout << ainv_a[k] << std::endl;
			return false;
		}
	}
	return true;
}


template <class T, int algo>
bool core::apparatus<T, algo>::isBInvertible(const tFuncMatrix& A, const tFuncMatrix& B, int r) const
{
	// TODO - What about the matrix B?
	// TODO - the product A*B must be invertible
	return (A.getNumRows() <= A.getNumCols()) && (r == A.getNumRows())
			&& (B.getNumRows() == A.getNumCols()) && (B.getNumCols() == A.getNumRows());
}

template <class T, int algo>
bool core::apparatus<T, algo>::isBInvertible(const tMatrixDiscretes& A, const tMatrixDiscretes& B, int r) const
{
	// TODO - What about the matrix B?
	// TODO - the product A*B must be invertible
	assert (m_di.K >= 0);
	return (A[0].getNumRows() <= A[0].getNumCols()) && (r == A[0].getNumRows())
			&& (B[0].getNumRows() == A[0].getNumCols()) && (B[0].getNumCols() == A[0].getNumRows());
}

template <class T, int algo>
bool core::apparatus<T, algo>::isBInvertible(const tFuncMatrix& A, int r) const
{
	// TODO - What about the matrix B?
	// TODO - the product A*B must be invertible
	return (A.getNumRows() <= A.getNumCols()) && (r == A.getNumRows());
}

template <class T, int algo>
bool core::apparatus<T, algo>::isBInvertible(const tMatrixDiscretes& A, int r) const
{
	// TODO - What about the matrix B?
	// TODO - the product A*B must be invertible
	assert (m_di.K >= 0);
	return (A[0].getNumRows() <= A[0].getNumCols()) && (r == A[0].getNumRows());
}

template <class T, int algo>
bool core::apparatus<T, algo>::getBInverse(const tFuncMatrix& A, tMatrixDiscretes& binv)
{
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	int r = this->getRank(discretes);
	return this->getBInverse(discretes, r, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getBInverse(const tFuncMatrix& A, int r, tMatrixDiscretes& binv)
{
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	return this->getBInverse(discretes, r, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getBInverse(const tMatrixDiscretes& A, tMatrixDiscretes& binv)
{
	return this->getBInverse(A, this->getRank(A), binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getBInverse(const tMatrixDiscretes& A, int r, tMatrixDiscretes& binv)
{
	tMatrixDiscretes B;
	this->chooseB(A, r, B);
	return this->getBInverse(A, r, B, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getBInverse(const tFuncMatrix& A, const tFuncMatrix& B, tMatrixDiscretes& binv)
{
	return this->getBInverse(A, this->getRank(A), B, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getBInverse(const tFuncMatrix& A, int r, const tFuncMatrix& B, tMatrixDiscretes& binv)
{
	if (!this->isBInvertible(A, B, r))
		throw algoException("(B)-Inversion: The matrix is not (B)-invertible");
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	tMatrixDiscretes bDiscs;
	this->applyDiffTrans(B, bDiscs);
	return this->getBInverse(discretes, r, bDiscs, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getBInverse(const tMatrixDiscretes& A, int r, const tMatrixDiscretes& B, tMatrixDiscretes& binv)
{
	if (!this->isBInvertible(A, B, r))
		throw algoException("(B)-Inversion: The matrix is not (B)-invertible");
	
	//std::cout << "A:\n";
	//std::copy(A.begin(), A.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//std::cout << "B:\n";
	//std::copy(B.begin(), B.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//
	// Calculate the product of A*B
	//
	tMatrixDiscretes AB;
	this->multDiscretes(A, B, AB);
	//std::cout << "AB:\n";
	//std::copy(AB.begin(), AB.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//
	// Calculate the inverse to it
	//
	tMatrixDiscretes AB1;
	this->getInverse(AB, r, AB1);
	//std::cout << "AB1:\n";
	//std::copy(AB1.begin(), AB1.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	////
	// Calculate the (B)-inverse
	//
	return this->multDiscretes(B, AB1, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::chooseB(const tMatrixDiscretes& A, int r, tMatrixDiscretes& B)
{
	assert (m_di.K >= 0);
	(void) r;
	B.resize(m_di.K+1, tMatrixDiscrete(A[0].getNumCols(), A[0].getNumRows(), 0));
	// Transpose the martix. Apply Gramm-Schmidt. Get the permutation matrix
	// Transpose the permutation matrix and multiply with the original to get
	// a row permutation such that the resulted matrix is invertible
	//tMatrixDiscrete Atrans = A[0].getTranspose();
	tMatrixDiscretes permMatrix;
	(void) this->impl_getRank(A, permMatrix);
	B[0] = permMatrix[0];
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::chooseQ(const tMatrixDiscretes& A, int r, tMatrixDiscretes& Q)
{
	assert (m_di.K >= 0);
	(void) r;
	Q.resize(m_di.K+1, tMatrixDiscrete(A[0].getNumCols(), A[0].getNumRows(), 0));
	// Transpose the martix. Apply Gramm-Schmidt. Get the permutation matrix
	// Transpose the permutation matrix and multiply with the original to get
	// a row permutation such that the resulted matrix is invertible
	tMatrixDiscrete Atrans = A[0].getTranspose();
	tMatrixDiscretes permMatrix;
	(void) this->impl_getRank(Atrans, permMatrix);
	Q[0] = permMatrix[0].getTranspose();
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::isQInvertible(const tFuncMatrix& A, const tFuncMatrix& Q, int r) const
{
	// TODO - What about the matrix Q?
	// TODO - the product Q*A must be invertible
	return (r == A.getNumCols()) && (A.getNumCols() <= A.getNumRows())
			&& (Q.getNumRows() == A.getNumCols()) && (Q.getNumCols() == A.getNumRows());
}

template <class T, int algo>
bool core::apparatus<T, algo>::isQInvertible(const tMatrixDiscretes& A, const tMatrixDiscretes& Q, int r) const
{
	// TODO - What about the matrix Q?
	// TODO - the product Q*A must be invertible
	assert (m_di.K >= 0);
	return (r == A[0].getNumCols()) && (A[0].getNumCols() <= A[0].getNumRows())
			&& (Q[0].getNumRows() == A[0].getNumCols()) && (Q[0].getNumCols() == A[0].getNumRows());
}

template <class T, int algo>
bool core::apparatus<T, algo>::isQInvertible(const tFuncMatrix& A, int r) const
{
	// TODO - What about the matrix Q?
	// TODO - the product Q*A must be invertible
	return (r == A.getNumCols()) && (A.getNumCols() <= A.getNumRows());
}

template <class T, int algo>
bool core::apparatus<T, algo>::isQInvertible(const tMatrixDiscretes& A, int r) const
{
	// TODO - What about the matrix Q?
	// TODO - the product Q*A must be invertible
	assert (m_di.K >= 0);
	return (r == A[0].getNumCols()) && (A[0].getNumCols() <= A[0].getNumRows());
}

template <class T, int algo>
bool core::apparatus<T, algo>::getQInverse(const tFuncMatrix& A, tMatrixDiscretes& binv)
{
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	int r = this->getRank(discretes);
	return this->getQInverse(discretes, r, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getQInverse(const tFuncMatrix& A, int r, tMatrixDiscretes& binv)
{
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	return this->getQInverse(discretes, r, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getQInverse(const tMatrixDiscretes& A, tMatrixDiscretes& binv)
{
	return this->getQInverse(A, this->getRank(A), binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getQInverse(const tMatrixDiscretes& A, int r, tMatrixDiscretes& binv)
{
	tMatrixDiscretes Q;
	this->chooseQ(A, r, Q);
	return this->getQInverse(A, r, Q, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getQInverse(const tFuncMatrix& A, const tFuncMatrix& Q, tMatrixDiscretes& binv)
{
	return this->getQInverse(A, this->getRank(A), Q, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getQInverse(const tFuncMatrix& A, int r, const tFuncMatrix& Q, tMatrixDiscretes& binv)
{
	if (!this->isQInvertible(A, Q, r))
		throw algoException("(Q)-Inversion: The matrix is not (Q)-invertible");
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	tMatrixDiscretes bDiscs;
	this->applyDiffTrans(Q, bDiscs);
	return this->getQInverse(discretes, r, bDiscs, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getQInverse(const tMatrixDiscretes& A, int r, const tMatrixDiscretes& Q, tMatrixDiscretes& binv)
{
	if (!this->isQInvertible(A, Q, r))
		throw algoException("(Q)-Inversion: The matrix is not (Q)-invertible");
	
	//std::cout << "A:\n";
	//std::copy(A.begin(), A.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//std::cout << "Q:\n";
	//std::copy(Q.begin(), Q.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//
	// Calculate the product of Q*A
	//
	tMatrixDiscretes QA;
	this->multDiscretes(Q, A, QA);
	//std::cout << "QA:\n";
	//std::copy(QA.begin(), QA.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//
	// Calculate the inverse to it
	//
	tMatrixDiscretes QA1;
	this->getInverse(QA, r, QA1);
	//std::cout << "QA1:\n";
	//std::copy(QA1.begin(), QA1.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	////
	// Calculate the (Q)-inverse
	//
	return this->multDiscretes(QA1, Q, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getInverseLowerTriangular(const tFuncMatrix& A, int r, tMatrixDiscretes& inv)
{
	if (!this->isInvertible(A, r))
		throw algoException("Inversion: the matrix is not invertible");
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	return this->getInverseLowerTriangular(discretes, r, inv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::isInvertible(const tFuncMatrix& A, int r) const
{
	if (!isSquare(A))
		return false;
	if (r != A.getNumRows())
		return false;
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::isInvertible(const tMatrixDiscretes& A, int r) const
{
	if (!isSquare(A))
		return false;
	assert (m_di.K >= 0);
	if (r != A[0].getNumRows())
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
	//for (int i = 1; i <= n; ++i)
	//{
	//	if (L[0][i][i] == 0)
	//		throw algoException("Inversion: the L Matrix is not invertible");
	//}
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
	//std::cout << "L1:\n";
	//std::copy(L1.begin(), L1.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::getInverseUpperTriangular(const tFuncMatrix& A, int r, tMatrixDiscretes& inv)
{
	if (!this->isInvertible(A, r))
		throw algoException("Inversion: the matrix is not invertible");
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
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
	//std::cout << "U1:\n";
	//std::copy(U1.begin(), U1.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::impl_getInverse(const tMatrixDiscretes& A, int r, tMatrixDiscretes& inv, eTriangularType et)
{
	switch (et)
	{
	case eLowerTriangular:
		return this->getInverseLowerTriangular(A, r, inv);
	break;
	case eUpperTriangular:
		return this->getInverseUpperTriangular(A, r, inv);
	break;
	default:
	break;
	}
	tMatrixDiscretes L;
	tMatrixDiscretes U;
	this->getLU(A, r, L, U);
	this->checkVal(L, 0);
	this->checkVal(U, 0);
	//std::cout << "L:\n";
	//std::copy(L.begin(), L.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//std::cout << "U:\n";
	//std::copy(U.begin(), U.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
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
bool core::apparatus<T, algo>::getInverse(const tMatrixDiscretes& A, int r, tMatrixDiscretes& inv)
{
	return this->impl_getInverse(A, r, inv, eNoneTriangular);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getInverse(const tFuncMatrix& A, tMatrixDiscretes& inv)
{
	return this->getInverse(A, getRank(A), inv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getInverse(const tFuncMatrix& A, int r, tMatrixDiscretes& inv)
{
	if (!this->isInvertible(A, r))
		throw algoException("Inversion: the matrix is not invertible");
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
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

	//std::cout << "After first, L:\n";
	//std::copy(L.begin(), L.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//std::cout << "After first, U:\n";
	//std::copy(U.begin(), U.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));

	for (int k = 0; k < K; ++k)
	{
		int li = 2;
		int ui = 2;
		while (li <= r && ui <= r)
		{
			// L
			for (int p = 2; p <= li; ++p)
			{
				double sum = 0;
				for (int j = 1; j <= p-1; ++j)
				{
					for (int l = 0; l <= k; ++l)
					{
						sum += L[l][li][j]*U[k-l][j][p];
					}
				}
				L[k][li][p] = discretes[k][li][p] - sum;
			}

			// U
			for (int p = ui; p <= r; ++p)
			{
				double sum1 = 0;
				for (int j = 1; j <= ui-1; ++j)
					for (int g = 0; g <= k; ++g)
						sum1 += L[g][ui][j]*U[k-g][j][p];

				double sum2 = 0;
				for (int l = 1; l <= k; ++l)
				{
					sum2 += U[k-l][ui][p]*L[l][ui][ui];
				}
				U[k][ui][p] = (discretes[k][ui][p] - sum1 - sum2)/L[0][ui][ui];
			}
			++li;
			++ui;
		}
	}

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
	//			L[k][i][p] = discretes[k][i][p] - sum;
	//		}
	//	}
	//}

	////std::cout << "After second, L:\n";
	////std::copy(L.begin(), L.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));

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
	//			U[k][i][p] = (discretes[k][i][p] - sum1 - sum2)/L[0][i][i];
	//		}
	//	}
	//}

	//std::cout << "L:\n";
	//std::copy(L.begin(), L.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//std::cout << "U:\n";
	//std::copy(U.begin(), U.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::getLU(const tFuncMatrix& A, int r, tMatrixDiscretes& L, tMatrixDiscretes& U)
{
	if (!this->isSquare(A))
		throw matrixException("Matrix is not square!"); 
	tMatrixDiscretes discretes;
	this->applyDiffTrans(A, discretes);
	return this->getLU(discretes, r, L, U);
}

template <class T, int algo>
bool core::apparatus<T, algo>::getLU(const tFuncMatrix& A, tMatrixDiscretes& L, tMatrixDiscretes& U)
{
	if (!this->isSquare(A))
		throw matrixException("Matrix is not square!"); 
	return this->getLU(A, this->getRank(A), L, U);
}

template <class T, int algo>
bool core::apparatus<T, algo>::applyDiffTrans(const tFuncMatrix& theMatrix, tMatrixDiscretes& discs) const
{
	std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();

	const int m = theMatrix.getNumRows(); 
	const int n = theMatrix.getNumCols();
	discs.resize(m_di.K+1, tMatrixDiscrete(m, n, 0));
	tFuncMatrix A = theMatrix;
	long long kfac = 1;
	for (int k = 0; k <= m_di.K; ++k)
	{
		kfac *= (k > 0 ? k : 1);
		for (int i = 1; i <= m; ++i)
		{
			for (int j = 1; j <= n; ++j)
			{
				A[i][j] = A[i][j]->derivative(k > 0 ? 1 : 0);
				discs[k][i][j] = std::pow(m_di.H, k)*(*A[i][j])(m_di.tv)/(T)kfac;
			}
		}
	}
	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count();
	std::cout << "ApplyDiffTrans duration = " << duration << std::endl;


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
bool core::apparatus<T, algo>::powerDiscretes(const tMatrixDiscretes& x, int n, tMatrixDiscretes& out)
{
	if (n == 0)
	{
		out.resize(m_di.K+1, tMatrixDiscrete(x[0].getNumRows(), x[0].getNumCols(), 0));
		this->makeIdentity(out[0]);
		return true;
	}
	if (n == 1)
	{
		out = x;
		return true;
	}
	if (n == 2)
	{
		this->multDiscretes(x, x, out);
		return true;
	}
	tMatrixDiscretes xn1;
	this->powerDiscretes(x, n-1, xn1);
	return this->multDiscretes(xn1, x, out);
}

template <class T, int algo>
bool core::apparatus<T, algo>::powerDiscretes(const tScalarDiscretes& x, int n, tScalarDiscretes& out)
{
	if (n == 0)
	{
		out.resize(m_di.K+1, T(0));
		out[0] = 1;
		return true;
	}
	if (n == 1)
	{
		out = x;
		return true;
	}
	if (n == 2)
	{
		this->multDiscretes(x, x, out);
		return true;
	}
	tScalarDiscretes xn1;
	this->powerDiscretes(x, n-1, xn1);
	return this->multDiscretes(xn1, x, out);
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
bool core::apparatus<T, algo>::multDiscretes(const tScalarDiscretes& x, const tMatrixDiscretes& y, tMatrixDiscretes& out)
{
	assert (m_di.K+1 == x.size());
	assert (m_di.K+1 == y.size());
	assert (m_di.K >= 0);
	out.resize(m_di.K+1, tMatrixDiscrete(y[0].getNumRows(), y[0].getNumCols(), 0));
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
bool core::apparatus<T, algo>::multDiscretes(const tScalarDiscretes& x, const tScalarDiscretes& y, tScalarDiscretes& out)
{
	assert (m_di.K+1 == x.size());
	assert (m_di.K+1 == y.size());
	assert (m_di.K >= 0);
	out.resize(m_di.K+1, T(0));
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
bool core::apparatus<T, algo>::inverseDiscrete(const tScalarDiscretes& x, tScalarDiscretes& out)
{
	out.resize(m_di.K+1, T(0));
	if (x[0] == 0)
		throw algoException("Discrete-Inversion: The discrete cannot be inverted");
	out[0] = 1./x[0];
	for (int k = 1; k <= m_di.K; ++k)
	{
		for (int l = 1; l <= k; ++l)
			out[k] += out[k-l]*x[l];
		out[k] = -out[k]/x[0];
	}
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::restoreTaylorSingle(const tMatrixDiscretes& theDiscretes, tFuncMatrix& out, int K)
{
	//assert (m_di.K+1 == theDiscretes.size());
	assert (K >= 0);
	const int m = out.getNumRows();
	const int n = out.getNumCols();
	for (int i = 1; i <= m; ++i)
	{
		for (int j = 1; j <= n; ++j)
		{
			out[i][j] = tFunctionPtr(new const_function<tArgType>(theDiscretes[0][i][j]));
			for (int k = 1; k <= K; ++k)
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

template <class T, int algo>
bool core::apparatus<T, algo>::restorePade(const tMatrixDiscretes& theDiscretes, int m, int n, tFuncMatrix& out, int K)
{
	assert (m <= n);
	assert (m+n <= K);
	for (int i = 1; i <= theDiscretes[0].getNumRows(); ++i)
	{
		for (int j = 1; j <= theDiscretes[0].getNumCols(); ++j)
		{
			std::cout << "i = " << i << ", j = " << j << std::endl;
			tMatrixDiscretes d(m_di.K+1, tMatrixDiscrete(n, n, 0));
			int dk = m;
			for (int di = 1; di <= n; ++di)
			{
				dk = m + di-1;
				for (int dj = 1; dj <= n; ++dj)
				{
					d[0][di][dj] = theDiscretes[dk--][i][j];
				}
			}
			std::cout << "Got the left for b-calculation\n";
			std::cout << "D: " << d[0] << std::endl;
			if (this->is_equal(d[0], 0))
			{
				out[i][j] = tFunctionPtr(new const_function<T>(0));
				continue;
			}
			// TODO
			//continue;
			tMatrixDiscretes inv;
			//this->getBQInverse(d, inv);
			tMatrixDiscrete C(d[0].getNumRows(), d[0].getNumCols());
			tMatrixDiscrete F(d[0].getNumRows(), d[0].getNumCols());
			(void) C;
			(void) F;
			int r = 0;
			this->getRankDecomposition(d[0], r, C, F);
			std::cout<< "r = " << r << std::endl;
			if (r == n)
				this->getInverse(d, r, inv);
			else 
				this->getBQInverse(d, inv);
			//std::cout << "R = " << r << std::endl;
			//if (this->isBInvertible(d, r))
			//{
			//	std::cout << "Getting B-inverse!\n";
			//	this->getBInverse(d, r, inv);
			//}
			//else if (this->isQInvertible(d, r))
			//{
			//	std::cout << "Getting Q-inverse!\n";
			//	this->getQInverse(d, r, inv);
			//}
			//else
			//{
			//	std::cout << "Getting B,Q-inverse!\n";
			//	this->getBQInverse(d, inv);
			//}
			//this->getInverse(d, n, inv);
			//this->getBQInverse(d, inv);
			//this->getBQInverse(d, inv);
			//int r = this->getRank(d);
			//std::cout << "rank = " << r << std::endl;
			//if (this->isInvertible(d, r))
			//{
			//try
			//{
			//	this->getInverse(d, n, inv);
			//}
			//catch (...)
			//{
			//	out[i][j] = tFunctionPtr(new const_function<T>(0));
			//	continue;
			//}
			//}
			//else if (this->isBInvertible(d, r))
			//	this->getBInverse(d, r, inv);
			//else if (this->isQInvertible(d, r))
			//	this->getQInverse(d, r, inv);
			//else if (this->isBQInvertible(d, r))
			//	this->getBQInverse(d, inv);
			std::cout << "Got the inverse = " << inv[0] << "\n";
			tMatrixDiscrete right(n, 1, 0);
			for (int ri = 1; ri <= n; ++ri)
				right[ri][1] = -theDiscretes[m+ri][i][j];
			std::cout << "Got the right:\n" << right << std::endl;
			tMatrixDiscrete b = inv[0]*right;

			tMatrixDiscrete left(m+1, m+1, 0);
			for (int li = 1; li <= m+1; ++li)
			{
				left[li][li] = 1;
				int bi = li-1;
				for (int lj = 1; lj < li; ++lj)
				{
					left[li][lj] = b[bi--][1];
				}
			}
			tMatrixDiscrete aright(m+1, 1, 0);
			for (int ai = 1; ai <= m+1; ++ai)
			{
				aright[ai][1] = theDiscretes[ai-1][i][j];
			}
			std::cout << "aright:\n" << aright << std::endl;
			tMatrixDiscrete a = left*aright;
			std::cout << "A:\n" << a << std::endl;
			std::cout << "B:\n" << b << std::endl;
			tFunctionPtr p(new subtract<T>(
								tFunctionPtr(new parameter<T>()),
								tFunctionPtr(new const_function<T>(m_di.tv))));
			tFunctionPtr ba = p;
			if (m_di.H != 1 && m_di.H != 0)
				ba = tFunctionPtr(new divide<T>(p, tFunctionPtr(new const_function<T>(m_di.H))));
					
			tFunctionPtr up(new const_function<T>(a[1][1]));
			for (int i = 2; i <= m+1; ++i)
			{
				tFunctionPtr mul(new multiply<T, eMultNum>(a[i][1], 
									tFunctionPtr(new power<T>(ba, i-1))));
				up = tFunctionPtr(new add<T>(up, mul));
			}
			tFunctionPtr down(new const_function<T>(1));

			for (int i = 1; i <= n; ++i)
			{
				tFunctionPtr mul(new multiply<T, eMultNum>(b[i][1], 
									tFunctionPtr(new power<T>(ba, i))));
				down = tFunctionPtr(new add<T>(down, mul));
			}
			out[i][j] = tFunctionPtr(new divide<T>(up, down));
		}
	}
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::getRankDecomposition(const tMatrixDiscrete& A, int& r, tMatrixDiscrete& C, tMatrixDiscrete& F)
{
	r = 0;
	const int m = A.getNumRows();
	const int n = A.getNumCols();

	int limit = std::min(m, n);
	tMatrixDiscrete B = A;
	for (int i = 1; i <= limit; ++i)
	{
		//int k = i;
		//for (int p = i+1; p <= limit; ++p)
		//{
		//	if (B[p][i] > B[k][i])
		//		k = p;
		//}
		//B.swap(i, k);
		//std::cout << "After swap: " << B << std::endl;
		//if (!this->is_equal(B[i][i], 0))
		//{
		//	for (int j = i+1; j <= n; ++j)
		//		B[i][j] = B[i][j] / B[i][i];
		//	B[i][i] = 1;
		//	++r;
		//}
		//else
		//	B[i][i] = 0;
		if (this->is_equal(B[i][i], 0))
		{
			B[i][i] = 0;
			continue;
		}

		//for (int p = 1; p <= limit; ++p)
		//{
		//	if (i == p)
		//		continue;
		for (int p = i+1; p <= limit; ++p)
		{
			if (this->is_equal(B[p][i], 0))
			{
				B[p][i] = 0;
				continue;
			}
			for (int j = i+1; j <= n; ++j)
			{
				B[p][j] -= (B[p][i]*B[i][j])/B[i][i];
			}
			B[p][i] = 0;
		}
		//std::cout << "B: " <<  B << std::endl;
	}
	std::cout << "B:" << B << std::endl;
	for (int i = 1; i <= limit; ++i)
	{
		for (int j = 1; j <= n; ++j)
		{
			if (!this->is_equal(B[i][j], 0))
			{
				++r;
				break;
			}
		}
	}

	return true;
}




#endif // __APPARATUS_H__
