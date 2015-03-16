#ifndef __APPARATUS_H__
#define __APPARATUS_H__

#include "function.h"
#include "constant_function.h"
#include "power_function.h"
#include "add.h"
#include "subtract.h"
#include "multiply.h"
#include "matrix.h"
#include "comparator.h"
#include <vector>
#include <cmath>
#include <memory>

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
		/// Restore the original @a out from the discretes @a theDiscretes with reverse single-point
		/// Taylor transformations
		//
		bool restoreTaylorSingle(const tMatrixDiscretes& theDiscretes, tFuncMatrix& out);

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
		/// Get the Drazin inverse
		//
		bool getDrazinInverseRecursive(const tFuncMatrix& A, tMatrixDiscretes& dinv);

		//
		/// Get the Drazin inverse
		//
		bool getDrazinInverseRecursive(const tMatrixDiscretes& A, tMatrixDiscretes& dinv);

		//
		/// Get the Drazin inverse
		//
		bool getDrazinInverseRecursive(const tFuncMatrix& A, int r, tMatrixDiscretes& dinv);

		//
		/// Get the Drazin inverse
		//
		bool getDrazinInverseRecursive(const tMatrixDiscretes& A, int r, tMatrixDiscretes& dinv);

	
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
	tMatrixDiscretes pm;
	return this->impl_getRank(discretes, pm);
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
			if (d[i].second == 0)
				continue;
			T scalar = this->getScalarProduct(aj, d[i].first);
			if (scalar == 0)
				continue;
			d[j].first -= (scalar/d[i].second)*d[i].first;
		}
		d[j].second = this->getNormNonSqrt(d[j].first);
		if (this->is_equal(d[j].second, 0))
		{
			d[j].second = 0;
			d[j].first.set(0);
		}
	}
	std::vector<int> theIndices(n, 0);
	int r = 0;
	int nn = n;
	for (int i = 1; i <= d.size(); ++i)
	{
		if (!is_equal(d[i].second, 0))
			theIndices[r++] = i;
		else
			theIndices[--nn] = i;
	}
	std::cout << "Orthogonalization was:\n";
	for (int i = 1; i <= d.size(); ++i)
	{
		std::cout << i << ":\n" << d[i].first << std::endl;
		std::cout << d[i].second << std::endl;
	}
	assert (r <= std::min(m, n));
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
	
	std::cout << "The rank is " << r << std::endl;
	std::cout << "The indices:\n";
	std::copy(theIndices.begin(), theIndices.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << "The permutation matrix is:\n";
	std::cout << permMatrix[0] << std::endl;

	return r;
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
	std::cout << "Here\n";
	assert (m_di.K >= 0);
	const int n = A[0].getNumRows();
	std::vector<tMatrixDiscretes> S(n+1, tMatrixDiscretes());
	S[0].resize(m_di.K+1, tMatrixDiscrete(n, n, 0));
	this->makeIdentity(S[0][0]);
	int j = 1;
	int p = 0;
	bool toStop = j > n;
	typedef std::vector<T> tScalarDiscretes;
	std::vector<tScalarDiscretes> betta(n, tScalarDiscretes(m_di.K+1, 0));
	tMatrixDiscrete I(n, n, 0);
	this->makeIdentity(I);
	std::cout << "Entering the loop\n";
	while (!toStop)
	{
		p = 0;
		tMatrixDiscretes AS;
		this->multDiscretes(A, S[j-1], AS);
		std::cout << "AS:\n";
		std::copy(AS.begin(), AS.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, " "));
		for (int k = 0; k <= m_di.K; ++k)
		{
			betta[n-j][k] = (-1./j) * trace(AS[k]);
			if (this->is_equal(betta[n-j][k], 0))
				betta[n-j][k] = 0;
		}
		std::cout << "Calculated bettas, now calculating S\n";
		std::cout << "betta[" << n-j << "]:\n";
		std::copy(betta[n-j].begin(), betta[n-j].end(), std::ostream_iterator<T>(std::cout, " "));
		this->multDiscretes(A, S[j-1], S[j]);
		for (int k = 0; k <= m_di.K; ++k)
		{
			if (!this->is_equal(betta[n-j][k], 0))
				S[j][k] += betta[n-j][k]*I;
			std::cout << "Was not null, checking the matrix\n";
			if (!this->is_equal(S[j][k], 0))
				p = -1;
		}
		std::cout << "Calculated S\n";
		std::cout << "p = " << p << std::endl;
		if (p != -1)
		{
			p = j;
			break;
		}
		++j;
		toStop = j > n;
	}
	std::cout << "p = " << p << std::endl;
	std::cout << "betta:\n";
	for (int i = 0; i < n; ++i)
	{
		std::cout << "betta[" << i << "]:\n";
		std::copy(betta[i].begin(), betta[i].end(), std::ostream_iterator<T>(std::cout, " "));
		std::cout << std::endl;
		std::cout << "S[" << i << "]:\n";
		std::copy(S[i].begin(), S[i].end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, " "));
		std::cout << std::endl;
	}
		std::cout << "S[" << n << "]:\n";
		std::copy(S[n].begin(), S[n].end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, " "));
		std::cout << std::endl;
	return true;
}

template <class T, int algo>
bool core::apparatus<T, algo>::isBQInvertible(const tFuncMatrix& A) const
{
	int r = this->getRank(A);
	return r < std::min(A.getNumRows(), A.getNumCols());
}

template <class T, int algo>
bool core::apparatus<T, algo>::isBQInvertible(const tFuncMatrix& A, int r) const
{
	return r < std::min(A.getNumRows(), A.getNumCols());
}


template <class T, int algo>
bool core::apparatus<T, algo>::isBQInvertible(const tMatrixDiscretes& A, int r) const
{
	assert (m_di.K >= 0);
	return r < std::min(A[0].getNumRows(), A[0].getNumCols());
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
		throw algoException("(Q)-Inversion: The matrix is not (B, Q)-invertible");
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

	tMatrixDiscretes Q;
	this->chooseQ(S, r, Q);
	std::cout << "Chose Q as:\n";
	std::copy(Q.begin(), Q.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes qinv;
	this->getQInverse(S, r, Q, qinv);
	std::cout << "Q-Inverse is:\n";
	std::copy(qinv.begin(), qinv.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes R;
	this->multDiscretes(qinv, A, R);
	std::cout << "R is:\n";
	std::copy(R.begin(), R.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes B;
	this->chooseB(R, r, B);
	std::cout << "Chose B as:\n";
	std::copy(B.begin(), B.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	tMatrixDiscretes rinv;
	this->getBInverse(R, r, B, rinv);
	std::cout << "B-Inverse is:\n";
	std::copy(rinv.begin(), rinv.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
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
		std::cout << "k = " << k << std::endl;
		std::cout << A[k] << std::endl;
		std::cout << a_mult[k] << std::endl;
		if (A[k] != a_mult[k])
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
		std::cout << "k = " << k << std::endl;
		std::cout << inv[k] << std::endl;
		std::cout << a_mult[k] << std::endl;
		if (inv[k] != a_mult[k])
		{
			std::cerr << "Discretes not equal at " << k << std::endl;
			return false;
		}
	}
	// Second
	tMatrixDiscretes a_ainv;
	this->multDiscretes(A, inv, a_ainv);
	for (int k = 0; k <= K; ++k)
	{
		std::cout << "k = " << k << std::endl;
		std::cout << a_ainv[k] << std::endl;
		std::cout << ainv_a[k] << std::endl;
		if (a_ainv[k] != ainv_a[k])
		{
			std::cerr << "Discretes not equal at " << k << std::endl;
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
	return (A.getNumRows() < A.getNumCols()) && (r == A.getNumRows())
			&& (B.getNumRows() == A.getNumCols()) && (B.getNumCols() == A.getNumRows());
}

template <class T, int algo>
bool core::apparatus<T, algo>::isBInvertible(const tMatrixDiscretes& A, const tMatrixDiscretes& B, int r) const
{
	// TODO - What about the matrix B?
	// TODO - the product A*B must be invertible
	assert (m_di.K >= 0);
	return (A[0].getNumRows() < A[0].getNumCols()) && (r == A[0].getNumRows())
			&& (B[0].getNumRows() == A[0].getNumCols()) && (B[0].getNumCols() == A[0].getNumRows());
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
	std::cout << "AB:\n";
	std::copy(AB.begin(), AB.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	//
	// Calculate the inverse to it
	//
	tMatrixDiscretes AB1;
	this->getInverse(AB, r, AB1);
	std::cout << "AB1:\n";
	std::copy(AB1.begin(), AB1.end(), std::ostream_iterator<tMatrixDiscrete>(std::cout, "\n"));
	////
	// Calculate the (B)-inverse
	//
	return this->multDiscretes(B, AB1, binv);
}

template <class T, int algo>
bool core::apparatus<T, algo>::chooseB(const tMatrixDiscretes& A, int r, tMatrixDiscretes& B)
{
	assert (m_di.K >= 0);
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
	return (r == A.getNumCols()) && (A.getNumCols() < A.getNumRows())
			&& (Q.getNumRows() == A.getNumCols()) && (Q.getNumCols() == A.getNumRows());
}

template <class T, int algo>
bool core::apparatus<T, algo>::isQInvertible(const tMatrixDiscretes& A, const tMatrixDiscretes& Q, int r) const
{
	// TODO - What about the matrix Q?
	// TODO - the product Q*A must be invertible
	assert (m_di.K >= 0);
	return (r == A[0].getNumCols()) && (A[0].getNumCols() < A[0].getNumRows())
			&& (Q[0].getNumRows() == A[0].getNumCols()) && (Q[0].getNumCols() == A[0].getNumRows());
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
