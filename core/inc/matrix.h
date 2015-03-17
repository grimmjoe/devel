#ifndef MATRIX_H__
#define MATRIX_H__

#include <vector>
#include <algorithm>
#include <ostream>
#include <iterator>
#include <cassert>
#include <memory>

#include "comparator.h"

namespace core
{

	class matrixException : public std::exception
	{
		std::string mErr;
	public:
		matrixException(const std::string& err)
			: mErr(err)
		{}
		const char* what() const noexcept
		{
			return mErr.c_str();
		}
	};

	template <class E>
	class matrix
	{
	public:
		template <class T, class Alloc = std::allocator<T> >
		class numVector : public std::vector<T, Alloc>
		{
		public:
			typedef std::vector<T, Alloc> tBase;
			using typename tBase::allocator_type;
			using typename tBase::pointer;
			using typename tBase::const_pointer;
			using typename tBase::reference;
			using typename tBase::const_reference;
			using typename tBase::iterator;
			using typename tBase::const_iterator;
			using typename tBase::reverse_iterator;
			using typename tBase::const_reverse_iterator;
			using typename tBase::value_type;
			using typename tBase::size_type;
			using typename tBase::difference_type;
			explicit numVector(const allocator_type& alloc = allocator_type())
				: tBase(alloc)
			{}
			explicit numVector(size_type n)
				: tBase(n)
			{}
			explicit numVector(size_type n, const value_type& val,
					const allocator_type& alloc = allocator_type())
				: tBase(n, val, alloc)
			{}
			template <class InputIterator>
			explicit numVector(InputIterator first, InputIterator last,
					const allocator_type& alloc = allocator_type())
				: tBase(first, last, alloc)
			{}
			numVector(const numVector& nv)
				: tBase(nv)
			{}
			numVector(const numVector& nv, const allocator_type& alloc)
				: tBase(nv, alloc)
			{}
			numVector(numVector&& nv)
				: tBase(nv)
			{}
			numVector(numVector&& nv, const allocator_type& alloc)
				: tBase(nv, alloc)
			{}
			numVector(std::initializer_list<value_type> il,
				const allocator_type& alloc = allocator_type())
				: tBase(il, alloc)
			{}

			reference operator[](size_type i)
			{
				return tBase::operator[](i-1);
			}
			const_reference operator[](size_type i) const
			{
				return tBase::operator[](i-1);
			}

			reference at(size_type i)
			{
				return tBase::at(i-1);
			}

			const_reference at(size_type i) const
			{
				return tBase::at(i-1);
			}

		};
		typedef numVector<E> tRow;
		typedef numVector<tRow> tMatrix;
	protected:
		tMatrix mMatrix;
		int mNumRows;
		int mNumCols;

		std::shared_ptr<comparator<E> > m_comparator;
	public:
		matrix(int m, int n, const E& e = E())
			: mNumRows(m)
			, mNumCols(n)
			, m_comparator(nullptr)
		{
			mMatrix.resize(mNumRows);
			for (int i = 1; i <= mNumRows; ++i)
			{
				mMatrix[i].resize(mNumCols, e);
			}
		}

		matrix(const tRow& r)
			: mNumRows(1)
			, mNumCols(r.size())
			, m_comparator(nullptr)
		{
			mMatrix.resize(mNumRows);
			mMatrix[1].reserve(mNumCols);
			std::copy(r.begin(), r.end(), std::back_inserter(mMatrix[1]));
		}

		typedef matrix<E> tMatrixRow;
		typedef matrix<E> tMatrixColumn;

		tMatrixRow getRow(int i) const
		{
			return tMatrixRow(mMatrix[i]);
		}

		tMatrixColumn getColumn(int j) const
		{
			assert (j >= 1 && j <= mNumCols);
			tMatrixColumn mc(mNumRows, 1);
			int i = 1;
			std::for_each(mMatrix.begin(), mMatrix.end(),
				[&](const tRow& r)
				{
					mc[i++][1]=r[j];
				}
				);
			return mc;
		}

		void setComparator(std::shared_ptr<comparator<E> > cm)
		{
			m_comparator = cm;
		}

		const std::shared_ptr<comparator<E> >& getComparator() const
		{
			return m_comparator;
		}

		int getNumRows() const
		{
			return mNumRows;
		}

		int getNumCols() const
		{
			return mNumCols;
		}

		//
		/// Get the transpose of the matrix
		//
		/*!
		 TODO - This is not efficient at all.
		 I can get use of the old matrix by wrapping it up and redirecting
		 i,j -> j,i
		 numRows->numCols
		 numCols->numRows
		*/
		matrix<E> getTranspose() const
		{
			matrix<E> m(this->getNumCols(), this->getNumRows());
			for (int i = 1; i <= this->getNumRows(); ++i)
				for (int j = 1; j <= this->getNumCols(); ++j)
					m[j][i]=(*this)[i][j];
			return m;
		}

		void set(const E& e)
		{
			for (int i = 1; i <= this->getNumRows(); ++i)
				for (int j = 1; j <= this->getNumCols(); ++j)
					(*this)[i][j] = e;
		}

		tRow& operator[](int i)
		{
			return mMatrix[i];
		}
		const tRow& operator[](int i) const
		{
			return mMatrix[i];
		}

		tRow& at(int i)
		{
			return mMatrix.at(i);
		}
		const tRow& at(int i) const
		{
			return mMatrix.at(i);
		}

		friend std::ostream& operator<<(std::ostream& out, const matrix<E>& m)
		{

			std::for_each(m.mMatrix.begin(), m.mMatrix.end(),
				[&](const typename matrix<E>::tRow& r)
				{
					std::copy(r.begin(), r.end(),
						std::ostream_iterator<E>(out, " "));
					out << std::endl;
				}
				);
			return out;
		}

		matrix<E>& operator +=(const matrix<E>& mat)
		{
			if (mat.getNumRows() != this->getNumRows() ||
				mat.getNumCols() != this->getNumCols())
				throw matrixException("Matrices are of inaddible boundaries");
			for (int i = 1; i <= this->getNumRows(); ++i)
				for (int j = 1; j <= this->getNumCols(); ++j)
					(*this)[i][j] += mat[i][j];
			return *this;
		}

		matrix<E>& operator -=(const matrix<E>& mat)
		{
			if (mat.getNumRows() != this->getNumRows() ||
				mat.getNumCols() != this->getNumCols())
				throw matrixException("Matrices are of insubtractible boundaries");
			for (int i = 1; i <= this->getNumRows(); ++i)
				for (int j = 1; j <= this->getNumCols(); ++j)
					(*this)[i][j] -= mat[i][j];
			return *this;
		}

		matrix<E>& operator=(const matrix<E>& mat)
		{
			if (mat.getNumRows() != this->getNumRows() ||
				mat.getNumCols() != this->getNumCols())
				throw matrixException("Matrices are of incorrect boundaries");
			for (int i = 1; i <= this->getNumRows(); ++i)
				for (int j = 1; j <= this->getNumCols(); ++j)
					(*this)[i][j] = mat[i][j];
			return *this;
		}

		matrix<E>& operator *=(const E& e)
		{
			for (int i = 1; i <= this->getNumRows(); ++i)
				for (int j = 1; j <= this->getNumCols(); ++j)
					(*this)[i][j] *= e;
			return *this;
		}
	};

	template <class Elem>
	matrix<Elem> operator*(const matrix<Elem>& mat1, const matrix<Elem>& mat2)
	{
		int m1 = mat1.getNumRows();
		int n1 = mat1.getNumCols();
		int m2 = mat2.getNumRows();
		int n2 = mat2.getNumCols();
		if (n1 != m2)
			throw matrixException("Matrix multiplication not possible, boundary error");
		matrix<Elem> result(m1, n2, Elem());
		for (int i = 1; i <= m1; ++i)
		{
			for (int j = 1; j <= n2; ++j)
			{
				Elem sum = Elem(0);
				for (int k = 1; k <= n1; ++k)
				{
					sum += mat1[i][k]*mat2[k][j];
				}
				result[i][j] = sum;
			}
		}
		return result;
	}

	template <class Elem>
	matrix<Elem> operator*(const matrix<Elem>& mat1, const Elem& e)
	{
		const int m = mat1.getNumRows();
		const int n = mat1.getNumCols();
		matrix<Elem> result(m, n, 0);
		for (int i = 1; i <= m; ++i)
		{
			for (int j = 1; j <= n; ++j)
			{
				result[i][j] = mat1[i][j]*e;
			}
		}
		return result;
	}


	template <class Elem>
	matrix<Elem> operator*(const Elem& e, const matrix<Elem>& mat1)
	{
		const int m = mat1.getNumRows();
		const int n = mat1.getNumCols();
		matrix<Elem> result(m, n, 0);
		for (int i = 1; i <= m; ++i)
		{
			for (int j = 1; j <= n; ++j)
			{
				result[i][j] = mat1[i][j]*e;
			}
		}
		return result;
	}

	template <class Elem>
	bool operator==(const matrix<Elem>& mat1, const matrix<Elem>& mat2)
	{
		int m1 = mat1.getNumRows();
		int n1 = mat1.getNumCols();
		int m2 = mat2.getNumRows();
		int n2 = mat2.getNumCols();
		if (m1 != m2 || n1 != n2)
			return false;
		std::shared_ptr<comparator<Elem> > theComp = mat1.getComparator();
		for (int i = 1; i <= m1; ++i)
		{
			for (int j = 1; j <= n1; ++j)
			{
				if ((theComp != nullptr 
						? !theComp->is_equal(mat1[i][j], mat2[i][j])
						: mat1[i][j] != mat2[i][j]))
					return false;
			}
		}
		return true;
	}

	template <class Elem>
	bool operator!=(const matrix<Elem>& mat1, const matrix<Elem>& mat2)
	{
		return ! (mat1 == mat2);
	}

	template <class Elem>
	Elem trace(const matrix<Elem>& m)
	{
		if (m.getNumRows() != m.getNumCols())
			throw matrixException("Matrix is not square!"); 
		Elem sum = Elem(0);
		for (int i = 1; i <= m.getNumRows(); ++i)
			sum += m[i][i];
		return sum;
	}

} // namespace core

#endif // MATRIX_H__
