#ifndef MATRIX_H__
#define MATRIX_H__

#include <vector>
#include <algorithm>
#include <ostream>
#include <iterator>

namespace core
{
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
	public:
		matrix(int m, int n, const E& e = E())
			: mNumRows(m)
			, mNumCols(n)
		{
			mMatrix.resize(mNumRows);
			for (int i = 1; i <= mNumRows; ++i)
			{
				mMatrix[i].resize(mNumCols, e);
			}
		}

		int getNumRows() const
		{
			return mNumRows;
		}

		int getNumCols() const
		{
			return mNumCols;
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

		template <class Elem>
		friend std::ostream& operator<<(std::ostream& out, const matrix<Elem>& m)
		{
			std::for_each(m.mMatrix.begin(), m.mMatrix.end(),
				[&](const typename matrix<Elem>::tRow& r)
				{
					std::copy(r.begin(), r.end(),
						std::ostream_iterator<Elem>(out, " "));
					out << std::endl;
				}
				);
			return out;
		}
	};

} // namespace core

#endif // MATRIX_H__
