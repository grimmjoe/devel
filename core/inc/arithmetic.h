#ifndef __ARITHMETIC_H__
#define __ARITHMETIC_H__

#include "function.h"
#include "add.h"
#include "subtract.h"
#include "multiply.h"
#include "divide.h"
#include "matrix.h"
#include "constant_function.h"

namespace core
{
	template <class T>
	typename function<T>::tFunctionPtr operator+(typename function<T>::tFunctionPtr f1,
												 typename function<T>::tFunctionPtr f2)
	{
		return typename function<T>::tFunctionPtr(new add<T>(f1, f2));
	}

	template <class T>
	typename function<T>::tFunctionPtr operator-(typename function<T>::tFunctionPtr f1,
												 typename function<T>::tFunctionPtr f2)
	{
		return typename function<T>::tFunctionPtr(new subtract<T>(f1, f2));
	}

	template <class T>
	typename function<T>::tFunctionPtr operator*(typename function<T>::tFunctionPtr f1,
												 typename function<T>::tFunctionPtr f2)
	{
		return typename function<T>::tFunctionPtr(new multiply<T>(f1, f2));
	}

	template <class T>
	typename function<T>::tFunctionPtr operator/(typename function<T>::tFunctionPtr f1,
												 typename function<T>::tFunctionPtr f2)
	{
		return typename function<T>::tFunctionPtr(new divide<T>(f1, f2));
	}

	template <class T>
	std::shared_ptr<function<T> > operator-(std::shared_ptr<function<T> > f1,
										    std::shared_ptr<function<T> > f2)
	{
		return typename function<T>::tFunctionPtr(new subtract<T>(f1, f2));
	}

	template <class T>
	std::shared_ptr<function<T> > operator+(std::shared_ptr<function<T> > f1,
										    std::shared_ptr<function<T> > f2)
	{
		return typename function<T>::tFunctionPtr(new add<T>(f1, f2));
	}

	template <class T>
	std::shared_ptr<function<T> > operator*(std::shared_ptr<function<T> > f1,
										    std::shared_ptr<function<T> > f2)
	{
		return typename function<T>::tFunctionPtr(new multiply<T>(f1, f2));
	}

	template <class T>
	std::shared_ptr<function<T> > operator/(std::shared_ptr<function<T> > f1,
										    std::shared_ptr<function<T> > f2)
	{
		return typename function<T>::tFunctionPtr(new divide<T>(f1, f2));
	}

	template <class T>
	matrix<std::shared_ptr<function<T> > > operator*(const matrix<std::shared_ptr<function<T> > >& mat1, const matrix<std::shared_ptr<function<T> > >& mat2)
	{
		int m1 = mat1.getNumRows();
		int n1 = mat1.getNumCols();
		int m2 = mat2.getNumRows();
		int n2 = mat2.getNumCols();
		if (n1 != m2)
			throw matrixException("Matrix multiplication not possible, boundary error");
		typedef std::shared_ptr<function<T> > Elem;
		matrix<Elem> result(m1, n2, Elem(nullptr));
		for (int i = 1; i <= m1; ++i)
		{
			for (int j = 1; j <= n2; ++j)
			{
				std::shared_ptr<function<T> > sum(new const_function<T>(0));
				for (int k = 1; k <= n1; ++k)
				{
					//sum += mat1[i][k]*mat2[k][j];
					sum = Elem(new add<T>(sum, mat1[i][k]*mat2[k][j]));
				}
				result[i][j] = sum;
			}
		}
		return result;
	}
} // namespace core

#endif // __ARITHMETIC_H__
