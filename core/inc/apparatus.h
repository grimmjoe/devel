#ifndef __APPARATUS_H__
#define __APPARATUS_H__

#include "function.h"
#include "matrix.h"
#include <vector>

namespace core
{
	namespace app
	{
		template <class T>
		class apparatus
		{
		public:
			typedef T tArgType;
			typedef typename function<T>::tFunction tFunction;
			typedef typename function<T>::tFunctionPtr tFunctionPtr;
			typedef core::matrix<tFunctionPtr> tFuncMatrix;
			typedef core::matrix<T> tDiscrete;
			typedef std::vector<tDiscrete> tDiscretes;

			//
			/// Get LU decomposition of the matrix
			//
			bool getLU(const tFuncMatrix& sourceMatrix, int K, tDiscrete& l, tDiscrete& u);
		};
	} // namespace app
} // namespace core

template <class T>
bool core::apparatus<T>::getLU(const tFuncMatrix& sourceMatrix, int K, tDiscrete& l, tDiscrete& u)
{
	const int m = sourceMatrix.getNumRows();
	const int n = sourceMatrix.getNumCols();
}

#endif // __APPARATUS_H__
