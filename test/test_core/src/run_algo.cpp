
#include "apparatus.h"
#include "arithmetic.h"
#include "function.h"
#include "matrix.h"
#include "parser.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <iterator>

#include <chrono>

#include "tbb/atomic.h"

enum eMulti
{
	eMultiNone = 0,
	eMultiBasic = 1,
	eMultiParallel = 2
};

template <class T>
class garsevan : public core::function<T>
{
public:
	typedef std::shared_ptr<core::function<T> > tFunctionPtr;
private:
	T begin;
	T end;
	tFunctionPtr m_base;
public:
	garsevan(T b, T e, tFunctionPtr mb)
		: begin(b)
		, end(e)
		, m_base(mb)
	{
	}

	virtual T operator()(const T& t)
	{
		if ((t > end) || (t < begin))
			return 0;
		if ((t == begin) || (t == end))
			return 0.5*(*m_base)(t);

		return (*m_base)(t);
	}
	virtual std::string toString() const
	{
		return "";
	}
	virtual tFunctionPtr derivative(int num)
	{
		return tFunctionPtr(nullptr);
	}
	virtual void optimize()
	{
	}
	virtual tFunctionPtr clone()
	{
		return tFunctionPtr(new garsevan<T>(begin, end, m_base));
	}
};

template <class APP, class D>
void apply(const std::string& method, APP& app, const D& theMatrixDiscretes, D& theDiscretes)
{
	if (method == "b_inverse")
		app.getBInverse(theMatrixDiscretes, theDiscretes);
	else if (method == "q_inverse")
		app.getQInverse(theMatrixDiscretes, theDiscretes);
	else if (method == "bq_inverse")
		app.getBQInverse(theMatrixDiscretes, theDiscretes);
	else if (method == "drazin_recursive")
		app.getDrazinInverseRecursive(theMatrixDiscretes, theDiscretes);
	else if (method == "drazin_skeleton")
		app.getDrazinInverseSkeleton(theMatrixDiscretes, theDiscretes);
	else if (method == "drazin_canonical")
		app.getDrazinInverseCanonical(theMatrixDiscretes, theDiscretes);
}

template <class DI, class MAT, class FUNC>
class multiTask : public tbb::task
{
	DI di;
	MAT mat;
	int degree;
	FUNC func;
public:
	multiTask(DI d, MAT m, int dg, FUNC f)
		: di(d)
		, mat(m)
		, degree(dg)
		, func(f)
	{
	}

	tbb::task* execute()
	{
		func(di, mat, degree);
		return NULL;
	}
};

template <class DI, class MAT, class FUNC>
tbb::task* makeMultiTask(tbb::task* parent, DI di, MAT theMatrix, int dg, FUNC f)
{
	return new(parent->allocate_child()) multiTask<DI, MAT, FUNC>(di, theMatrix, dg, f);
}


template <int algo>
void run_algo_multi(const std::string& file, const std::string& method, eMulti em = eMultiBasic)
{
	typedef core::apparatus<double, algo> tApp;
	typedef typename tApp::tMatrixDiscrete tDiscrete;
	typedef typename tApp::tMatrixDiscretes tDiscretes;
	typedef typename tApp::tFunction tFunction;
	typedef typename tApp::tFunctionPtr tFunctionPtr;
	typedef typename tApp::tFuncMatrix tMatrix;
	typedef typename tApp::tDiffInfo tDiffInfo;


	std::ifstream fin(file);
	assert (fin.good());
	int m, n;
	fin >> m >> n;
	int k;
	int degree;
	double a;
	double b;
	double h;
	double epsilon;
	fin >> k >> degree >> a >> b >> h  >> epsilon;
	tMatrix theMatrix(m, n, tFunctionPtr(nullptr));
	core::parser<double> p;
	for (int i = 1; i <= m; ++i)
	{
		for (int j = 1; j <= n; ++j)
		{
			std::string s;
			fin >> s;
			theMatrix[i][j] = p.parse(s);
		}
	}
	tbb::atomic<long long> duration;
	duration = 0;
	std::chrono::high_resolution_clock::time_point tbegin = std::chrono::high_resolution_clock::now();
	if (em == eMultiBasic)
	{
		tDiffInfo di(1, 1, k);
		tApp app(di);
		if (epsilon != 0)
		{
			app.setComparator(std::shared_ptr<core::comparator<double> >(new core::comparator<double>(epsilon)));
		}
		tMatrix realOrig(n, m, tFunctionPtr(nullptr));
		for (double begin = a; begin <= b; begin += h)
		{
			//std::cout << "********************************\n";
			//std::cout << "[" << begin << ", " << begin+h << "]\n";
			//std::cout << "********************************\n";
			tDiffInfo di(begin, h, k);
			app.setDiffInfo(di);
			tDiscretes theMatrixDiscretes;
			app.applyDiffTrans(theMatrix, theMatrixDiscretes, degree);
			tDiscretes theDiscretes;
			apply(method, app, theMatrixDiscretes, theDiscretes);


			tMatrix origMatrix(n, m, tFunctionPtr(nullptr));
			app.restoreTaylorSingle(theDiscretes, origMatrix, k-1);
			if (realOrig[1][1] == nullptr)
			{
				for (int i = 1; i <= n; ++i)
					for (int j = 1; j <= m; ++j)
						realOrig[i][j] = tFunctionPtr(new garsevan<double>(begin-h, begin+h, origMatrix[i][j]));
			}
			else
			{
				for (int i = 1; i <= n; ++i)
					for (int j = 1; j <= m; ++j)
						realOrig[i][j] = tFunctionPtr(new core::add<double>(realOrig[i][j],
												tFunctionPtr(new garsevan<double>(begin, begin+h, origMatrix[i][j]))));
			}
		}
	}
	else
	{
		std::shared_ptr<core::comparator<double> > cmp(new core::comparator<double>(epsilon));
		(void) cmp;
		tbb::task* dummy = new(tbb::task::allocate_root()) tbb::empty_task;
		dummy->increment_ref_count();
		for (double begin = a; begin <= b; begin += h)
		{
			tDiffInfo di(begin, h, k);
			tbb::task* ch = makeMultiTask(dummy, di, theMatrix, degree, 				[&](tDiffInfo& d, tMatrix& theM, int dg)
				{
					tApp app(d);
					app.setComparator(cmp);
					tDiscretes theMatrixDiscretes;
					app.applyDiffTrans(theM, theMatrixDiscretes, dg);
					tDiscretes theDiscretes;
					apply(method, app, theMatrixDiscretes, theDiscretes);
				}
			);
			dummy->increment_ref_count();
			dummy->spawn(*ch);
		}
		dummy->wait_for_all();
		dummy->destroy(*dummy);
	}
	std::chrono::high_resolution_clock::time_point tend = std::chrono::high_resolution_clock::now();
	duration += std::chrono::duration_cast<std::chrono::microseconds>(tend-tbegin).count();
	std::cout << "Duration = " << duration << std::endl;
}

template <int algo>
void run_algo(const std::string& file, const std::string& method, eMulti em = eMultiNone)
{
	if (em != eMultiNone)
		return run_algo_multi<algo>(file, method, em);
	typedef core::apparatus<double, algo> tApp;
	typedef typename tApp::tMatrixDiscrete tDiscrete;
	typedef typename tApp::tMatrixDiscretes tDiscretes;
	typedef typename tApp::tFunction tFunction;
	typedef typename tApp::tFunctionPtr tFunctionPtr;
	typedef typename tApp::tFuncMatrix tMatrix;
	typedef typename tApp::tDiffInfo tDiffInfo;
	std::ifstream fin(file);
	assert (fin.good());
	int m, n;
	fin >> m >> n;
	int k;
	double tv;
	fin >> k >> tv;
	int degree;
	fin >> degree;
	double epsilon;
	fin >> epsilon;
	tDiffInfo di(tv, 1, k);
	tApp app(di);
	if (epsilon != 0)
	{
		app.setComparator(std::shared_ptr<core::comparator<double> >(new core::comparator<double>(epsilon)));
	}
	tDiscretes theDiscretes;
	tDiscretes theMatrixDiscretes;
	{
	tMatrix theMatrix(m, n, tFunctionPtr(nullptr));
	core::parser<double> p;
	for (int i = 1; i <= m; ++i)
	{
		for (int j = 1; j <= n; ++j)
		{
			std::string s;
			fin >> s;
			theMatrix[i][j] = p.parse(s);
		}
	}
	//std::cout << "The matrix:\n";
	//printFuncMatrix(theMatrix);
	app.applyDiffTrans(theMatrix, theMatrixDiscretes, degree);
	}
	//std::cout << "The matrix discretes:\n";
	//std::copy(theMatrixDiscretes.begin(), theMatrixDiscretes.end(), std::ostream_iterator<tDiscrete>(std::cout, "\n"));
	std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
	apply(method, app, theMatrixDiscretes, theDiscretes);
	//app.getDrazinInverseRecursive(theMatrixDiscretes, theDiscretes);
	//app.getDrazinInverseRecursive(theMatrix, theDiscretes);
	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count();
	std::cout << "Duration = " << duration << std::endl;
}

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		std::cerr << "Please specify the file, the algorithm and parallel/serial\n";
		return 1;
	}

	if (strcmp(argv[3], "basic") == 0)
	{
		run_algo<core::eAlgoBasic>(argv[1], argv[2]);
	}
	else if (strcmp(argv[3], "mult_parallel_for") == 0)
	{
		run_algo<core::eAlgoParallelMultParallelFor>(argv[1], argv[2]);
	}
	else if (strcmp(argv[3], "mult_parallel_task") == 0)
	{
		run_algo<core::eAlgoParallelMultTask>(argv[1], argv[2]);
	}
	else if (strcmp(argv[3], "multi_basic") == 0)
	{
		run_algo<core::eAlgoBasic>(argv[1], argv[2], eMultiBasic);
	}
	else if (strcmp(argv[3], "multi_parallel_basic") == 0)
	{
		run_algo<core::eAlgoBasic>(argv[1], argv[2], eMultiParallel);
	}
	else if (strcmp(argv[3], "multi_parallel_mult_parallel_for") == 0)
	{
		run_algo<core::eAlgoParallelMultParallelFor>(argv[1], argv[2], eMultiParallel);
	}
	else if (strcmp(argv[3], "multi_parallel_mult_parallel_task") == 0)
	{
		run_algo<core::eAlgoParallelMultTask>(argv[1], argv[2], eMultiParallel);
	}
	else if (strcmp(argv[3], "multi_basic_mult_parallel_for") == 0)
	{
		run_algo<core::eAlgoParallelMultParallelFor>(argv[1], argv[2], eMultiBasic);
	}
	else if (strcmp(argv[3], "multi_basic_mult_parallel_task") == 0)
	{
		run_algo<core::eAlgoParallelMultTask>(argv[1], argv[2], eMultiBasic);
	}
	else
	{
		std::cerr << "UNKNOWN type of algorithm will use Basic instead!\n";
		run_algo<core::eAlgoBasic>(argv[1], argv[2]);
	}
	
	return 0;
}
