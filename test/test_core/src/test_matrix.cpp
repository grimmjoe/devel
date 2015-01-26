#include "matrix.h"
#include <iostream>

int main()
{
	std::cout << "Making the matrix" << std::endl;
	core::matrix<int> m(3, 4, 5);
	std::cout << "Made the matrix, now printing" << std::endl;
	//m.print();
	std::cout << m;
	//for (int i = 1; i <= m.getNumRows(); ++i)
	//{
	//	for (int j = 1; j <= m.getNumCols(); ++j)
	//		std::cout << m[i][j] << " ";
	//	std::cout << std::endl;
	//}
	std::cout << "Printed\n";

	return 0;
}
