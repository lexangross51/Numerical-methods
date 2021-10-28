#include <iostream>
#include <chrono>
#include "slau.h"

using namespace std;
using namespace chrono;

const string dense = "dense.txt";
const string directory = "C:\\Users\\lexan\\OneDrive\\Рабочий стол\\НГТУ\\3 курс\\Численные методы\\CM3\\CM3\\";

int main()
{
	//slau s(directory);
	size_t itersCount = 0;	// Получившееся число итераций 
	real residual = 0;		// Относительная невязка
	
	//size_t size = 0;		// Размер матрицы Гильберта



	for (size_t size = 1; size < 14; size++)
	{
		slau s(size, 10000, 1e-15, 0);
		s.createHilbertMatrix(size);
		auto start = high_resolution_clock::now();
		//s.MSGLU(&residual, &itersCount);
		//s.MSGdiag(&residual, &itersCount);
		//s.LOS(&residual, &itersCount);
		s.LOSLU(&residual, &itersCount);
		//s.LOSdiag(&residual, &itersCount);
		auto end = high_resolution_clock::now();
		auto ms = duration_cast<std::chrono::microseconds>(end - start);

		s.writeResult("result.txt");
		cout << "Iters = " << itersCount << "   Residual = " << residual
			<< "   Time = " << ms.count() << endl;
	}
	return 0;
}