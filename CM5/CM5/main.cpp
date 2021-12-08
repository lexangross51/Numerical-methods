#include "matrix.h"
#include <windows.h>

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	matrix m("matrix.txt");

	m.eigenval_max();
	m.eigenval_min();

	return 0;
}