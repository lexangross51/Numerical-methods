#include "snu.h"
#include <algorithm>

// Конструктор
snu::snu(const std::string directory)
{
	n = m = maxIter = eps1 = eps2 = 0;
	beta = 1;
	std::ifstream params(directory + "params.txt");

	if (params.is_open())
	{
		params >> m >> n >> maxIter >> eps1 >> eps2;
		params.close();

		x.resize(n);
		F.resize(m);
		A.resize(m);
		for (size_t i = 0; i < m; i++)
			A[i].resize(n);
	}

	std::ifstream initApprox(directory + "initApprox.txt");
	if (initApprox.is_open())
	{
		for (size_t i = 0; i < n; i++)
			initApprox >> x[i];
		initApprox.close();
	}
}

// Посчитать норму вектора
real snu::norm(Vector& v)
{
	real scalar = 0;
	for (size_t i = 0; i < v.size(); i++)
		scalar += v[i] * v[i];
	return sqrt(scalar);
}

// Посчитать значение каждой из функций
real snu::f(size_t funcNumber)
{
	switch (funcNumber)
	{
	case 0:
		return x[0] * x[0] + x[1] * x[1] - 1;
	case 1:
		return (x[0] - 2) * (x[0] - 2) + x[1] * x[1] - 1;
	}
}

// Получить вектор правой части
void snu::f()
{
	for (size_t i = 0; i < m; i++)
		F[i] = -f(i);
}

// Производная функции
// funcNumber - номер функции, которую будем дифференцировать
// varNumber - номер переменной, по которой будем дифференцировать
real snu::df(size_t funcNumber, size_t varNumber)
{
	switch (funcNumber)
	{
	case 0:
		switch (varNumber)
		{
		case 0:
			return 2 * x[0];
		case 1:
			return 2 * x[1];
		}
	case 1:
		switch (varNumber)
		{
		case 0:
			return 2 * (x[0] - 2);
		case 1:
			return 2 * x[1];
		}
	}
}

// Найти минимальные |F|
void snu::findMinF()
{
	for (size_t i = 0; i < m - n; i++)
	{
		size_t min = 0;
		for (size_t j = 1; j < m - i; j++)
		{
			if (abs(F[j]) < abs(F[min]))
				min = j;
		}
		minF[i] = min;
	}
	std::sort(minF.begin(), minF.end());
}

// Проверяем, есть ли уравнени с номером equationNum
// в списке уравнений, которые нужно исключить
bool snu::isInMinF(size_t equationNum)
{
	for (size_t i = 0; i < minF.size(); i++)
		if (equationNum == i)
			return true;
	return false;
}

// Построить матрицу Якоби
void snu::createJacobiMatrix()
{
	// Если число уравнений совпадает с числом неизвестных 
	if (m == n)
	{
		// то заполняем матрицу Якоби
		for (size_t i = 0; i < m; i++)
			for (size_t j = 0; j < n; j++)
				A[i][j] = df(i, j);
	}
	// Иначе нужно исключить m - n уравнений с минимальными
	// |F|
	else
	{
		for (size_t l = 0, i = 0, k = 0; l < m; l++)
		{
			if (l != minF[k])
			{
				for (size_t j = 0; j < n; j++)
					A[i][j] = df(i, j);
				F[i++] = F[l];
			}
			else
				k = k + 1 < m - n ? k + 1 : k;
		}
	}
}