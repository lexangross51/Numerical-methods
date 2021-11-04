﻿#include "snu.h"
#include <algorithm>
#include <iomanip>

using std::setw;
using std::setprecision;
using std::left;

// Конструктор
snu::snu(const std::string directory)
{
	n = m = maxIter = eps1 = eps2 = beta = 0;
	normF = currNormF = prevNormF = 0;
	std::ifstream params(directory + "params.txt");

	if (params.is_open())
	{
		params >> m >> n >> maxIter >> eps1 >> eps2;
		params.close();

		x.resize(n);
		F.resize(m);
		A.resize(m);
		deltaX.resize(n, 0);
		for (size_t i = 0; i < m; i++)
			A[i].resize(n);

		if (m > n)
			minF.resize(m - n);
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
real snu::f(tests test, size_t funcNumber)
{
	switch (test)
	{
	case tests::TEST1:
	{
		switch (funcNumber)
		{
		case 0:
			return (x[0] + 2) * (x[0] + 2) + x[1] * x[1] - 1;
		case 1:
			return (x[0] - 2) * (x[0] - 2) + x[1] * x[1] - 1;
		}
	}

	case tests::TEST2:
	{
		switch (funcNumber)
		{
		case 0:
			return x[0] * x[0] + x[1] * x[1] - 1;
		case 1:
			return (x[0] - 2) * (x[0] - 2) + x[1] * x[1] - 1;
		}
	}

	case tests::TEST3:
	{
		switch (funcNumber)
		{
		case 0:
			return (x[0] - 1) * (x[0] - 1) + x[1] * x[1] - 1;
		case 1:
			return (x[0] - 2) * (x[0] - 2) + x[1] * x[1] - 1;
		}
	}

	case tests::TEST4:
	{
		switch (funcNumber)
		{
		case 0:
			return x[0] * x[0] + x[1] * x[1] - 1;
		case 1:
			return (x[0] - 2) * (x[0] - 2) + x[1] * x[1] - 1;
		case 2:
			return x[0] - x[1] - 1;
		}
	}

	case tests::TEST5:
	{
		switch (funcNumber)
		{
		case 0:
			return x[0] - x[1] + 1;
		case 1:
			return x[0] + x[1];
		case 2:
			return 0.5 * x[0] + x[1];
		}
	}

	case tests::TEST6:
	{
		switch (funcNumber)
		{
		case 0:
			return sin(x[0]) - x[1];
		case 1:
			return x[0] - x[1];
		}
	}
	}	
}

// Получить вектор правой части
void snu::f(tests test)
{
	for (size_t i = 0; i < m; i++)
		F[i] = -f(test, i);
}

// Производная функции
// funcNumber - номер функции, которую будем дифференцировать
// varNumber - номер переменной, по которой будем дифференцировать
real snu::df(tests test, size_t funcNumber, size_t varNumber)
{
	switch(test)
	{
	case tests::TEST1:
	{
		switch (funcNumber)
		{
		case 0:
			switch (varNumber)
			{
			case 0:
				return 2 * (x[0] + 2);
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

	case tests::TEST2:
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

	case tests::TEST3:
	{
		switch (funcNumber)
		{
		case 0:
			switch (varNumber)
			{
			case 0:
				return 2 * (x[0] - 1);
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
	
	case tests::TEST4:
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
		case 2:
			switch (varNumber)
			{
			case 0:
				return 1;
			case 1:
				return -1;
			}
		}
	}

	case tests::TEST5:
	{
		switch (funcNumber)
		{
		case 0:
			switch (varNumber)
			{
			case 0:
				return 1;
			case 1:
				return -1;
			}
		case 1:
			switch (varNumber)
			{
			case 0:
				return 1;
			case 1:
				return 1;
			}
		case 2:
			switch (varNumber)
			{
			case 0:
				return 0.5;
			case 1:
				return 1;
			}
		}
	}

	case tests::TEST6:
	{
		switch (funcNumber)
		{
		case 0:
			switch (varNumber)
			{
			case 0:
				return cos(x[0]);
			case 1:
				return -1;
			}
		case 1:
			switch (varNumber)
			{
			case 0:
				return 1;
			case 1:
				return -1;
			}
		}
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

// Построить матрицу Якоби
void snu::createJacobiMatrix(tests test)
{
	// Если число уравнений совпадает с числом неизвестных 
	if (m == n)
	{
		// то заполняем матрицу Якоби
		for (size_t i = 0; i < m; i++)
			for (size_t j = 0; j < n; j++)
				A[i][j] = df(test, i, j);
	}
	// Иначе нужно исключить (m - n) уравнений с минимальными |F|
	else
	{
		size_t ii = 0, k = 0;
		for (size_t i = 0; i < m; i++)
		{
			if (i != minF[k])
			{
				for (size_t j = 0; j < n; j++)
					A[ii][j] = df(test, i, j);
				F[ii++] = F[i];
			}
			else
				k = k + 1 < m - n ? k + 1 : k;
		}
	}
}

// Метод Гаусса
void snu::gauss()
{	
	real verySmall = 1e-15;
	real max = 0;

	for (size_t k = 0; k < n; k++)
	{
		// Находим строку с максимальным элементом
		max = abs(A[k][k]);
		size_t ind = k;
		for (size_t i = k + 1; i < n; i++)
			if (abs(A[i][k]) > max)
			{
				max = abs(A[i][k]);
				ind = i;
			}

		// Меняем местами текущую строку со строкой, в которой
		// находится максимальный элемент
		std::swap(A[k], A[ind]);
		std::swap(F[k], F[ind]);

		// Нормализуем уравнения
		for (size_t i = k; i < n; i++)
		{
			real m = A[i][k];
			if (abs(m) < verySmall)
				throw exceptions::VERY_SMALL_DIAG;

			for (size_t j = 0; j < n; j++)
				A[i][j] /= m;
			F[i] /= m;

			if (i != k)
			{
				for (size_t j = 0; j < n; j++)
					A[i][j] -= A[k][j];
				F[i] -= F[k];
			}
		}
	}

	// Обратный ход
	for (int k = n - 1; k >= 0; k--)
	{
		deltaX[k] = F[k];
		for (size_t i = 0; i < k; i++)
			F[i] = F[i] - A[i][k] * deltaX[k];
	}
}

// Ищем x(k+1)
void snu::calc_xk1(tests test)
{
	beta = 1;
	prevX = x;
	for (size_t i = 0; i < x.size(); i++)
		x[i] += deltaX[i];

	f(test);
	currNormF = norm(F);
	for (size_t v = 0; beta > eps1 && currNormF > prevNormF; v++)
	{
		beta /= 2.0;
		x = prevX + beta * deltaX;
		f(test);
		currNormF = norm(F);
	}
}

// Метод Ньютона
void snu::newton(tests test)
{
	if (n > m)
		throw exceptions::CANT_SOLVE;

	std::ofstream coords("coords.txt");

	FILE* table;
	fopen_s(&table, "result.txt", "w");
	if (table)
	{
		fprintf(table, "| ITERS |     BETA     |       X       |       Y       |     NORM     |\n");
		fprintf(table, "|_____________________________________________________________________|\n");
	}

	f(test);
	normF = norm(F);
	currNormF = normF;

	if (table)
		fprintf(table, "|%-7llu|%-14f|%-15f|%-15f|%-14f|\n", size_t(0), beta, x[0], x[1], currNormF);

	coords << x[0] << " " << x[1] << std::endl;
	for (size_t i = 1; i < maxIter && currNormF / normF > eps2; i++)
	{
		prevNormF = norm(F);
		if (m > n)
			findMinF();
		createJacobiMatrix(test);

		gauss();

		calc_xk1(test);

		if (table)
			fprintf(table, "|%-7llu|%-14f|%-15f|%-15f|%-14f|\n", i, beta, x[0], x[1], prevNormF);
		coords << x[0] << " " << x[1] << std::endl;
	}
	if (table)
		fclose(table);

	coords.close();
}