#include "snu.h"
#include <algorithm>

// Конструктор
snu::snu(const std::string directory)
{
	n = m = maxIter = eps1 = eps2 = 0;
	beta = 1;
	currNormF = prevNormF = 0;
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

// Записать результат в файл
void snu::writeToFile(const std::string file)
{
	std::ofstream out(file);
	if (out.is_open())
	{
		for (size_t i = 0; i < x.size(); i++)
			out << x[i] << std::endl;
		out.close();
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
	// Иначе нужно исключить (m - n) уравнений с минимальными |F|
	else
		for (size_t l = 0, i = 0; l < m; l++)
			if (!isInMinF(l))
			{
				for (size_t j = 0; j < n; j++)
					A[i][j] = df(i, j);
				F[i++] = F[l];
			}
}

// Метод Гаусса
void snu::gauss()
{	
	deltaX.resize(n);
	real verySmall = 1e-10;

	for (size_t i = 0; i < n; i++)
	{
		// Ищем максимальный элемент в столбце
		size_t max = i;
		for (size_t j = i + 1; j < n; j++)
			if (fabs(A[j][i]) > fabs(A[max][i]))
				max = j;

		// Меняем текущую строку со строкой, в которой
		// находится максимальный элемент
		for (size_t j = i; j < n; j++)
			std::swap(A[i][j], A[max][j]);
		std::swap(F[i], F[max]);

		// Вычитаем из всех строк, расположенных ниже
		// текущей, строку с максимальным элементом,
		// умноженную на вычисленный коэффициент.
		// То же самое делаем для вектора
		for (size_t j = i + 1; j < n; j++)
		{
			if (abs(A[i][i]) < verySmall)
				throw snu::exceptions::VERY_SMALL_DIAG;

			real m = A[j][i] / A[i][i];
			F[j] -= m * F[i];
			for (size_t k = i; k < n; k++)
				A[j][k] -= m * A[i][k];
		}

		if (abs(A[i][i]) < verySmall)
			throw snu::exceptions::VERY_SMALL_DIAG;

		real m = A[i][i];
		for (size_t j = i; j < n; j++)
			A[i][j] /= m;
		F[i] /= m;
	}

	// Обратный ход
	for (int i = n - 1; i >= 0; i--)
	{
		real prod = F[i];
		for (int j = n - 1; j >= i; j--)
			prod -= A[i][j] * deltaX[j];

		if (abs(A[i][i]) < verySmall)
			throw snu::exceptions::VERY_SMALL_DIAG;

		prod /= A[i][i];
		deltaX[i] = prod;
	}
}

// Ищем x(k+1)
void snu::calc_xk1()
{
	beta = 1;
	prevX = x;
	for (size_t i = 0; i < x.size(); i++)
		x[i] += deltaX[i];

	f();
	currNormF = norm(F);
	for (size_t v = 0; beta > eps1 && currNormF > prevNormF; v++)
	{
		beta /= 2;
		x = prevX + beta * deltaX;
		f();
		currNormF = norm(F);
	}
}

// Метод Ньютона
void snu::newton()
{
	if (n > n)
		throw snu::exceptions::CANT_SOLVE;

	f();

	currNormF = norm(F);
	prevNormF = currNormF;

	for (size_t i = 0; i < maxIter && currNormF / prevNormF > eps2; i++)
	{
		if (m > n)
			findMinF();
		createJacobiMatrix();

		gauss();

		calc_xk1();
	}

	writeToFile("result.txt");
}