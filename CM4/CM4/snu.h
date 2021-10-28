#ifndef SNU_H
#define SNU_H

#include <vector>
#include <string>
#include <fstream>

typedef double real;
typedef std::vector<real> Vector;
typedef std::vector<std::vector<real>> Matrix;

class snu
{

public:
	snu(const std::string directory);

	real f(size_t funcNumber);
	void f();
	real df(size_t funcNumber, size_t varNumber);

	void createJacobiMatrix();

private:
	size_t m;		// Кол-во уравнений
	size_t n;		// Кол-во неизвестных
	size_t maxIter;	// Максимальное число итераций
	real eps1;		// Точность для бетты
	real eps2;		// Точность для невязки

	Vector F;		// Вектор правой части
	Vector x;		// Вектор неизвестных
	real beta;		// Параметра бетта

	Matrix A;		// Матрица Якоби

	std::vector<size_t> minF;	// Массив с индексами, которые
								// соответсвуют минимальным F
	real norm(Vector &v);

	void findMinF();

	bool isInMinF(size_t equationNum);
};

#endif

