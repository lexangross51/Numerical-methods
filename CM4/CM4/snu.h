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
	enum class exceptions { VERY_SMALL_DIAG, CANT_SOLVE };

public:
	snu(const std::string directory);

	real f(size_t funcNumber);
	real df(size_t funcNumber, size_t varNumber);

	void f();
	void createJacobiMatrix();

	void gauss();
	void newton();

	void calc_xk1();

	void writeToFile(const std::string file);

private:
	size_t m;		// Кол-во уравнений
	size_t n;		// Кол-во неизвестных
	size_t maxIter;	// Максимальное число итераций
	real eps1;		// Точность для бетты
	real eps2;		// Точность для невязки

	Vector F;		// Вектор правой части
	real currNormF;	// Текущая норма F
	real prevNormF;	// Предыдущая норма F
	Vector x;		// Вектор неизвестных
	Vector prevX;	// Предыдущее значение x
	Vector deltaX;	// Направление
	real beta;		// Параметра бетта

	Matrix A;		// Матрица Якоби

	std::vector<size_t> minF;	// Массив с индексами, которые
								// соответсвуют минимальным F
	real norm(Vector &v);

	void findMinF();

	bool isInMinF(size_t equationNum);
};

inline Vector operator +(const Vector& a, const Vector& b)
{
	Vector res = a;
	for (size_t i = 0; i < res.size(); i++)
		res[i] += b[i];
	return res;
}

inline Vector operator *(const real& a, const Vector& b)
{
	Vector res = b;
	for (size_t i = 0; i < res.size(); i++)
		res[i] *= a;
	return res;
}

#endif

