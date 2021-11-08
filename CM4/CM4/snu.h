#ifndef SNU_H
#define SNU_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

typedef double real;
typedef std::vector<real> Vector;
typedef std::vector<std::vector<real>> Matrix;

class snu
{
public:
    enum class exceptions { VERY_SMALL_DIAG, CANT_SOLVE };

    // Тесты:
    // TEST1 - две непересекающиеся окружности
    // TEST2 - две окружности, пересекающиеся в одной точке
    // TEST3 - две окружности, пересекающиеся в двух точках
    // TEST4 - к окружностям добавлена прямая
    // TEST5 - три попарно пересекающиеся прямые
    // TEST6 - синусоида и прямая
    enum class tests { TEST1, TEST2, TEST3, TEST4, TEST5, TEST6 };

public:
    snu(const std::string directory);

    real f(tests test, size_t funcNumber);
    real df(tests test, size_t funcNumber, size_t varNumber);

    void f(tests test);
    void createJacobiMatrix(tests test);

    void gauss();
    void newton(tests test);

    void calc_xk1(tests test);

private:
    size_t m;		// Кол-во уравнений
    size_t n;		// Кол-во неизвестных
    size_t maxIter;	// Максимальное число итераций
    real eps1;		// Точность для бетты
    real eps2;		// Точность для невязки

    Vector F;		// Вектор правой части
    real normF;		// Норма вектора F до начала итераций 
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
};

inline Vector operator+(const Vector& a, const Vector& b)
{
    Vector res = a;
    for (size_t i = 0; i < res.size(); i++)
        res[i] += b[i];
    return res;
}

inline Vector operator*(const real& a, const Vector& b)
{
    Vector res = b;
    for (size_t i = 0; i < res.size(); i++)
        res[i] *= a;
    return res;
}

#endif

