#include "matrix.h"
#include <iomanip>

matrix::matrix(const std::string fileMatrix)
{
    dim = 0;
    max_eig = min_eig = 0;

    std::ifstream mat(fileMatrix);
    if (mat.is_open())
    {
        mat >> dim;
        A.resize(dim);

        for (size_t i = 0; i < dim; i++)
        {
            A[i].resize(dim);
            for (size_t j = 0; j < dim; j++)
                mat >> A[i][j];
        }
        mat.close();
    }
}

// Скалярное произведение
double matrix::sp(const std::vector<double>& a, const std::vector<double>& b)
{
    double scalar = 0;
    for (size_t i = 0; i < a.size(); i++)
        scalar += a[i] * b[i];
    return scalar;
}

// Умножение матрицы на вектор
void matrix::mult(const std::vector<double>& a, std::vector<double>& b)
{
    b.clear();
    b.resize(dim);
    for (size_t i = 0; i < dim; i++)
        for (size_t j = 0; j < dim; j++)
            b[i] += A[i][j] * a[j];
}

// Найти максимальное собственное значение 
void matrix::eigenval_max()
{
    double eps = 1e-15;
    size_t maxIter = 1000, iters = 0;

    double prev_max_eig = max_eig;
    std::vector<double> xk, xk_1(dim, 1);

    for (iters = 0; iters < maxIter; iters++)
    {
        mult(xk_1, xk);

        // Каждые 5 итераций делаем нормировку вектора xk
        if (iters % 5 == 0 && iters != 0)
        {
            double normX = sqrt(sp(xk, xk));
            for (size_t i = 0; i < dim; i++)
                xk[i] /= normX;
        }

        prev_max_eig = max_eig;
        max_eig = sp(xk, xk_1) / sp(xk_1, xk_1);
        xk_1 = xk;

        if (abs((max_eig - prev_max_eig) / max_eig) < eps)
            break;
    }

    std::cout << "Собственное значение " << max_eig << 
                 " найдено за " << iters << " итераций" << std::endl;
}

// LU разложение матрицы
void matrix::LU()
{
    for (size_t i = 0; i < dim; i++)
    {
        for (size_t j = 0; j < dim; j++)
        {
            if (i < j)
            {
                double sum = 0;
                for (size_t k = 0; k < i; k++)
                    sum += A[i][k] * A[k][j];
                A[i][j] -= sum;
                A[i][j] /= A[i][i];
            }
            else
            {
                double sum = 0;
                for (size_t k = 0; k < j; k++)
                    sum += A[i][k] * A[k][j];
                A[i][j] -= sum;
            }
        }
    }
}

// Решить СЛАУ с помощью LU разложения
void matrix::solve_by_LU(const std::vector<double>& b, std::vector<double>& res)
{
    res.resize(dim, 0);

    // Прямой ход
    for (size_t i = 0; i < dim; i++)
    {
        double sum = 0;
        for (size_t j = 0; j < i; j++)
            sum += A[i][j] * res[j];

        res[i] = b[i] - sum;
        res[i] /= A[i][i];
    }

    // Обратный ход
    for (__int64 i = dim - 1; i >= 0; i--)
    {
        double sum = 0;
        for (size_t j = i + 1; j < dim; j++)
            sum += A[i][j] * res[j];
        res[i] -= sum;
    }
}

// Найти мининмальное собственное значение 
void matrix::eigenval_min()
{
    double eps = 1e-15;
    size_t maxIter = 1000, iters = 0;

    double prev_min_eig = min_eig;
    std::vector<double> xk, xk_1(dim, 1);

    LU();

    for (iters = 0; iters < maxIter; iters++)
    {
        solve_by_LU(xk_1, xk);

        // Каждые 5 итераций делаем нормировку вектора xk
        if (iters % 5 == 0 && iters != 0)
        {
            double normX = sqrt(sp(xk, xk));
            for (size_t i = 0; i < xk.size(); i++)
                xk[i] /= normX;
        }

        prev_min_eig = min_eig;
        min_eig = sp(xk, xk_1) / sp(xk_1, xk_1);
        xk_1 = xk;

        if (abs((min_eig - prev_min_eig) / min_eig) < eps)
            break;
    }

    min_eig = 1 / min_eig;
    std::cout << "Собственное значение " << min_eig <<
        " найдено за " << iters << " итераций" << std::endl;
}