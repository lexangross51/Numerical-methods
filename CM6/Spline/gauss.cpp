#include "gauss.h"

gauss::gauss()
{
    // Точки Гаусса
    gauss_points.resize(3);
    gauss_points[0] = 0.0;
    gauss_points[1] = sqrt(3.0 / 5);
    gauss_points[2] = -sqrt(3.0 / 5);

    // Квадратурные веса Гаусса
    weights.resize(3);
    weights[0] = 8.0 / 9;
    weights[1] = 5.0 / 9;
    weights[2] = 5.0 / 9;
}

// Решение СЛАУ методом Гаусса
gauss::complete gauss::solve_slae(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& q)
{
    double verySmall = 1e-16;

    for (size_t i = 0; i < A.size(); i++)
    {
        size_t max = i;
        for (size_t j = i + 1; j < A.size(); j++)
            if (fabs(A[j][i]) > fabs(A[max][i]))
                max = j;

        for (size_t j = i; j < A.size(); j++)
            std::swap(A[i][j], A[max][j]);
        std::swap(b[i], b[max]);

        for (size_t j = i + 1; j < A.size(); j++)
        {
            if (fabs(A[i][i]) < verySmall)
                return complete::FAIL;

            double m = A[j][i] / A[i][i];
            b[j] -= m * b[i];

            for (size_t k = i; k < A.size(); k++)
                A[j][k] -= m * A[i][k];
        }

        double m = A[i][i];
        if (fabs(m) < verySmall)
            return complete::FAIL;

        for (size_t j = i; j < A.size(); j++)
            A[i][j] /= m;
        b[i] /= m;
    }

    for (int i = A.size() - 1; i >= 0; i--)
    {
        double prod = b[i];
        for (int j = A.size() - 1; j >= i; j--)
            prod -= A[i][j] * q[j];

        if (fabs(A[i][i]) < verySmall)
            return complete::FAIL;

        prod /= A[i][i];
        q[i] = prod;
    }
    return complete::SUCCESS;
}

// Численное интегрирование 3-х узловым методом Гаусса
double gauss::integrate(ScalFunc f, double xk, double xk1)
{
    double sum = 0, qi, h, pi;

    for (uint8_t i = 0; i < 3; i++)
    {
        qi = weights[i];
        h = fabs(xk1 - xk);
        pi = (xk + xk1 + gauss_points[i] * h) / 2;
        sum += qi * f(pi);
    }

    return h * sum / 2.0;
}
