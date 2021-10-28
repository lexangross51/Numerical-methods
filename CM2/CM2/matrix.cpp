#include "matrix.h"
#include <fstream>

// Читаем матрицу из файла
void matrix::read(std::string name)
{
    std::ifstream fin(name);
    if (fin.is_open())
    {		
        fin >> n >> m;

        di.resize(n);
        for (size_t i = 0; i < di.size(); i++)
            fin >> di[i];

        l1.resize(n - 1);
        for (size_t i = 0; i < l1.size(); i++)
            fin >> l1[i];

        l2.resize(n - m - 2);
        for (size_t i = 0; i < l2.size(); i++)
            fin >> l2[i];

        l3.resize(n - m - 3);
        for (size_t i = 0; i < l3.size(); i++)
            fin >> l3[i];

        u1.resize(n - 1);
        for (size_t i = 0; i < u1.size(); i++)
            fin >> u1[i];

        u2.resize(n - m - 2);
        for (size_t i = 0; i < u2.size(); i++)
            fin >> u2[i];

        u3.resize(n - m - 3);
        for (size_t i = 0; i < u3.size(); i++)
            fin >> u3[i];

        fin.close();
    }
}

// Умножение i-й строки матрицы на вектор v. mode указывает, как нужно умножать 
// строку на вектор: 1-всю строку на вектор (для Якоби), 2 - от начала строки 
// до диагонали(не включая), 3 - от диагонали до конца строки (для Гаусса-Зейделя)
real matrix::multLineOnVector(size_t i, const std::vector<real> &v, int mode)
{
    real sum = 0.0;

    if (mode == 1 || mode == 2)
    {
        // Если не в 1-й строке,
        if (i > 0)
        {
            // то точно попадется элемент 1-й диагонали 
            // из нижнего треугольника
            sum += l1[i - 1] * v[i - 1];

            // если захватываем и 2-ю нижнюю диагональ
            if (i > m + 1)
                sum += l2[i - m - 2] * v[i - m - 2];
            // если захватываем и 3-ю нижнюю диагональ
            if (i > m + 2)
                sum += l3[i - m - 3] * v[i - m - 3];
        }
    }
    if (mode == 1 || mode == 3)
    {
        // Умножаем диагональный на соответствующий 
        // элемент вектора
        sum += di[i] * v[i];
        
        // Если не в последней строке,
        if (i < n - 1)
        {
            // то точно попадется элемент 1-й диагонали
            // верхнего треугольника
            sum += u1[i] * v[i + 1];
            // если захватываем и 2-ю верхнюю диагональ
            if (i < n - m - 2)
                sum += u2[i] * v[i + m + 2];
            // если захватываем и 3-ю верхнюю диагональ
            if (i < n - m - 3)
                sum += u3[i] * v[i + m + 3];
        }
    }
    return sum; 
}

// Умножение матрицы на вектор v. 
void matrix::mult(const std::vector<real> &v, std::vector<real>& res)
{
    res.resize(n);
    for (size_t i = 0; i < n; i++)
        res[i] = multLineOnVector(i, v, 1);
}

#pragma region Перевод в плотный формат
void matrix::toDense(std::vector<std::vector<real>>& mat)
{
    mat.clear();
    mat.resize(n);

    for (size_t i = 0; i < n; i++)
    {
        mat[i].resize(n, 0);
        mat[i][i] = di[i];
    }

    size_t j = 1;
    for (size_t i = 0; i < l1.size(); i++, j++)
    {
        mat[i][j] = u1[i];
        mat[j][i] = l1[i];
    }

    j = m + 2;
    for (size_t i = 0; i < l2.size(); i++, j++)
    {
        mat[i][j] = u2[i];
        mat[j][i] = l2[i];
    }

    j = m + 3;
    for (size_t i = 0; i < l3.size(); i++, j++)
    {
        mat[i][j] = u3[i];
        mat[j][i] = l3[i];
    }
}

#pragma endregion