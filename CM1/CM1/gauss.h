#pragma once

#include <cmath>

void slau_by_gauss(real* a, real* b, real*& x, size_t n)
{
    x = new real[n];
    for (size_t i = 0; i < n; i++)
        x[i] = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        // Ищем максимальный элемент в столбце
        size_t max = i;
        for (size_t j = i + 1; j < n; j++)
            if (fabs(a[j * n + i]) > fabs(a[max * n + i]))
                max = j;

        // Меняем текущую строку со строкой, в которой
        // находится максимальный элемент
        for (size_t j = i; j < n; j++)
            std::swap(a[i * n + j], a[max * n + j]);
        std::swap(b[i], b[max]);

        // Вычитаем из всех строк, расположенных ниже
        // текущей, строку с максимальным элементом,
        // умноженную на вычисленный коэффициент.
        // То же самое делаем для вектора
        for (size_t j = i + 1; j < n; j++)
        {
            real m = a[j * n + i] / a[i * n + i];
            b[j] -= m * b[i];
            for (size_t k = i; k < n; k++)
                a[j * n + k] -= m * a[i * n + k];
        }
        real m = a[i * n + i];
        for (size_t j = i; j < n; j++)
            a[i * n + j] /= m;
        b[i] /= m;
    }

    // Обратный ход
    for (int i = n - 1; i >= 0; i--)
    {
        real prod = b[i];
        for (int j = n - 1; j >= i; j--)
            prod -= a[i * n + j] * x[j];
        prod /= a[i * n + i];
        x[i] = prod;
    }
}