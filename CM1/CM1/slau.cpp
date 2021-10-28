#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "lu.h"
#include "ax_b.h"
#include "hilbert.h"
#include "multiplication.h"
#include "toProfile.h"
#include "gauss.h"

typedef double real;

// Функция считывания данных из файлов
void read(size_t& size, real* matrix[4], real*& vector)
{
    std::ifstream file_size("size.txt");
    if (file_size.is_open())
        // Читаем размер матрицы и вектора правой части
        file_size >> size;
    file_size.close();

    std::ifstream file_matrix("matrix.txt");
    if (file_matrix.is_open())
    {
        std::string line;
        std::stringstream ss;

        // Выделяем память под массив ia и читаем его
        matrix[0] = new real[size + 1];
        getline(file_matrix, line);
        ss << line;
        for (size_t i = 0; i < size + 1; i++)
            ss >> matrix[0][i];

        // Выделяем память под массив di и читаем его
        matrix[1] = new real[size];
        getline(file_matrix, line);
        ss << line;
        for (size_t i = 0; i < size; i++)
            ss >> matrix[1][i];

        // Размер массивов al и au равен ia[last] - 1
        size_t n = matrix[0][size] - 1;
        // Выделяем память под массив al и читаем его
        matrix[2] = new real[n];
        getline(file_matrix, line);
        ss << line;
        for (size_t i = 0; i < n; i++)
            ss >> matrix[2][i];

        // Выделяем память под массив au и читаем его
        matrix[3] = new real[n];
        getline(file_matrix, line);
        ss << line;
        for (size_t i = 0; i < n; i++)
            ss >> matrix[3][i];
    }
    file_matrix.close();

    std::ifstream file_vector("vector.txt");
    if (file_vector.is_open())
    {
        std::string line;
        std::stringstream ss;
        // Выделяем память под вектор и читаем его
        vector = new real[size];
        getline(file_vector, line);
        ss << line;
        for (size_t i = 0; i < size; i++)
            ss >> vector[i];
    }
    file_vector.close();
}

// Записываем результат в файл
void write(real* x, size_t n)
{
    std::ofstream file_res;
    file_res.open("x.txt", std::ios::app);
    if (file_res.is_open())
        for (size_t i = 0; i < n; i++)
            file_res << x[i] << std::endl;
    file_res.close();
}

int main()
{
    real* a[4];           // Матрица: a[0] - массив ia
                          //          a[1] - массив di  
                          //          a[2] - массив al      
                          //          a[3] - массив au

    real* x;              // Вектор неизвестных
    real* b;              // Вектор правой части
    //size_t n;			  // Размер

    size_t n = 10;			  // Размер

    real* m;
    //for (size_t k = 1; k <= 13; k++)
    //{
        //real* hm = hilbert_matrix(k);
    readd(m, n);
    toProfile(m, a, n);
    multiply(a[0], a[1], a[2], a[3], b, n);

    //slau_by_gauss(hm, b, x, k);

   // read(n, a, b);
    LU(n, a[0], a[1], a[2], a[3]);
    x_solve(a[0], a[3], y(a[0], a[1], a[2], b, n), x, n);
    write(x, n);
    //}

    for (size_t i = 0; i < 4; i++)
        delete[] a[i];
    delete[] x;
    delete[] b;

    return 0;
}