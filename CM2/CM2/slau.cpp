#include "slau.h"
#include <fstream>
#include <iomanip>

// Считываем максимальное число итераций и точность решения
void slau::readParameters(std::string name)
{
    std::ifstream fin(name);
    if (fin.is_open())
    {
        fin >> maxIter >> eps >> blockSize;
    }
    fin.close();
}

// Считываем вектор правой части
void slau::readVector(std::string name)
{
    std::ifstream fin(name);
    if (fin.is_open())
    {
        F.resize(matrix::n);
        for (size_t i = 0; i < F.size(); i++)
            fin >> F[i];
    }
    fin.close();
}

// Считываем все, что необходимо для решения СЛАУ:
// матрицу, вектор, параметры для вычислений
void slau::readSLAU()
{
    readParameters("parameters.txt");
    matrix::read("matrix.txt");
    readVector("vector.txt");
}

// Задаем начальное приближение в виде нулевого вектора
void slau::setInitialApprox()
{
    xk.clear();
    xk.resize(matrix::n, 1);
    normF = norm(F);
}

// Записываем результат в файл
void slau::writeResult(std::string name)
{
    std::ofstream fout(name, std::ios::app);
    if (fout.is_open())
    {
        for (size_t i = 0; i < xk.size(); i++)
            fout << std::setprecision(15) << xk[i] << std::endl;
        fout << "______________" << std::endl;
    }
    fout.close();
}

// Метод Якоби с параметром релаксации
void slau::Jacobi(real w)
{ 
    xk1.clear();
    xk1.resize(n, 0);

    for (size_t i = 0; i < xk.size(); i++)
    {
        real sum = matrix::multLineOnVector(i, xk, 1);
        xk1[i] = xk[i] + w * (F[i] - sum) / di[i];
    }
    xk = xk1;
}

// Метод Гаусса-Зейделя с параметром релаксации
void slau::Gauss_Seidel(real w)
{
    xk1.clear();
    xk1.resize(n, 0);

    real suml = 0;
    real sumu = 0;

    for (size_t i = 0; i < xk.size(); i++)
    {
        suml = matrix::multLineOnVector(i, xk1, 2);
        sumu = matrix::multLineOnVector(i, xk, 3);
        xk1[i] = xk[i] + w * (F[i] - suml - sumu) / di[i];
    }
    xk = xk1;
}

// Решаем СЛАУ одним из итерационных методов:
// 0 - Якоби, 1 - Гаусс-Зейдель;
// выводим получившееся число итераций
void slau::solve_by_iterative(real w, int method)
{
    size_t i;
    for (i = 0; i < maxIter && relative_discrepancy() >= eps; i++)
    {
        if (method == 0)
            Jacobi(w);
        else
            Gauss_Seidel(w);
    }
    std::cout << "Number of iterations: " << i << std::endl;
    writeTable("result.txt", w, i);
}

// LU-разложение главных блоков
void slau::LU()
{
    size_t blocksCount = n / blockSize; // Кол-во блоков

    // Проходим по всем главным блокам
    for (size_t i = 0; i < blocksCount; i++)
    {
        size_t k = blockSize - 1 + i * blockSize;
        for (size_t j = i * blockSize; j < k; j++)
        {
            u1[j] /= di[j];
            di[j + 1] -= l1[j] * u1[j];
        }
    }
}

// Решение СЛАУ прямым и обратным ходом
void slau::LUSolution(size_t blockNumber, real w)
{
    size_t j = (blockNumber + 1) * blockSize;
    for (size_t i = blockNumber * blockSize; i < j; i++)
        xk[i] *= w;

    size_t k = blockNumber * blockSize;

    // Прямой ход
    xk[k] /= di[k];
    for (size_t i = k + 1; i < j; i++)
    {
        xk[i] -= xk[i - 1] * l1[i - 1];
        if (di[i] == 0)
        {
            std::cerr << "Cannot be devided by zero!" << std::endl;
            break;
        }
        xk[i] /= di[i];
    }

    // Обратный ход
    __int64 end = j - blockSize;
    for (__int64 i = j - 2; i >= end; i--)
        xk[i] -= u1[i] * xk[i + 1];
}

// Умножение блока на вектор
void slau::multBlockOnVector(size_t blockNumber)
{
    size_t j = blockSize * blockNumber;
    for (int i = j; i < blockSize + j; i++)
        xk[i] = 0;

    if (j > 0)
        xk[j] += l1[j - 1] * xk[j - 1];

    for (int l = 0; l < blockSize; l++, j++)
    {
        //умножение нижнего треугольника матрицы
        if (j > m + 1)
        {
            xk[j] += l2[j - m - 2] * xk[j - m - 2];
            if (j > m + 2)
                xk[j] += l3[j - m - 3] * xk[j - m - 3];
        }

        // умножение верхнего треугольника матрицы
        if (j < n - m - 2)
        {
            xk[j] += u2[j] * xk[j + m + 2];
            if (j < n - m - 3)
                xk[j] += u3[j] * xk[j + m + 3];
        }
    }

    // Домножаем оставшийся элемент из 1-й диагонали верхнего треугольника
    // (если он есть)
    j--;
    if (j < n - 1)
        xk[j] += u1[j] * xk[j + 1];
}
 
// Метод блочной релаксации
void slau::Block_Relaxation(real w, size_t blocknumber, real &error)
{
    // Выполняем цикл столько раз, сколько главных блоков
    for (size_t i = 0; i < blocknumber; i++)
    {
        xk1 = xk;
        // Умножаем блок на вектор
        multBlockOnVector(i);

        // Вычитаем из вектора правой части посчитанную сумму
        size_t end = (i + 1) * blockSize;
        for (size_t k = i * blockSize; k < end; k++)
            xk[k] = F[k] - xk[k];

        // Решаем СЛАУ (прямой и обратный ходы)
        LUSolution(i, w);

        for (size_t k = i * blockSize; k < end; k++)
            xk[k] += (1 - w) * xk1[k];
    }

    // Считаем погрешность (нужно для выходы из итерацинного процесса)
    error = norm(xk);
    for (size_t i = 0; i < n; i++)
        xk1[i] -= xk[i];
    error = norm(xk1) / error;
}

// Итерации для метода блочной релаксации
void slau::BR_Iterations(real w)
{
    if (n % blockSize != 0)
    {
        std::cerr << "Cannot be devided into equal blocks!" << std::endl;
        return;
    }

    // Считаем число главных блоков
    size_t blocknumber = n / blockSize;

    real error = 1; 
    size_t i;

    LU();
    for (i = 0; i < maxIter && error >= eps; i++)
        Block_Relaxation(w, blocknumber, error);

    writeTable("result.txt", w, i);
    std::cout << "Iters = " << i << std::endl;
    std::cout << "Error = " << error << std::endl;
}

// Считаем норму вектора
real slau::norm(std::vector<real> v)
{
    real sum = 0;
    for (size_t i = 0; i < v.size(); i++)
        sum += v[i] * v[i];
    return sqrt(sum);
}

// Считаем относительную невязку
real slau::relative_discrepancy()
{
    std::vector<real> v(matrix::n, 0);
    // v = A * xk
    matrix::mult(xk, v);
    // v = A * xk - F
    for (size_t i = 0; i < v.size(); i++)
        v[i] -= F[i];
    
    // || F - A * xk || / || F ||
    return norm(v) / normF;
}

// Выводим табличку с результатами
void slau::writeTable(std::string name, real w, size_t iter)
{
    std::ofstream fout(name, std::ios::app);
    if (fout.is_open())
    {
        std::vector<real> x(matrix::n);
        for (size_t i = 0; i < x.size(); i++)
            x[i] = i + 1;

        fout << "---------------------------------------------------------------" << std::endl;
        fout <<std::left << std::setw(5) << "w" << std::setw(20) << "|x" << std::setw(20) 
             << " ||x*-x|" << std::setw(10) << "  |iterations" <<  std::setw(4) << "|cond" << std::endl;
        fout << "---------------------------------------------------------------" << std::endl;
        fout << std::left << std::setprecision(14) << std::setw(5) << w << '|' 
             << std::setw(20) << xk[0] << '|' << std::setw(20) << fabs(xk[0] - x[0]) << '|' 
             << std::setw(10)  << iter << '|' << std::setw(3) << cond() << std::endl;
        for (size_t i = 1; i < x.size(); i++)
        {
            fout << std::left << std::setprecision(14) << std::setw(5) << "" << '|' 
                 << std::setw(20) << xk[i] << '|' << std::setw(20) << fabs(xk[i] - x[i]) 
                 << '|' << std::setw(10) << "" << '|' << std::setw(6) << "" << std::endl;
        }
        fout << "---------------------------------------------------------------" << std::endl << std::endl << std::endl;
    }
    fout.close();

}

// Расчет числа обусловленности матрицы
double slau::cond()
{
    std::vector<real> error(matrix::n, 0);

    for (size_t i = 0; i < error.size(); i++)
        error[i] = fabs(i + 1 - xk[i]);
    
    // ||x*|| = sqrt(385)
    return round((norm(error) / sqrt(385.0)) / relative_discrepancy());
}