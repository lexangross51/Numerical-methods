#include "slau.h"
#include <fstream>
#include <iostream>
#include <iomanip>

// Конструктор: считываем данные из файлов и сразу 
// создаем объект типа slau
slau::slau(const std::string directoryName)
{
    n = maxIter = 0;
    eps = 0.0;
    normF = 0.0;

    // Читаем размер матрицы, максимальное число итераций, 
    // и необходимую точность решения
    std::ifstream kuslau(directoryName + "kuslau.txt");
    if (kuslau.is_open())
    {
        kuslau >> n >> maxIter >> eps;
        kuslau.close();
    }
    else
        std::cerr << "Can't open file!" << std::endl;

    // Читаем главную диагональ
    di.resize(n, 0);
    std::ifstream fdi(directoryName + "di.txt");
    if (fdi.is_open())
    {
        for (size_t i = 0; i < di.size(); i++)
            fdi >> di[i];
        fdi.close();
    }
    else
        std::cerr << "Can't open file!" << std::endl;

    // Читаем указатели на строки матрицы
    ig.resize(n + 1, 0);
    std::ifstream fig(directoryName + "ig.txt");
    if (fig.is_open())
    {
        for (size_t i = 0; i < ig.size(); i++)
        {
            fig >> ig[i];
            ig[i]--;
        }
        fig.close();
    }
    else
        std::cerr << "Can't open file!" << std::endl;

    // Читаем указатели на столбцы
    jg.resize(ig.back(), 0);
    std::ifstream fjg(directoryName + "jg.txt");
    if (fjg.is_open())
    {
        for (size_t i = 0; i < jg.size(); i++)
        {
            fjg >> jg[i];
            jg[i]--;
        }
        fjg.close();
    }
    else
        std::cerr << "Can't open file!" << std::endl;

    // Читаем нижний треугольник
    ggl.resize(ig.back(), 0);
    std::ifstream fggl(directoryName + "ggl.txt");
    if (fggl.is_open())
    {
        for (size_t i = 0; i < ggl.size(); i++)
            fggl >> ggl[i];
        fggl.close();
    }
    else
        std::cerr << "Can't open file!" << std::endl;

    // Читаем верхний треугольник
    ggu.resize(ig.back(), 0);
    std::ifstream fggu(directoryName + "ggu.txt");
    if (fggu.is_open())
    {
        for (size_t i = 0; i < ggu.size(); i++)
            fggu >> ggu[i];
        fggu.close();
    }
    else
        std::cerr << "Can't open file!" << std::endl;

    // Читаем вектор правой части
    F.resize(n, 0);
    std::ifstream fpr(directoryName + "pr.txt");
    if (fpr.is_open())
    {
        for (size_t i = 0; i < F.size(); i++)
            fpr >> F[i];
        normF = sqrt(F * F);
        fpr.close();
    }
    else
        std::cerr << "Can't open file!" << std::endl;
}

// Записать результат в файл
void slau::writeResult(const std::string file)
{
    std::ofstream out(file, std::ios::app);
    if (out.is_open())
    {
        for (size_t i = 0; i < x.size(); i++)
            out << std::setprecision(15) << x[i] << std::endl;
        out << "--------------" << std::endl;
        out.close();
    }
    else
        std::cerr << "Can't open file!" << std::endl;
}

// Задать начальное приближение
void slau::setInitialApprox()
{
    // x = (0, ..., 0)
    x.clear();
    x.resize(n, 0);
}

// Сгенерировать вектор (1, 2, ... , size)
std::vector<real> slau::generateVector(size_t size)
{
    std::vector<real> v(size);
    for (size_t i = 0; i < v.size(); i++)
        v[i] = i + 1;
    return v;
}

// Перевод матрицы в плотный формат
void slau::toDense(const std::string _dense)
{
    mat.resize(n);
    for (size_t i = 0; i < mat.size(); i++)
        mat[i].resize(n, 0);

    for (size_t i = 0; i < mat.size(); i++)
    {
        mat[i][i] = di[i];
        for (size_t j = ig[i]; j < ig[i + 1]; j++)
        {
            mat[i][jg[j]] = ggl[j];
            mat[jg[j]][i] = ggu[j];
        }
    }

    std::ofstream dense(_dense);
    if (dense.is_open())
    {
        for (size_t i = 0; i < mat.size(); i++)
        {
            for (size_t j = 0; j < mat[i].size(); j++)
                dense << std::setw(12) << mat[i][j];
            dense << std::endl << std::endl;
        }
    }
}

// Умножение матрицы на вектор
std::vector<real> slau::mult(const std::vector<real>& v)
{
    std::vector<real> res(v.size(), 0);
    for (size_t i = 0; i < v.size(); i++)
    {
        res[i] = di[i] * v[i];
        for (size_t j = ig[i]; j < ig[i + 1]; j++)
        {
            res[i] += ggl[j] * v[jg[j]];
            res[jg[j]] += ggu[j] * v[i];
        }
    }
    return res;
}

// Умножение транспонированной матрицы на вектор
std::vector<real> slau::multT(const std::vector<real>& v)
{
    std::vector<real> res(v.size(), 0);
    for (size_t i = 0; i < v.size(); i++)
    {
        res[i] = di[i] * v[i];
        for (size_t j = ig[i]; j < ig[i + 1]; j++)
        {
            res[i] += ggu[j] * v[jg[j]];
            res[jg[j]] += ggl[j] * v[i];
        }
    }
    return res;
}

//==================================================================
//================== ЛОКАЛЬНО-ОПТИМАЛЬНАЯ СХЕМА ====================
//==================================================================
// Локально-оптимальная схема
void slau::LOS(real* residual, size_t* itersCount)
{
    // Задаем начальное прилижениие
    setInitialApprox();
    r.resize(n);

    real alpha, beta;
    std::vector<real> Ar;

    // r0 = f - Ax0
    r = F - mult(x);

    // z0 = r0;
    z = r;

    //p0 = Az0
    p = mult(z);

    // Выполняем цикл, пока не достигнем максимального числа итераций
    // либо невязка не станет меньше заданного значения
    size_t k = 0;
    real rr = (r * r);
    for (; k < maxIter && rr >= eps; k++)
    {
        // alpha = ( p(k-1), r(k-1) ) / ( p(k-1), p(k-1) )
        real pp = (p * p);
        alpha = (p * r) / pp;

        // x(k) = x(k-1) + alpha * z(k-1)
        x = x + alpha * z;

        // (r(k), r(k)) = (r(k-1), r(k-1)) - alpha^2 * (p(k-1), p(k-1))
        rr = (r * r) - alpha * alpha * pp;

        // r(k) = r(k-1) - alpha * p(k-1)
        r = r - alpha * p;

        // beta = -( p(k-1), Ar(k) ) / ( p(k-1), p(k-1) )
        Ar = mult(r);
        beta = -(p * Ar) / pp;

        // z(k) = r(k) + beta * z(k-1)
        z = r + beta * z;

        // p(k) = Ar(k) + beta * p(k-1)
        p = Ar + beta * p;
    }

    // Запоминаем итоговую невязку
    (*residual) = rr;
    // и число итераций
    (*itersCount) = k;
}

//---------------------------------------------------------------------------
// LU-факторизация матрицы
void slau::LU()
{
    dif.clear();
    gglf.clear();
    gguf.clear();

    dif = di;
    gglf = ggl;
    gguf = ggu;

    // Проходим по всем строкам
    for (size_t i = 0; i < n; i++)
    {
        real sumdi = 0.0;		// Сумма для изменения диагонального элемента

        size_t i0 = ig[i];		// индекс 1-ого ненулевого элемента в строке
        size_t i1 = ig[i + 1];	// индекс последнего ненулевого элемента в строке

        // Идем по профилю
        for (size_t k = i0; k < i1; k++)
        {
            size_t j = jg[k];		// текущий столбец
            size_t j0 = ig[j];		// строка, с которой начинаются 
                                    // элементы профиля столбца
            size_t j1 = ig[j + 1];	// строка, которой заканчиваются 
                                    // элементы профиля столбца

            size_t ik = i0;			// ik-й элемент нижнего треугольника
            size_t kj = j0;			// kj-й элемент верхнего треугольника

            real suml = 0.0;		// Сумма для нижнего треугольника
            real sumu = 0.0;		// Сумма для верхнего треугольника

            while (ik < k && kj < j1)
            {
                // Если элементы на перемножаемых позициях
                if ( jg[ik] == jg[kj])
                {
                    // Накапливаем их произведения в суммы
                    suml += gglf[ik] * gguf[kj];
                    sumu += gguf[ik] * gglf[kj];
                    ik++;
                    kj++;
                }
                // Иначе сдвигаем позиции
                else
                    jg[ik] > jg[kj] ? kj++ : ik++;
            }

            // И пересчитываем элементы
            gglf[k] -= suml;
            gguf[k] = (gguf[k] - sumu) / dif[j];
            sumdi += gglf[k] * gguf[k];
        }

        // Пересчитываем диагональный элемент			
        dif[i] -= sumdi;
    }
}

// Решение СЛАУ путем LU-разложения: прямой ход
std::vector<real> slau::LUDirect(const std::vector<real>& b)
{
    std::vector<real> y = b;

    for (size_t i = 0; i < y.size(); i++)
    {
        real sum = 0.0;
        for (size_t j = ig[i]; j < ig[i + 1]; j++)
            sum += gglf[j] * y[jg[j]];

        y[i] -= sum;
        y[i] /= dif[i];
    }

    return y;
}

// Решение СЛАУ путем LU-разложения: обратный ход
std::vector<real> slau::LUReverse(const std::vector<real>& y)
{
    std::vector<real> res = y;

    for (__int64 i = n - 1; i >= 0; i--)
    {
        for (size_t j = ig[i]; j < ig[i + 1]; j++)
            res[jg[j]] -= gguf[j] * res[i];
    }

    return res;
}

// Локально-оптимальная схема + LU-факторизация
void slau::LOSLU(real* residual, size_t* itersCount)
{
    // Задаем начальное приближение
    setInitialApprox();
    r.resize(n);

    real alpha, beta;
    std::vector<real> LAU, U;

    // Выполняем LU-факторизацию
    LU();

    // r0 = L^-1 * (f - Ax0)
    r = LUDirect(F - mult(x));

    // z0 = U^-1 * r0
    z = LUReverse(r);

    // p0 = L^-1 * Az0
    p = LUDirect(mult(z));

    // Выполняем цикл, пока не достигнем максимального числа итераций
    // либо невязка не станет меньше заданного значения
    size_t k = 0;
    real rr = (r * r);
    for ( ; k < maxIter && rr >= eps; k++)
    {
        // alpha = ( p(k-1), r(k-1) ) / ( p(k-1), p(k-1) )
        real pp = p * p;
        alpha = (p * r) / pp;

        // x(k) = x(k-1) + alpha * z(k-1)
        x = x + alpha * z;

        // (r(k), r(k)) = (r(k-1), r(k-1)) - alpha^2 * (p(k-1), p(k-1))
        rr = (r * r) - alpha * alpha * pp;

        // r(k) = r(k-1) - alpha * p(k-1)
        r = r - alpha * p;

        // L^-1 * A * U^-1 ~ LAU
        // beta = -( p(k-1), LAU * r(k) ) / ( p(k-1), p(k-1) )
        LAU = LUReverse(r);
        LAU = mult(LAU);
        LAU = LUDirect(LAU);
        beta = -(p * LAU) / pp;

        // z(k) = U^-1 * r(k) + beta * z(k-1)
        U = LUReverse(r);
        z = U + beta * z;

        // p(k) = LAU * r(k) + beta * p(k-1)
        p = LAU + beta * p;
    }

    // Запоминаем итоговую невязку
    (*residual) = rr;
    // и число итераций
    (*itersCount) = k;
}

//---------------------------------------------------------------------------
// Диагональное предобуславливание
void slau::diagFact()
{
    dif.clear();
    dif = di;
    for (size_t i = 0; i < n; i++)
        dif[i] = 1.0 / sqrt(dif[i]);
}

// Умножение диагональной матрицы на вектор
std::vector<real> slau::multD(const std::vector<real>& v)
{
    std::vector<real> res(n, 0);
    for (size_t i = 0; i < n; i++)
        res[i] = dif[i] * v[i];

    return res;
}
 
// Локально-оптимальная схема с диагональным предобуславливанием
void slau::LOSdiag(real* residual, size_t* itersCount)
{
    // Задаем начальное приближение
    setInitialApprox();
    r.resize(n);

    real alpha, beta;
    std::vector<real> LAU, U;

    // Выполняем диагональную факторизацию
    diagFact();

    // r0 = L^-1 * (f - Ax0)
    r = multD(F - mult(x));

    // z0 = U^-1 * r0
    z = multD(r);

    // p0 = L^-1 * Az0
    p = multD(mult(z));

    // Выполняем цикл, пока не достигнем максимального числа итераций
    // либо невязка не станет меньше заданного значения
    size_t k = 0;
    real rr = (r * r);
    for (; k < maxIter && rr >= eps; k++)
    {
        // alpha = ( p(k-1), r(k-1) ) / ( p(k-1), p(k-1) )
        real pp = (p * p);
        real alpha = (p * r) / pp;

        // x(k) = x(k-1) + alpha * z(k-1)
        x = x + alpha * z;

        // (r(k), r(k)) = (r(k-1), r(k-1)) - alpha^2 * (p(k-1), p(k-1))
        rr = (r * r) - alpha * alpha * pp;

        // r(k) = r(k-1) - alpha * p(k-1)
        r = r - alpha * p;

        // L^-1 * A * U^-1 ~ LAU
        // beta = -( p(k-1), LAU * r(k) ) / ( p(k-1), p(k-1) )
        LAU = multD(r);
        LAU = mult(LAU);
        LAU = multD(LAU);
        real beta = -(p * LAU) / pp;

        // z(k) = U^-1 * r(k) + beta * z(k-1)
        U = multD(r);
        z = U + beta * z;

        // p(k) = LAU * r(k) + beta * p(k-1)
        p = LAU + beta * p;
    }

    // Запоминаем итоговую невязку
    (*residual) = rr;
    // и число итераций
    (*itersCount) = k;
}

//---------------------------------------------------------------------------
// Переввести матрицу в разреженный формат
void slau::toSparse()
{
    di.resize(mat.size(), 0);
    ig.resize(mat.size() + 1, 0);

    // Записываем диагональ
    for (size_t i = 0; i < di.size(); i++)
        di[i] = mat[i][i];

    // Считаем размеры профилей строк
    ig[0] = ig[1] = 0;
    for (size_t i = 1; i < mat.size(); i++)
    {
        size_t count = 0;
        for (size_t j = 0; j < i; j++)
        {
            if (mat[i][j] != 0.0 && mat[j][i] != 0.0 ||
                mat[i][j] == 0.0 && mat[j][i] != 0.0 ||
                mat[i][j] != 0.0 && mat[j][i] == 0.0)
                count++;
        }
        ig[i + 1] = ig[i] + count;
    }

    // Определяем столбцы, в которых есть 
    // ненулевые элементы. (Идем построчно
    // по нижнему треугольнику)
    jg.resize(ig.back(), 0);
    size_t k = 0;
    for (size_t i = 0; i < mat.size(); i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            if (mat[i][j] != 0.0 && mat[j][i] != 0.0 ||
                mat[i][j] == 0.0 && mat[j][i] != 0.0 ||
                mat[i][j] != 0.0 && mat[j][i] == 0.0)
            {
                jg[k] = j;
                k++;
            }
        }
    }

    // Заполняем нижний и верхний треугольники
    ggl.resize(ig.back(), 0);
    ggu.resize(ig.back(), 0);
    k = 0;
    for (size_t i = 0; i < mat.size(); i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            if (mat[i][j] != 0.0 && mat[j][i] != 0.0 ||
                mat[i][j] == 0.0 && mat[j][i] != 0.0 ||
                mat[i][j] != 0.0 && mat[j][i] == 0.0)
            {
                ggl[k] = mat[i][j];
                ggu[k] = mat[j][i];
                k++;
            }
        }
    }
}

// Сгенерировать матрицу Гильберта и вектор заданного размера, 
// перевести матрицу в разреженный формат и 
// получить вектор правой части
void slau::createHilbertMatrix(size_t size)
{
    mat.clear();
    mat.resize(size);
    for (size_t i = 0; i < mat.size(); i++)
        mat[i].resize(size);

    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            mat[i][j] = 1.0 / real(i + j + 1);

    toSparse();

    F = mult(generateVector(size));
    normF = sqrt(F * F);
}

//==================================================================
//================= МЕТОД СОПРЯЖЕННЫХ ГРАДИЕНТОВ ===================
//==================================================================

// Решение СЛАУ путем LU-разложения: прямой ход для транспонированной матрицы
std::vector<real> slau::LUDirectT(const std::vector<real>& b)
{
    std::vector<real> y = b;

    for (size_t i = 0; i < y.size(); i++)
    {
        real sum = 0.0;
        for (size_t j = ig[i]; j < ig[i + 1]; j++)
            sum += gguf[j] * y[jg[j]];

        y[i] -= sum;
    }

    return y;
}

// Решение СЛАУ путем LU-разложения: обратный ход для транспонированной матрицы
std::vector<real> slau::LUReverseT(const std::vector<real>& y)
{
    std::vector<real> res = y;

    for (__int64 i = n - 1; i >= 0; i--)
    {
        res[i] /= dif[i];
        real sum = 0.0;
        for (size_t j = ig[i]; j < ig[i + 1]; j++)
            res[jg[j]] -= gglf[j] * res[i];
    }

    return res;
}

// Метод сопряженных градиентов + LU-факторизация
void slau::MSGLU(real* residual, size_t* itersCount)
{
    setInitialApprox();
    r.resize(n);

    std::vector<real> tmp;
    real alpha, beta;

    // Выполняем LU-факторизацию
    LU();

    // r = U^-T * A^T * L^-T * L^-1 * (F - Ax)
    r = F - mult(x);
    r = LUDirect(r);
    r = LUReverseT(r);
    r = multT(r);
    r = LUDirectT(r);

    z = r;

    real rr = (r * r);
    size_t i = 0;
    for ( ; i < maxIter && rr / normF >= eps; i++)
    {
        // tmp = U^-T * A^T * L^-T * L^-1 * A * U^-1 * z(k-1)
        tmp = LUReverse(z);
        tmp = mult(tmp);
        tmp = LUDirect(tmp);
        tmp = LUReverseT(tmp);
        tmp = multT(tmp);
        tmp = LUDirectT(tmp);

        // alpha = ( r(k-1), r(k-1) ) / ( tmp, z(k-1) ) 
        alpha = (r * r) / (tmp * z);

        // x(k) = x(k-1) + alpha * z(k-1)
        x = x + alpha * z;

        // r(k) = r(k-1) - alpha * tmp
        r = r - alpha * tmp;

        // beta = ( r(k), (k) ) / ( r(k-1), r(k-1) )
        beta = 1.0 / rr;
        rr = r * r;
        beta *= rr;

        // z(k) = r(k) + beta * z(k-1)
        z = r + beta * z;

    }

    x = LUReverse(x);

    // Запоминаем итоговую невязку
    (*residual) = rr / normF;
    // и число итераций
    (*itersCount) = i;
}

void slau::MSGdiag(real* residual, size_t* itersCount)
{
    setInitialApprox();
    r.resize(n);

    std::vector<real> tmp;
    real alpha, beta;

    // Выполняем диагональную факторизацию
    diagFact();

    // r = U^-T * A^T * L^-T * L^-1 * (F - Ax)
    r = F - mult(x);
    r = multD(r);
    r = multD(r);
    r = multT(r);
    r = multD(r);

    z = r;

    real rr = (r * r);
    size_t i = 0;
    for (; i < maxIter && rr / normF >= eps; i++)
    {
        // tmp = U^-T * A^T * L^-T * L^-1 * A * U^-1 * z(k-1)
        tmp = multD(z);
        tmp = mult(tmp);
        tmp = multD(tmp);
        tmp = multD(tmp);
        tmp = multT(tmp);
        tmp = multD(tmp);

        // alpha = ( r(k-1), r(k-1) ) / ( tmp, z(k-1) ) 
        alpha = (r * r) / (tmp * z);

        // x(k) = x(k-1) + alpha * z(k-1)
        x = x + alpha * z;

        // r(k) = r(k-1) - alpha * tmp
        r = r - alpha * tmp;

        // beta = ( r(k), (k) ) / ( r(k-1), r(k-1) )
        beta = 1.0 / rr;
        rr = r * r;
        beta *= rr;

        // z(k) = r(k) + beta * z(k-1)
        z = r + beta * z;

    }

    x = multD(x);

    // Запоминаем итоговую невязку
    (*residual) = rr / normF;
    // и число итераций
    (*itersCount) = i;
}