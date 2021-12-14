#include "spline.h"

using namespace std;

spline::spline()
{
    m = n = 0;

    alpha_reg = beta_reg = 0.0;

    Gauss = new gauss();
}

// Считать узлы
void spline::read_nodes(const char* fileN)
{
    std::ifstream in(fileN);
    if (in.is_open())
    {
        in >> m;
        nodes.resize(m);
        for (size_t i = 0; i < m; i++)
            in >> nodes[i];
        in.close();

        slae = new SLAE(2 * m);
    }
    else
        throw ("File can't be opened!");
}

// Считать точки
void spline::read_points(const char* fileP)
{
    std::ifstream in(fileP);
    if (in.is_open())
    {
        in >> n;
        points.resize(n);
        for (size_t i = 0; i < n; i++)
            in >> points[i].x >> points[i].y;
        in.close();
    }
    else
        throw ("File can't be opened!");
}

// Значение базисной функции в точке
double spline::psi(size_t func, double x, double h)
{
    switch (func)
    {
    case 0:
        return 1 - 3*x*x + 2*x*x*x;

    case 1:
        return h * (x - 2*x*x + x*x*x);

    case 2:
        return 3*x*x - 2*x*x*x;

    case 3:
        return h * (-x*x + x*x*x);
    }
}

// Первая производная от базисной функции
double spline::d_psi(size_t func, double x, double h)
{
    switch (func)
    {
    case 0:
        return -6 * (x - x*x) / h;

    case 1:
        return 1 - 4*x + 3*x*x;

    case 2:
        return 6 * (x - x*x) / h;

    case 3:
        return -2*x + 3*x*x;
    }
}

// Вторая производная от базисной функции
double spline::dd_psi(size_t func, double x, double h)
{
    switch (func)
    {
    case 0:
        return -6 * (1 - 2 * x) / (h * h);

    case 1:
        return (-4 + 6 * x) / h;

    case 2:
        return 6 * (1 - 2 * x) / (h * h);

    case 3:
        return (-2 + 6 * x) / h;
    }
}

// Регуляризирующий коэффициент альфа
double spline::alpha()
{
    return alpha_reg;
}

// Регуляризирующий коэффициент бета
double spline::beta()
{
    return beta_reg;
}

// Плсчитать коэффициенты разложения
void spline::calc_coefs()
{
    size_t currPos = 0;

    // Проходим по КЭ
    for (size_t k = 0; k < m - 1; k++)
    {
        // Определяем его длину
        double h_i = nodes[k + 1] - nodes[k];

        // Если еще не прошли по всем точкам
        if (currPos < n && points[currPos].x >= nodes[k])
            // Проходим по всем точкам КЭ, начиная с currPos
            for (; currPos < n && points[currPos].x < nodes[k + 1]; currPos++)
            {
                double ksi = (points[currPos].x - nodes[k]) / h_i;

                // Идем по всем базисным функциям
                for (size_t psi_i = 0; psi_i < 4; psi_i++)
                {
                    double value = psi(psi_i, ksi, h_i);

                    // Заполняем матрицу СЛАУ
                    for (size_t psi_j = 0; psi_j < 4; psi_j++)
                        slae->A[2 * k + psi_i][2 * k + psi_j] += value * psi(psi_j, ksi, h_i);

                    // И вектор правой части
                    slae->b[2 * k + psi_i] += value * points[currPos].y;
                }
            }

        // Учтем регуляризирующие коэффициента альфа и бета
        for (size_t psi_i = 0; psi_i < 4; psi_i++)
            for (size_t psi_j = 0; psi_j < 4; psi_j++)
            {
                ScalFunc f = [this, psi_i, psi_j, h_i](double x)
                {
                    double d_psi_i = d_psi(psi_i, x, h_i);
                    double d_psi_j = d_psi(psi_j, x, h_i);

                    double dd_psi_i = dd_psi(psi_i, x, h_i);
                    double dd_psi_j = dd_psi(psi_j, x, h_i);

                    return (alpha() * d_psi_i * d_psi_j + beta() * dd_psi_i * dd_psi_j);
                };
                slae->A[2 * k + psi_i][2 * k + psi_j] += Gauss->integrate(f, nodes[k], nodes[k + 1]);
            }
    }

    Gauss->solve_slae(slae->A, slae->b, slae->q);
}

// Выдать значение сплайна в точке
double spline::calc_spline_in_point(double x)
{
    if (x < nodes[0] || x > nodes[m - 1])
        throw ("out of range");

    size_t elem = 0;
    bool find = false;
    for (size_t i = 0; i < m - 1 && !find; i++)
        if (x >= nodes[i] && x <= nodes[i + 1])
        {
            elem = i;
            find = true;
        }

    double h = nodes[elem + 1] - nodes[elem];
    double ksi = (x - nodes[elem]) / h;
    double value = 0;

    for (size_t psi_i = 0; psi_i < 4; psi_i++)
    {
        value += slae->q[2 * elem + psi_i] * psi(psi_i, ksi, h);
    }

    return value;
}

// Построить сплайн
void spline::build()
{
    calc_coefs();

    for (double x = nodes[0]; x <= nodes[m - 1]; x += 0.1)
        splineVals.push_back({ x, calc_spline_in_point(x) });
     splineVals.push_back({ points[n - 1].x, calc_spline_in_point(points[n - 1].x) });

    std::ofstream coords("coords.txt.");
    if (coords.is_open())
    {
        for (size_t i = 0; i < splineVals.size(); i++)
            coords << std::setw(4) << splineVals[i].x << std::setw(20) << splineVals[i].y << std::endl;
        coords.close();
    }
}
