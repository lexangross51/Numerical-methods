#pragma once
#ifndef SPLINE_H
#define SPLINE_H

#include "head.h"
#include "gauss.h"

class spline
{
    friend class MainWindow;
    friend class SndWin;
    friend class TrdWin;
public:
    spline();

    void read_nodes(const char* fileN);
    void read_points(const char* fileP);

    double psi(size_t func, double x, double h);
    double d_psi(size_t func, double x, double h);
    double dd_psi(size_t func, double x, double h);

    double alpha();
    double beta();

    void build();

    void calc_coefs();

    double calc_spline_in_point(double x);

private:
    size_t n;					// Кол-во точек
    size_t m;					// Кол-во КЭ
    std::vector<point> points;	// Координаты точек
    std::vector<double> nodes;	// Границы КЭ
    std::vector<point> splineVals;
    double alpha_reg, beta_reg; // Регуляризирующие коэффициенты

    SLAE* slae;					// СЛАУ

    gauss* Gauss;				// Гаусс (решение СЛАУ, интегрирование)
};

#endif
