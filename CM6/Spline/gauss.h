#pragma once
#ifndef GAUSS_H
#define GAUSS_H

#include "head.h"

class gauss
{
    enum class complete { SUCCESS, FAIL };

public:
    gauss();

    complete solve_slae(std::vector<std::vector<double>>& A,
        std::vector<double>& b,
        std::vector<double>& q);

    double integrate(ScalFunc f, double xk, double xk1);

private:
    std::vector<double> gauss_points;
    std::vector<double> weights;
};

#endif
