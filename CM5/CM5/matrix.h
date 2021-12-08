#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

class matrix
{
public:
    matrix(const std::string fileMatrix);

    void eigenval_max();

    void eigenval_min();

private:
    size_t dim;
    std::vector<std::vector<double>> A;

    double max_eig, min_eig;

    double sp(const std::vector<double>& a, const std::vector<double>& b);
    
    void mult(const std::vector<double>& a, std::vector<double>& b);

    void LU();
    void solve_by_LU(const std::vector<double>& b, std::vector<double>& res);
};

#endif