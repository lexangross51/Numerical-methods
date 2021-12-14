#ifndef HEAD_H
#define HEAD_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <functional>
#include <iomanip>
#include <cmath>

typedef std::function<double(double)> ScalFunc;

struct SLAE
{
    size_t size;
    std::vector<std::vector<double>> A;
    std::vector<double> q;
    std::vector<double> b;

    SLAE(size_t _size) : size(_size)
    {
        A.resize(size);
        for (size_t i = 0; i < size; i++)
            A[i].resize(size);
        q.resize(size);
        b.resize(size);
    }

    ~SLAE() {};
};

struct point
{
    double x;
    double y;
};

#endif
