#pragma once
#include "matrix.h"

typedef double real;

class slau : public matrix
{
public:
    slau() { maxIter = eps = blockSize = normF = 0; }

    void readSLAU();
    void solve_by_iterative(real w, int method);
    void setInitialApprox();
    void Jacobi(real w);
    void Gauss_Seidel(real w);

    void BR_Iterations(real w);
    void Block_Relaxation(real w, size_t blocknumber, real& error);
    void LU();
    void LUSolution(size_t blockNumber, real w);
    void multBlockOnVector(size_t blockNumber);


    void writeResult(std::string);
    void writeTable(std::string, real, size_t);
    double cond();

private:
    std::vector<real> F;    // вектор правой части
    real normF;             // норма вектора правой части
    std::vector<real> xk;   // вектор неизвестных
    std::vector<real> xk1;  // вектор неизвестных (на k+1-й итерации)
    size_t maxIter;		    // макс число итераций
    real eps;			    // точность решения
    size_t blockSize;       // размер блока

    void readParameters(std::string);
    void readVector(std::string);
    real norm(std::vector<real> v);
    real relative_discrepancy();
};