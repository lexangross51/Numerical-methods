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
    std::vector<real> F;    // ������ ������ �����
    real normF;             // ����� ������� ������ �����
    std::vector<real> xk;   // ������ �����������
    std::vector<real> xk1;  // ������ ����������� (�� k+1-� ��������)
    size_t maxIter;		    // ���� ����� ��������
    real eps;			    // �������� �������
    size_t blockSize;       // ������ �����

    void readParameters(std::string);
    void readVector(std::string);
    real norm(std::vector<real> v);
    real relative_discrepancy();
};