#pragma once
#include <vector>
#include <iostream>
typedef double real;

class matrix
{
public:
    matrix() { m = n = 0; }

    void read(std::string);

    real multLineOnVector(size_t i, const std::vector<real> &v, int mode);

    void mult(const std::vector<real> &v, std::vector<real>& res);

    void toDense(std::vector<std::vector<real>> &mat);

protected:
    std::vector<real> di;		  // ������� ���������
    std::vector<real> l1, l2, l3; // 3 ��������� ��� �������
    std::vector<real> u1, u2, u3; // 3 ��������� ��� �������
    size_t m;					  // ���������� ������� ����������
    size_t n;					  // ������ �������
};