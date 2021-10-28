#include "matrix.h"
#include <fstream>

// ������ ������� �� �����
void matrix::read(std::string name)
{
    std::ifstream fin(name);
    if (fin.is_open())
    {		
        fin >> n >> m;

        di.resize(n);
        for (size_t i = 0; i < di.size(); i++)
            fin >> di[i];

        l1.resize(n - 1);
        for (size_t i = 0; i < l1.size(); i++)
            fin >> l1[i];

        l2.resize(n - m - 2);
        for (size_t i = 0; i < l2.size(); i++)
            fin >> l2[i];

        l3.resize(n - m - 3);
        for (size_t i = 0; i < l3.size(); i++)
            fin >> l3[i];

        u1.resize(n - 1);
        for (size_t i = 0; i < u1.size(); i++)
            fin >> u1[i];

        u2.resize(n - m - 2);
        for (size_t i = 0; i < u2.size(); i++)
            fin >> u2[i];

        u3.resize(n - m - 3);
        for (size_t i = 0; i < u3.size(); i++)
            fin >> u3[i];

        fin.close();
    }
}

// ��������� i-� ������ ������� �� ������ v. mode ���������, ��� ����� �������� 
// ������ �� ������: 1-��� ������ �� ������ (��� �����), 2 - �� ������ ������ 
// �� ���������(�� �������), 3 - �� ��������� �� ����� ������ (��� ������-�������)
real matrix::multLineOnVector(size_t i, const std::vector<real> &v, int mode)
{
    real sum = 0.0;

    if (mode == 1 || mode == 2)
    {
        // ���� �� � 1-� ������,
        if (i > 0)
        {
            // �� ����� ��������� ������� 1-� ��������� 
            // �� ������� ������������
            sum += l1[i - 1] * v[i - 1];

            // ���� ����������� � 2-� ������ ���������
            if (i > m + 1)
                sum += l2[i - m - 2] * v[i - m - 2];
            // ���� ����������� � 3-� ������ ���������
            if (i > m + 2)
                sum += l3[i - m - 3] * v[i - m - 3];
        }
    }
    if (mode == 1 || mode == 3)
    {
        // �������� ������������ �� ��������������� 
        // ������� �������
        sum += di[i] * v[i];
        
        // ���� �� � ��������� ������,
        if (i < n - 1)
        {
            // �� ����� ��������� ������� 1-� ���������
            // �������� ������������
            sum += u1[i] * v[i + 1];
            // ���� ����������� � 2-� ������� ���������
            if (i < n - m - 2)
                sum += u2[i] * v[i + m + 2];
            // ���� ����������� � 3-� ������� ���������
            if (i < n - m - 3)
                sum += u3[i] * v[i + m + 3];
        }
    }
    return sum; 
}

// ��������� ������� �� ������ v. 
void matrix::mult(const std::vector<real> &v, std::vector<real>& res)
{
    res.resize(n);
    for (size_t i = 0; i < n; i++)
        res[i] = multLineOnVector(i, v, 1);
}

#pragma region ������� � ������� ������
void matrix::toDense(std::vector<std::vector<real>>& mat)
{
    mat.clear();
    mat.resize(n);

    for (size_t i = 0; i < n; i++)
    {
        mat[i].resize(n, 0);
        mat[i][i] = di[i];
    }

    size_t j = 1;
    for (size_t i = 0; i < l1.size(); i++, j++)
    {
        mat[i][j] = u1[i];
        mat[j][i] = l1[i];
    }

    j = m + 2;
    for (size_t i = 0; i < l2.size(); i++, j++)
    {
        mat[i][j] = u2[i];
        mat[j][i] = l2[i];
    }

    j = m + 3;
    for (size_t i = 0; i < l3.size(); i++, j++)
    {
        mat[i][j] = u3[i];
        mat[j][i] = l3[i];
    }
}

#pragma endregion