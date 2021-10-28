#ifndef SLAU_H
#define SLAU_H

#include <string>
#include <vector>

typedef double real;

class slau
{
public:
    slau(size_t _n, size_t _maxIter, real _eps, real _normF) :
        n(_n), maxIter(_maxIter), eps(_eps), normF(_normF) {};

    slau(const std::string directoryName);

    ~slau() {};

    void toDense(const std::string _dense);

    void writeResult(const std::string file);

    void setInitialApprox();

    std::vector<real> mult(const std::vector<real>& v);
    std::vector<real> multT(const std::vector<real>& v);
    std::vector<real> multD(const std::vector<real>& v);

    void LU();
    std::vector<real> LUDirect(const std::vector<real>& b);
    std::vector<real> LUReverse(const std::vector<real>& y);

    std::vector<real> LUDirectT(const std::vector<real>& b);
    std::vector<real> LUReverseT(const std::vector<real>& y);

    void diagFact();

    void LOS(real* residual, size_t* itersCount);
    void LOSLU(real* residual, size_t* itersCount);
    void LOSdiag(real* residual, size_t* itersCount);

    void MSGLU(real* residual, size_t* itersCount);
    void MSGdiag(real* residual, size_t* itersCount);
    
    void createHilbertMatrix(size_t size);

private:
    size_t n;					// ������ �������
    size_t maxIter;				// ���� ����� ��������
    real eps;					// �������� �������
    std::vector<size_t> ig;		// ��������� �� ������
    std::vector<size_t> jg;		// ��������� �� �������
    std::vector<real> ggl;		// ������ �����������
    std::vector<real> ggu;		// ������� �����������
    std::vector<real> di;		// ������� ���������
    std::vector<real> F;		// ������ ������ �����
    real normF;                 // ����� ������� ������ �����
    std::vector<real> x;		// ������ �����������

    std::vector<real> gglf;		// ������ ����������� ����������� �������
    std::vector<real> gguf;		// ������� ����������� ����������� �������
    std::vector<real> dif;		// ������� ��������� ����������� �������

    std::vector<real> r;        // ������� 
    std::vector<real> z;        // ����������� �����������
    std::vector<real> p;        // ���. ������

    std::vector<real> generateVector(size_t size);
    
    void toSparse();

    std::vector<std::vector<real>> mat; // ������� � ������� �������
};

// ���������� ��������� �������� ��������
inline std::vector<real> operator+(const std::vector<real>& a, const std::vector<real>& b)
{
    std::vector<real> res = a;
    for (size_t i = 0; i < res.size(); i++)
        res[i] += b[i];
    return res;
}

    // ���������� ��������� ��������� ��������
    inline std::vector<real> operator-(const std::vector<real>& a, const std::vector<real>& b)
    {
        std::vector<real> res = a;
        for (size_t i = 0; i < res.size(); i++)
            res[i] -= b[i];
        return res;
    }

    // ���������� ��������� ��������� �������� (��������� ������������)
    inline real operator*(const std::vector<real>& a, const std::vector<real>& b)
    {
        real scalar = 0.0;
        for (size_t i = 0; i < a.size(); i++)
            scalar += a[i] * b[i];
        return scalar;
    }

    // ���������� ��������� ��������� ������� �� ���������
    inline std::vector<real> operator*(real c, const std::vector<real>& a)
    {
        std::vector<real> res = a;
        for (size_t i = 0; i < res.size(); i++)
            res[i] *= c;
        return res;
    }

#endif