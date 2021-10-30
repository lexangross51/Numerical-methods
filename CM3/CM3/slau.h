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
    size_t n;					// Размер матрицы
    size_t maxIter;				// Макс число итераций
    real eps;					// Точность решения
    std::vector<size_t> ig;		// Указатели на строки
    std::vector<size_t> jg;		// Указатели на столбцы
    std::vector<real> ggl;		// Нижний треугольник
    std::vector<real> ggu;		// Верхний треугольник
    std::vector<real> di;		// Главная диагональ
    std::vector<real> F;		// Вектор правой части
    real normF;                 // Норма вектора правой части
    std::vector<real> x;		// Вектор неизвестных

    std::vector<real> gglf;		// Нижний треугольник разложенной матрицы
    std::vector<real> gguf;		// Верхний треугольник разложенной матрицы
    std::vector<real> dif;		// Главная диагональ разложенной матрицы

    std::vector<real> r;        // Невязка 
    std::vector<real> z;        // Сопряженное направление
    std::vector<real> p;        // Доп. массив

    std::vector<real> generateVector(size_t size);
    
    void toSparse();

    std::vector<std::vector<real>> mat; // Матрица в плотном формате
};

// Перегрузка оператора сложения векторов
inline std::vector<real> operator+(const std::vector<real>& a, const std::vector<real>& b)
{
    std::vector<real> res = a;
    for (size_t i = 0; i < res.size(); i++)
        res[i] += b[i];
    return res;
}

    // Перегрузка оператора вычитания векторов
    inline std::vector<real> operator-(const std::vector<real>& a, const std::vector<real>& b)
    {
        std::vector<real> res = a;
        for (size_t i = 0; i < res.size(); i++)
            res[i] -= b[i];
        return res;
    }

    // Перегрузка оператора умножения векторов (скалярное произведение)
    inline real operator*(const std::vector<real>& a, const std::vector<real>& b)
    {
        real scalar = 0.0;
        for (size_t i = 0; i < a.size(); i++)
            scalar += a[i] * b[i];
        return scalar;
    }

    // Перегрузка оператора умножения вектора на константу
    inline std::vector<real> operator*(real c, const std::vector<real>& a)
    {
        std::vector<real> res = a;
        for (size_t i = 0; i < res.size(); i++)
            res[i] *= c;
        return res;
    }

#endif