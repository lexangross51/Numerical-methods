#ifndef SNU_H
#define SNU_H

#include <vector>
#include <string>
#include <fstream>

typedef double real;
typedef std::vector<real> Vector;
typedef std::vector<std::vector<real>> Matrix;

class snu
{

public:
	snu(const std::string directory);

	real f(size_t funcNumber);
	void f();
	real df(size_t funcNumber, size_t varNumber);

	void createJacobiMatrix();

private:
	size_t m;		// ���-�� ���������
	size_t n;		// ���-�� �����������
	size_t maxIter;	// ������������ ����� ��������
	real eps1;		// �������� ��� �����
	real eps2;		// �������� ��� �������

	Vector F;		// ������ ������ �����
	Vector x;		// ������ �����������
	real beta;		// ��������� �����

	Matrix A;		// ������� �����

	std::vector<size_t> minF;	// ������ � ���������, �������
								// ������������ ����������� F
	real norm(Vector &v);

	void findMinF();

	bool isInMinF(size_t equationNum);
};

#endif

