#include <iostream>
#include <windows.h>
#include "snu.h"

using namespace std;

const string directory = "C:\\Users\\lexan\\OneDrive\\Рабочий стол\\НГТУ\\3 курс\\Численные методы\\CM4\\CM4\\";

int main()
{
    setlocale(LC_ALL, "");

    try {
        snu s(directory);

        s.newton(snu::tests::TEST5, snu::methods::SYMMETRIZATION);

        return 0;
    }
    catch (snu::exceptions error)
    {
        switch (error)
        {
            case snu::exceptions::VERY_SMALL_DIAG:
            {
                cout << "Нельзя решить СЛАУ методом Гаусса. Диагональынй элемент равен 0." << endl;
                break;
            }
            case snu::exceptions::CANT_SOLVE:
            {
                cout << "Невозможно решить СЛАУ. Число неизвестных больше числа уравнений." << endl;
                break;
            }
        }
    }

    return 0;
}
