#include <iostream>
#include <windows.h>
#include "snu.h"

using namespace std;

const string directory = "C:\\Users\\lexan\\OneDrive\\������� ����\\����\\3 ����\\��������� ������\\CM4\\CM4\\";

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
                cout << "������ ������ ���� ������� ������. ������������ ������� ����� 0." << endl;
                break;
            }
            case snu::exceptions::CANT_SOLVE:
            {
                cout << "���������� ������ ����. ����� ����������� ������ ����� ���������." << endl;
                break;
            }
        }
    }

    return 0;
}
