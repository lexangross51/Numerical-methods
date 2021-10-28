#include "slau.h"

typedef double real;

int main()
{
    enum { Jacobi, Gauss_Seidel };

    slau s;
    s.readSLAU();

    //for (double w = 0; w < 2 ; w += 0.1)
    {
        s.setInitialApprox();

        //s2.setInitialApprox();
        //s.solve_by_iterative(w, Gauss_Seidel);
        s.BR_Iterations(1.58);

        s.writeResult("x.txt");
    }

    return 0;
}