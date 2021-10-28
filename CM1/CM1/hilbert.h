#pragma once
#include <fstream>
#include <iomanip>

real* hilbert_matrix(size_t k)
{
    real* hm = new real[k * k];
    for (size_t i = 0; i < k; i++)
        for (size_t j = 0; j < k; j++)
            hm[i * k + j] = 1.0 / real((i + 1) + (j + 1) - 1);

    std::ofstream out("hilbert.txt");
    if (out.is_open())
    {
        for (size_t i = 0; i < k; i++)
        {
            for (size_t j = 0; j < k; j++)
                out << std::left << std::setw(10) << hm[i * k + j];            
            out << std::endl;
        }
    }
    return hm;
}