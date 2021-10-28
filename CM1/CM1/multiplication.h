#pragma once

void multiply(real* ia, real* di, real* al, real* au, real *&res, size_t n)
{
	res = new real[n];
	real *v = new real[n];
	for (size_t i = 0; i < n; i++)
		v[i] = i + 1;


	for (size_t i = 0; i < n; i++)
		res[i] = 0.0;

	for (size_t i = 0; i < n; i++)
	{
		res[i] = di[i] * v[i];
		size_t l = ia[i + 1] - ia[i];
		size_t k = i - l;
		for (size_t j = 0; j < l; j++)
		{
			res[i] += al[size_t(ia[i] + j - 1)] * v[k];
			res[k] += au[size_t(ia[i] + j - 1)] * v[i];
			k++;
		}
	}
}
