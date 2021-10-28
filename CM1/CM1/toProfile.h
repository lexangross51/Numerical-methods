#pragma once
#include <iostream>
#include <fstream>
#include <windows.h>

void readd(real*& mat, size_t n)
{
	std::ifstream fileMatrix("hilbert.txt");

	mat = new real[n * n];

	if (fileMatrix.is_open())
	{
		for (size_t i = 0; i < n; i++)
		{
			for (size_t j = 0; j < n; j++)
			{
				fileMatrix >> mat[i * n + j];
			}
		}
	}
	fileMatrix.close();
}

void toProfile(real* mat, real * a[4], size_t n)
{
	size_t count = 0;
	bool first = true;
	size_t alSize = 0;
	size_t k = 0;

	a[0] = new real[n + 1];
	a[1] = new real[n];

	a[0][0] = a[0][1] = 1;

	for (size_t i = 0; i < n; i++)
	{
		a[1][i] = mat[i * n + i];
		for (size_t j = 0; j < i; j++)
		{
			if (mat[i * n + j] == 0.0 && first && mat[j * n + i] != 0.0 ||
				mat[i * n + j] != 0.0 && first && mat[j * n + i] == 0.0 ||
				mat[i * n + j] != 0.0 && first && mat[j * n + i] != 0.0)
			{
				count++;

				first = false;
			}
			else if (!first)
				count++;
		}
		a[0][i + 1] = a[0][i] + count;
		count = 0;
		first = true;
	}

	a[2] = new real[size_t(a[0][n] - 1)];
	a[3] = new real[size_t(a[0][n] - 1)];

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (mat[i * n + j] == 0.0 && first && mat[j * n + i] != 0.0 ||
				mat[i * n + j] != 0.0 && first && mat[j * n + i] == 0.0 ||
				mat[i * n + j] != 0.0 && first && mat[j * n + i] != 0.0)
			{
				a[2][k] = mat[i * n + j];
				a[3][k] = mat[j * n + i];
				k++;

				first = false;
			}
			else if (!first)
			{
				a[2][k] = mat[i * n + j];
				a[3][k] = mat[j * n + i];
				k++;
			}
		}
		first = true;
	}
}