#pragma once
#include <cmath>
#include <iostream>


void SequentialAlgorithm(double* seqMatrix, double* result, int size)
{
	int Iteration;
	int PivotRow;
	int* pPivotPosition; // Номера строк матрицы, выбираемых в качестве ведущих по итерациям прямого хода
	int* pPivotIteration; // Номера итераций прямого хода, на которых строки выбирались в качестве ведущих

	pPivotPosition = new int[size];
	pPivotIteration = new int[size];

	for (int i = 0; i < size; i++)
	{
		pPivotPosition[i] = 0;
		pPivotIteration[i] = -1;
	}

	for (Iteration = 0; Iteration < size; Iteration++)
	{
		double MaxValue = 0.0;
		for (int i = 0; i < size; i++)
		{
			if ((pPivotIteration[i] == -1) && (fabs(seqMatrix[i * (size + 1) + Iteration]) > MaxValue))
			{
				PivotRow = i;
				MaxValue = seqMatrix[i * (size + 1) + Iteration];
			}
		}

		pPivotPosition[Iteration] = PivotRow;
		pPivotIteration[PivotRow] = Iteration;

		double PivotValue, PivotFactor;
		PivotValue = seqMatrix[PivotRow * (size + 1) + Iteration];

			for (int i = 0; i < size; i++)
			{
				if (pPivotIteration[i] == -1)
				{
					PivotFactor = seqMatrix[i * (size + 1) + Iteration] / PivotValue;
					for (int j = Iteration; j < size; j++)
					{
						seqMatrix[i * (size + 1) + j] -= PivotFactor * seqMatrix[PivotRow * (size + 1) + j];
					}

					seqMatrix[i * (size + 1) + size] -= PivotFactor * seqMatrix[PivotRow * (size + 1) + size];
				}
			}
	}

	int RowIndex, Row;
	for (int i = size - 1; i >= 0; --i)
	{
		RowIndex = pPivotPosition[i];
		result[i] = seqMatrix[RowIndex * (size + 1) + size] / seqMatrix[RowIndex * (size + 1) + i];
		seqMatrix[RowIndex * (size + 1) + i] = 1;

		for (int j = 0; j < i; ++j)
		{
			Row = pPivotPosition[j];
			seqMatrix[Row * (size + 1) + size] -= seqMatrix[Row * (size + 1) + i] * result[i];
			seqMatrix[Row * (size + 1) + i] = 0;
		}
	}

	if (size < 11)
	{
		std::cout << "\nResult is: " << std::endl;
		for (int i = 0; i < size; i++)
		{
			std::cout.precision(5);
			std::cout << std::fixed << (result[i]) << " ";
			std::cout << std::endl;
		}
	}
}