#pragma once
#include <mpi.h>
#include <cstddef>
#include <cmath>
#include <iostream>

void ParallelAlgorithm(double* parMatrix, double* result, int size, int ProcNum, int ProcRank)
{
	int PivotRow;
	int* pPivotPositionp = new int[size]; // ������ ����� �������, ���������� � �������� ������� �� ��������� ������� ����
	int* pPivotIterationp = new int[size]; // ������ �������� ������� ����, �� ������� ������ ���������� � �������� �������
	

	for (int i = 0; i < size; i++)
		pPivotIterationp[i] = -1;

	double* pSubMatr = nullptr;

	typedef struct
	{
		int Row;
		double Value;

	}sPivotRow;

	MPI_Datatype TPivot;
	int blocklens[2] = { 1 , 1 }; // ������ ���� ������� �����;
	MPI_Aint indices[2] = { offsetof(sPivotRow, Row), offsetof(sPivotRow, Value) }; // ���, ��������������� ��� �������� �������, 
	// � ����� �������� ����� ���������� �������� � ������.���������� ��� �������� �������������
		// ���, ������ �������� �������� ����������� ��� �������� ������ ������
	MPI_Datatype oldtypes[2] = { MPI_INT, MPI_DOUBLE };
	MPI_Type_create_struct(2, blocklens, indices, oldtypes, &TPivot); // ����� ������, ����� ��������� � ������ �����, 
	// �������� ������� ����� �� ������ ����, �������� ��� ������, ����������� ��� ������;
	MPI_Type_commit(&TPivot); // ����������� ������ ���� ������;

	int k = size / ProcNum;
	int p = k + size % ProcNum;
	int Iteration;
	int *pPivotIterL;

	if (ProcRank == 0)
	{
		pPivotIterL = new int[p];
		for (int i = 0; i < p; i++)
		{
			pPivotIterL[i] = -1;
		}
	}
	else {
		pPivotIterL = new int[k];
		for (int i = 0; i < k; i++)
		{
			pPivotIterL[i] = -1;
		}
	}

	int *send_counts, *displs, *send_countsg, *displsg;
	send_counts = new int[ProcNum]; // ����� ���������, ������������ ������� ��������;
	displs = new int[ProcNum]; // ������ �������� (� ���������) �� ������ ������ ��������;
	send_countsg = new int[ProcNum];
	displsg = new int[ProcNum];

	send_counts[0] = p * (size + 1);
	for (int i = 1; i < ProcNum; i++)
		send_counts[i] = k * (size + 1);

	displs[0] = 0;
	for (int i = 1; i < ProcNum; i++)
		displs[i] = p * (size + 1) + (i - 1) * k * (size + 1);

	send_countsg[0] = p * (size + 1);
	for (int i = 1; i < ProcNum; i++)
		send_countsg[i] = k * (size + 1);

	displsg[0] = 0;
	for (int i = 1; i < ProcNum; i++)
		displsg[i] = p * (size + 1) + (i - 1) * k * (size + 1);

	double* pPivotString = new double[size + 1];
	double* pPivotStr = new double[size + 1];

	if (ProcRank == 0)
	{
		pSubMatr = new double[(size + 1) * p];
	}
	else {
		pSubMatr = new double[(size + 1) * k];
	}

	// ������� MPI_Scatterv �������� ������� �������� ��������� ���������� ���������. ������ ������������ ��������� �����, ����������� i-�� ��������,
	// �������� � ������� �������� displs, � ����� ���������� ��������� -  � ������� sendcounts.
	MPI_Scatterv(parMatrix, send_counts, displs, MPI_DOUBLE, pSubMatr, send_counts[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (Iteration = 0; Iteration < size; Iteration++)
	{
		sPivotRow* recvtmp = new sPivotRow[ProcNum];
		sPivotRow A; 
		A.Row = -1; A.Value = 0;

		if (ProcRank == 0)
		{
			for (int i = 0; i < p; i++)
			{
				if ((pPivotIterL[i] == -1) && fabs(pSubMatr[i * (size + 1) + Iteration]) > A.Value)
				{
					A.Value = pSubMatr[i * (size + 1) + Iteration];
					A.Row = i;
				}
			}
		}
		else {

			for (int i = 0; i < k; i++)
			{
				if ((pPivotIterL[i] == -1) && fabs(pSubMatr[i * (size + 1) + Iteration]) > A.Value)
				{
					A.Value = pSubMatr[i * (size + 1) + Iteration];
					A.Row = p + (ProcRank - 1) * k + i;

				}
			}
		}

		// ���������� ������ ������ ������, ���������� ����� ���������� ������, � ���� ������ �������� � ������� root. 
		// ����� ������ �������������� ����������. ����������� ���������� � ������� ���������� ������� ���������-������������. 
		// �� ���� ������, ��������� ��������� i �� ������ ������ sendbuf, ���������� � i-� ������ ������ recvbuf �������� root. 
		// ����� �������, � ������� ���������� ������, ������ ���� ����������� ��� �� ����������.

		MPI_Gather(&A, 1, TPivot, recvtmp, 1, TPivot, 0, MPI_COMM_WORLD);

		if (ProcRank == 0)
		{
			double MaxValue = 0.0;

			for (int i = 0; i < ProcNum; i++)
			{
				if (fabs(recvtmp[i].Value) > MaxValue && recvtmp[i].Row != -1)
				{
					PivotRow = recvtmp[i].Row;
					MaxValue = recvtmp[i].Value;
				}
			}
		}
		// ������� � ������� root ��������� ��������� �� ������ ������ �������� ���� ��������� ������� ����� ������������� COMM;
		MPI_Bcast(&PivotRow, 1, MPI_INT, 0, MPI_COMM_WORLD); // ����� ������ ������������ � ������ ����������� ������, ����� ���������� ���������,
		// ��� ���������� ���������, ����� �������� - �����������, ������������;

		if (ProcRank == 0)
		{
			pPivotPositionp[Iteration] = PivotRow;
		}

		int r;
		if (PivotRow < p)
		{
			r = 0;
		}
		else {

			r = (PivotRow - p) / k + 1;
		}

		int tmp;

		if (r == 0)
		{
			tmp = PivotRow;
		}
		else {
			tmp = (PivotRow - p - k * (ProcRank - 1)) % k;
		}

		if (ProcRank == r)
		{
			for (int i = 0; i < size + 1; i++)
			{
				pPivotStr[i] = pSubMatr[tmp * (size + 1) + i];
			}
			pPivotIterL[tmp] = Iteration;
		}

		MPI_Bcast(pPivotStr, size + 1, MPI_DOUBLE, r, MPI_COMM_WORLD);

		double PivotValue, PivotFactor;
		PivotValue = pPivotStr[Iteration];

		if (ProcRank == 0)
		{
			for (int i = 0; i < p; i++)
			{
				if (pPivotIterL[i] == -1)
				{
					PivotFactor = pSubMatr[i * (size + 1) + Iteration] / PivotValue;
					for (int j = Iteration; j < size; j++)
					{
						pSubMatr[i * (size + 1) + j] -= PivotFactor * pPivotStr[j];
					}
					pSubMatr[i * (size + 1) + size] -= PivotFactor * pPivotStr[size];
				}
			}
		}
		else {
			for (int i = 0; i < k; i++)
			{
				if (pPivotIterL[i] == -1)
				{
					PivotFactor = pSubMatr[i * (size + 1) + Iteration] / PivotValue;
					for (int j = Iteration; j < size; j++)
					{
						pSubMatr[i * (size + 1) + j] -= PivotFactor * pPivotStr[j];
					}
					pSubMatr[i * (size + 1) + size] -= PivotFactor * pPivotStr[size];
				}
			}
		}
	}

	MPI_Gatherv(pSubMatr, send_countsg[ProcRank], MPI_DOUBLE, parMatrix, send_countsg, displsg, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(pPivotPositionp, size, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = size - 1; i >= 0; --i)
	{
		int strpos = pPivotPositionp[i];
		int rank;

		if (strpos < p)
		{
			rank = 0;
		}
		else {
			rank = (strpos - p) / k + 1;
		}

		if (ProcRank == rank)
		{
			if (ProcRank == 0)
			{
				int strposL = strpos;
				result[i] = pSubMatr[strposL * (size + 1) + size] / pSubMatr[(size + 1) * strposL + i];
				pSubMatr[strposL * (size + 1) + i] = 1;

				for (int m = 0; m < p; m++)
				{
					if (m != strposL && i > pPivotIterL[m])
					{
						pSubMatr[m * (size + 1) + size] -= pSubMatr[m * (size + 1) + i] * result[i];
						pSubMatr[m * (size + 1) + i] = 0;
					}
				}
			}
			else {
				int strposL = (strpos - p - k * (ProcRank - 1)) % k;
				result[i] = pSubMatr[strposL * (size + 1) + size] / pSubMatr[(size + 1) * strposL + i];
				pSubMatr[strposL * (size + 1) + i] = 1;

				for (int m = 0; m < k; m++)
				{
					if (m != strposL && i > pPivotIterL[m])
					{
						pSubMatr[m * (size + 1) + size] -= pSubMatr[m * (size + 1) + i] * result[i];
						pSubMatr[m * (size + 1) + i] = 0;
					}
				}
			}
		}

		MPI_Bcast(&result[i], 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);

		if (ProcRank != rank)
		{
			if (ProcRank == 0)
			{
				for (int m = 0; m < p; m++)
					if (i > pPivotIterL[m])
					{
						pSubMatr[m * (size + 1) + size] -= pSubMatr[m * (size + 1) + i] * result[i];
						pSubMatr[m * (size + 1) + i] = 0;
					}
			}
			else {
				for (int m = 0; m < k; m++)
					if (i > pPivotIterL[m])
					{
						pSubMatr[m * (size + 1) + size] -= pSubMatr[m * (size + 1) + i] * result[i];
						pSubMatr[m * (size + 1) + i] = 0;
					}

			}
		}
	}

	if (ProcRank == 0)
	{
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
}