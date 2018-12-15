#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <ctime>
using namespace std;
#define endLine cout << endl
void quicksort(int *mas, int first, int last)
{
	int mid, count;
	int f = first, l = last;
	mid = mas[(f + l) / 2];
	do
	{
		while (mas[f]<mid) f++;
		while (mas[l]>mid) l--;
		if (f <= l)
		{
			count = mas[f];
			mas[f] = mas[l];
			mas[l] = count;
			f++;
			l--;
		}

	} while (f<l);
	if (first<l) quicksort(mas, first, l);
	if (f<last) quicksort(mas, f, last);

}
int *merge(int *a, int n, int *b, int m)
{
	int *result = new int[n + m];
	int i = 0, j = 0;
	int index = 0;

	while (i<n && j<m)
	{
		if (a[i] < b[j])
		{
			result[index] = a[i];
			i++;
		}
		else
		{
			result[index] = b[j];
			j++;
		}

		index++;
	}

	while (i < n)
	{
		result[index] = a[i];
		index++;
		i++;
	}

	while (j < m)
	{
		result[index] = b[j];
		index++;
		j++;
	}

	return result;
}

void main(int argc, char *argv[])
{
	int Size = 1000000;
	int *mass , *sendCounts = nullptr, *displs = nullptr, *paralMass = nullptr, *rez = nullptr;
	int first, last, paralPart;
	double timeStartOfSequential, timeEndOfSequential, timeStartOfParallel, timeEndOfParallel;
	int ProcNum, ProcRank;
	MPI_Status status;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int s = ProcNum;


	if (argc > 1)
		Size = atoi(argv[1]);

	mass = new int[Size];

	if (ProcRank == 0)
	{
		for (int i = 0; i < Size; i++)
			mass[i] = rand() % 100;

		if (Size < 11)
		{
			cout << "Array:" << endl;
			for (int i = 0; i < Size; i++)
			{
				cout << mass[i] << " ";
			}
		}
	}
	// Массивы sendCounts и displs для передачи массива через Scatterv 
	if (ProcNum>1)
	{
		paralPart = Size / (ProcNum - 1);

		int given = 0;
		sendCounts = new int[ProcNum];
		sendCounts[0] = 0;
		for (int i = 1; i<ProcNum - 1; i++)
		{
			sendCounts[i] = paralPart;
			given += paralPart;
		}
		sendCounts[ProcNum - 1] = Size - given;

		displs = new int[ProcNum];
		for (int i = 0; i<ProcNum; i++)
			displs[i] = 0;
		for (int i = 2; i<ProcNum; i++)
			displs[i] += (i - 1)*sendCounts[i - 1];

		paralMass = new int[sendCounts[ProcRank]];
		for (int i = 0; i<sendCounts[ProcRank]; i++)
			paralMass[i] = 0;

		rez = new int[Size];
		for (int i = 0; i<Size; i++)
			rez[i] = 0;
		

		MPI_Scatterv(mass, sendCounts, displs, MPI_INT, paralMass, sendCounts[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);
	}


	//Последовательный алгоритм
	if (ProcRank == 0)
	{
		timeStartOfSequential = MPI_Wtime() * 1000;
		quicksort(mass, 0, Size - 1);
		timeEndOfSequential = MPI_Wtime() * 1000;

		/*for(int i=0; i<Size; i++)
		cout<<mass[i]<<" ";
		endLine;*/

		if (ProcNum>1)
		{
			MPI_Send(mass, Size, MPI_INT, ProcNum - 1, 8, MPI_COMM_WORLD);
			MPI_Send(&timeStartOfSequential, 1, MPI_DOUBLE, ProcNum - 1, 9, MPI_COMM_WORLD);
			MPI_Send(&timeEndOfSequential, 1, MPI_DOUBLE, ProcNum - 1, 10, MPI_COMM_WORLD);
		}
		endLine;
		cout << "Time of sequential algorithm: " << timeEndOfSequential - timeStartOfSequential << "ms" << endl;
		
	}
	else
	{
		if (ProcRank == ProcNum - 1)
			timeStartOfParallel = MPI_Wtime() * 1000;

		quicksort(paralMass, 0, sendCounts[ProcRank] - 1);
		/*
		for(int i=0; i<sendCounts[ProcRank]; i++)
		cout<<paralMass[i]<<" ";
		endLine;
		*/
		int m = 1;

		int k = sendCounts[ProcRank];

		while (s > 1)
		{

			s = s / 2 + s % 2;
			
			if ((ProcRank - 1 - m) % (2 * m) == 0)
			{
				MPI_Send(&k, 1, MPI_INT, ProcRank - m, 0, MPI_COMM_WORLD);
				MPI_Send(paralMass, k, MPI_INT, ProcRank - m, 0, MPI_COMM_WORLD);
				
			}

			if (((ProcRank - 1) % (2 * m) == 0) && (ProcNum - ProcRank > m))
			{
				int k1;
				int *paralMass2;
				MPI_Recv(&k1, 1, MPI_INT, ProcRank + m, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				paralMass2 = new int[k1];
				MPI_Recv(paralMass2, k1, MPI_INT, ProcRank + m, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				paralMass = merge(paralMass, k, paralMass2, k1);
				if (k + k1 == Size)
					MPI_Send(paralMass, Size, MPI_INT, ProcNum - 1, 11, MPI_COMM_WORLD);
				k = k + k1;
			}

			
			m = 2 * m;
				
		}


		if (ProcRank == ProcNum - 1)
		{
			timeEndOfParallel = MPI_Wtime() * 1000;
			MPI_Recv(mass, Size, MPI_INT, 0, 8, MPI_COMM_WORLD, &status);
			MPI_Recv(&timeStartOfSequential, 1, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD, &status);
			MPI_Recv(&timeEndOfSequential, 1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, &status);
			MPI_Recv(rez, Size, MPI_INT, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, &status);

			//Сравниваем поэлементно результаты двух алгоритмов
			int same = 0;
			for (int i = 0; i<Size; i++)
				if (rez[i] == mass[i])
					same++;

			
				if (Size < 11)
				{
					cout << "Sorted array: " << endl;
					for (int i = 0; i < Size; i++)
					{
						cout << mass[i] << " ";
					}
				}
			
				endLine;
			cout << "Time of parallel algorithm: " << timeEndOfParallel - timeStartOfParallel << "ms" << endl;

			if (same == Size)
				cout << "The results of algorithms are the same" << endl;
			else
				cout << "The results of algorithms are different" << endl;

			cout << "Acceleration: " << (timeEndOfSequential - timeStartOfSequential) / (timeEndOfParallel - timeStartOfParallel) << endl;
		}
	}
	MPI_Finalize();
	//if (ProcNum == 1)
		//system("pause");

}