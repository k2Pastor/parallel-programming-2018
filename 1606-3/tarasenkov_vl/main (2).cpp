#include "SequentialRealisation.h"
#include "ParallelRealisation.h"
#include <mpi.h>
#include <ctime>
#include <iostream>
#include <cmath>

int main(int argc, char* argv[])
{
	int ProcNum, ProcRank;
	double timeStartOfSequential, timeStartOfParallel;
	double timeEndOfSequential, timeEndOfParallel;
	int size;
	double *seqMatrix = nullptr, *parMatrix = nullptr;
	double *seqResult = nullptr, *parResult = nullptr;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0)
	{
		std::cout << "Add matrix size: " << std::endl;
		std::cin >> size;
		seqMatrix = new double[size * (size + 1)];
		srand((unsigned)time(NULL));


		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size + 1; j++)
			{
				seqMatrix[(size + 1) * i + j] = rand() % 100 + 50;
			}
		}

		if (size < 11)
		{
			for (int i = 0; i < size; i++)
			{
				for (int j = 0; j < size + 1; j++)
				{
					std::cout << seqMatrix[(size + 1) * i + j] << " ";
				}

				std::cout << std::endl;
			}
		}

		parMatrix = new double[size * (size + 1)];
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size + 1; j++)
			{
				parMatrix[i * (size + 1) + j] = seqMatrix[i * (size + 1) + j];
			}
		}

		seqResult = new double[size];
		timeStartOfSequential = MPI_Wtime();
		SequentialAlgorithm(seqMatrix, seqResult, size);
		timeEndOfSequential = MPI_Wtime();

		std::cout.precision(5);
		std::cout << "Time of the sequential algorithm is: " << std::fixed << timeEndOfSequential - timeStartOfSequential << std::endl;

	}

	timeStartOfParallel = MPI_Wtime();
	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	parResult = new double[size];
	ParallelAlgorithm(parMatrix, parResult, size, ProcNum, ProcRank);
	timeEndOfParallel = MPI_Wtime();

	if (ProcRank == 0)
	{
		std::cout.precision(5);
		std::cout << "Time of the parallel algorithm is: " << std::fixed << timeEndOfParallel - timeStartOfParallel << std::endl;
		std::cout << "Acceleration is :" << ((timeEndOfSequential - timeStartOfSequential) - (timeEndOfParallel - timeStartOfParallel)) << std::endl;

		bool f = true;

		for (int i = 0; i < size; i++)
		{
			if (round(parResult[i] * 1000) / 1000 != round(seqResult[i] * 1000) / 1000)
			{
				f = false;
			}
		}

		if (f == true)
		{
			std::cout << "Solutions are equal" << std::endl;
		}
		else {

			std::cout << "Solutions are not equal" << std::endl;
		}
	}

	
	MPI_Finalize();
	return 0;
}