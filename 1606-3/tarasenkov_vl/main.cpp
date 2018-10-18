#include <mpi.h>
#include<iostream>
#include <ctime>
using namespace std;

void InitializeOfMatrix(int *_Matrix, int &_size)
{
	srand((unsigned)time(NULL));

	for ( int i = 0; i < _size * _size; i++ )
	{
		_Matrix[i] = rand() % 100;
	}

}

void OutputOfMatrix(int *_Matrix, int &_size)
{
	for (int i = 0; i < _size; i++)
	{
		for (int j = 0; j < _size; j++)
			cout << _Matrix[i * _size + j] << " ";

		cout << endl;
	}
}

void main(int argc, char **argv)
{
	int *Matrix = nullptr;
	int *resOfSequential = nullptr, *resOfParallel;
	int *Data;
	int size;
	int target_size; // целевой размер, с которым будут работать процессы;
	int ProcNum, ProcRank;
	int minSequential, minParallel;
	int minIndexSequential, minIndexParallel;
	double timeStartOfSequential, timeEndOfSequential, timeStartOfParallel, timeEndOfParallel;
	double Acceleration;

	MPI_Init(&argc, &argv);
	/*Sequential Algorithm*/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum); // получение количества процессов;
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank); // определение ранга процесса;

	

	if (ProcRank == 0)
	{
		cout << "Add Matrix size: ";
		cin >> size;

		Matrix = new int[size * size];
		InitializeOfMatrix(Matrix, size);

		if (size < 11)
			OutputOfMatrix(Matrix, size);

		

		timeStartOfSequential = MPI_Wtime(); // получение текущего момента выполнения программы;
		resOfSequential = new int[size];

		for ( int i = 0; i < size; i++ )
		{
			minSequential = Matrix[i * size];
			minIndexSequential = 0;

			for ( int j = 0; j < size; j++ )
			{
				if (minSequential > Matrix[j + i * size])
				{
					minSequential = Matrix[j + i * size];
					minIndexSequential = j;
				}
			}

			resOfSequential[i] = minIndexSequential;
		}

		timeEndOfSequential = MPI_Wtime();

		if (size < 11)
		{
			for (int i = 0; i < size; i++)
				cout << "String number: " << i << " Minimum " << Matrix[i * size + resOfSequential[i]] << " at "
				<< resOfSequential[i] + 1 << " position" << endl;
		}

	

		cout << endl;
		cout << "TIme of Sequential Algorithm: " << timeEndOfSequential - timeStartOfSequential << endl;
	}
	/*Parallel Algorithm*/
	
		MPI_Barrier(MPI_COMM_WORLD);
		

		MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD); // передача всем процессам размера матрицы;
		target_size = size / ProcNum;
		Data = new int[size * target_size];
		
		timeStartOfParallel = MPI_Wtime();

		if (ProcRank == 0)
		{
			resOfParallel = new int[size];
		}
		else
		{
			resOfParallel = new int[target_size];
		}

		if (target_size > 0)
		{
			MPI_Scatter(Matrix, size * target_size, MPI_INT, Data, size * target_size, MPI_INT, 0, MPI_COMM_WORLD); // рассылка всем процессам соотвествующих блоков;

			for (int i = 0; i < target_size; i++)
			{
				minParallel = Data[i * size];
				minIndexParallel = 0;
				for (int j = 0; j < size; j++)
				{
					if (minParallel > Data[j + i * size])
					{
						minParallel = Data[j + i * size];
						minIndexParallel = j;
					}
				}

				resOfParallel[i] = minIndexParallel;
			}

			for (int i = 1; i < ProcNum; i++) // отправление результатов из всех процессов в 0;
			{
				if ( ProcRank == i )
					MPI_Send(resOfParallel, target_size, MPI_INT, 0, MPI_TAG_UB, MPI_COMM_WORLD);
				if (ProcRank == 0)
					MPI_Recv(&(resOfParallel[i * target_size]), target_size, MPI_INT, i, MPI_TAG_UB, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

			}
		}
		else
		{
			Data = new int[size];
			resOfParallel = new int[1];
		}

		int RestStrings = size - target_size * ProcNum;

		if (ProcRank == 0)
		{
			for (int i = 1; i < RestStrings + 1; i++)
				MPI_Send(&Matrix[(i + target_size * ProcNum - 1) * size], size, MPI_INT, i, MPI_TAG_UB, MPI_COMM_WORLD);
		}
		else if (ProcRank < RestStrings + 1)
		{
			MPI_Recv(Data, size, MPI_INT, 0, MPI_TAG_UB, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

			minParallel = Data[0];
			minIndexParallel = 0;
			for (int i = 0; i < size; i++)
			{
				if (minParallel > Data[i])
				{
					minParallel = Data[i];
					minIndexParallel = i;
				}
			}

			MPI_Send(&minIndexParallel, 1, MPI_INT, 0, MPI_TAG_UB, MPI_COMM_WORLD);
		}

		if (ProcRank == 0) // прием из всех процессов;
		{
			for (int i = 1; i < RestStrings + 1; i++)
			{
				MPI_Recv(&(resOfParallel[size - RestStrings + i - 1]), 1, MPI_INT, i, MPI_TAG_UB, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			}

			if (size < 11)
				for ( int i = 0; i < size; i++ )
				cout << "String number:  " << i << " Minimum: " << Matrix[i * size + resOfParallel[i]] << " at " << resOfParallel[i] + 1
				<< " position" << endl;
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (ProcRank == 0)
		{
			timeEndOfParallel = MPI_Wtime();
			cout << "Time of Parallel Algorithm: " << timeEndOfParallel - timeStartOfParallel << endl;

			bool flag = true;

		/*	for (int i = 0; i < size; i++)
				flag = flag && (resOfSequential[i] == resOfParallel[i]);
			if (flag)
				cout << "Results are equal" << endl;
			else cout << "Results are not equal" << endl; */

			//Acceleration = (timeEndOfParallel - timeStartOfParallel) / (timeEndOfSequential - timeStartOfSequential);
			//cout << "Acceleration: " << Acceleration << endl;
		}

		/*Блок освобождения памяти*/
		if (ProcRank == 0)
			delete[] Matrix;
		delete[] Data;
		delete[] resOfSequential;
	MPI_Finalize();
}
