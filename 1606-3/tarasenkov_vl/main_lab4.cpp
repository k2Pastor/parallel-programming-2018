#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <random>
#include <ctime>



size_t GenerateRandomNumber(std::mt19937& generator)
{
	std::uniform_int_distribution<int> distribution(0, 10000);
	return static_cast<size_t>(distribution(generator));
}

void FillArray(int* arr, uint32_t size)
{
	std::mt19937 gen(time(0));
	std::random_device random_device;
	std::mt19937 generator(random_device());

	for (int i = 0; i < size; i++) {
		arr[i] = GenerateRandomNumber(generator);
	}
}

void ShowArray(int* arr, uint32_t size)
{
	for (int i = 0; i < size; i++) {
		std::cout << "Arr[ " << i << " ]" << " = " <<  arr[i] << std::endl;
	}
}

uint32_t BinSearch(int* arr, uint32_t l, uint32_t r, int x) {
	if (l == r) {
		return l;
	}
	if (l + 1 == r) {
		if (x < arr[l]) {
			return l;
		}
		else {
			return r;
		}
	}
	uint32_t m = (l + r) / 2;
	if (x < arr[m]) {
		r = m;
	}
	else {
		if (x > arr[m]) {
			l = m;
		}
		else {
			return m;
		}
	}
	return BinSearch(arr, l, r, x);
}

void ShellSort(int* arr, uint32_t size)
{
	uint32_t i, j, d;
	int t;
	for (d = size / 2; d > 0; d /= 2)
	{
		for (i = d; i < size; i++)
		{
			t = arr[i];
			for (j = i; j >= d; j -= d)
			{
				if (t < arr[j - d]) {
					arr[j] = arr[j - d];
				}
				else break;
			}
			arr[j] = t;
		}
	}
}

void Merge(int* arr1, int* arr2, int* res, uint32_t size1, uint32_t size2)
{
	uint32_t a = 0;
	uint32_t b = 0;
	uint32_t i = 0;

	while ((a != size1) && (b != size2)) {
		if (arr1[a] <= arr2[b]) {
			res[i] = arr1[a];
			a++;
		}
		else {
			res[i] = arr2[b];
			b++;
		}
		i++;
	}

	if (a == size1) {
		uint32_t j = b;
		for (; j < size2; j++, i++) {
			res[i] = arr2[j];
		}
	}
	else {
		uint32_t j = a;
		for (; j < size1; j++, i++) {
			res[i] = arr1[j];
		}

	}
}

void MergeArrays(int* arr1, int* arr2, int* res, uint32_t size1, uint32_t size2)
{
	uint32_t MedIndex = BinSearch(arr2, 0, size2, (arr1[size1 / 2]));
	uint32_t tmp1Size = MedIndex + size1 / 2;
	int* tmp1 = new int[tmp1Size];

	uint32_t tmp2Size = size1 + size2 - tmp1Size;
	int* tmp2 = new int[tmp2Size];

	Merge(arr1, arr2, tmp1, size1 / 2, MedIndex);
	Merge(arr1 + size1 / 2, arr2 + MedIndex, tmp2, size1 - size1 / 2, size2 - (MedIndex));

	std::copy(tmp1, tmp1 + tmp1Size, res);
	std::copy(tmp2, tmp2 + tmp2Size, res + tmp1Size);

	delete[] tmp1;
	delete[] tmp2;
}

void ShellMerge(const int* arr, int* res, uint32_t size)
{
	uint32_t tmpArr1Size = size / 2;
	int* tmpArr1 = new int[tmpArr1Size];
	std::copy(arr, arr + tmpArr1Size, tmpArr1);

	uint32_t tmpArr2Size = size - size / 2;
	int* tmpArr2 = new int[tmpArr2Size];
	std::copy(arr + tmpArr1Size, arr + size, tmpArr2);

	ShellSort(tmpArr1, tmpArr1Size);
	ShellSort(tmpArr2, tmpArr2Size);

	MergeArrays(tmpArr1, tmpArr2, res, tmpArr1Size, tmpArr2Size);

	delete[] tmpArr1;
	delete[] tmpArr2;
}

int main()
{
	uint32_t size;
	std::cout << "Add size: " << std::endl;
	std::cin >> size;
	int* arr = new int[size];
	int* copyArr = new int[size];
	int* resArr = new int[size];

	FillArray(arr, size);

	if (size < 50) {
		ShowArray(arr, size);
	}

	std::copy(arr, arr + size, copyArr);
	std::cout << std::endl;
	std::cout << "Let`s start merging!" << std::endl;
	std::cout << std::endl;
	ShellMerge(arr, resArr, size);
	std::cout << "Let`s start sorting!" << std::endl;
	ShellSort(copyArr, size);

	if (size < 51) {
		ShowArray(copyArr, size);
	}

	return 0;
}