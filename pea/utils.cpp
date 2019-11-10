#include "pch.h"
#include "utils.h"

int * utils::randomSolution(int size)
{
#if _DEBUG
	cout << "Generator losowych rozwiazan" << endl;
#endif
	int *solution = new int[size];
	for (int i = 0; i < size; i++) {
		solution[i] = -1;
	}
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	uniform_int_distribution<int> distribution(0, size - 1);
	int i = 0;
	while (i < size) {
		int nextIndex = distribution(generator);
		bool solutionAlreadyContains = false;
#if _DEBUG
		cout << "Wylosowano " << nextIndex << endl;
#endif
		for (int j = 0; j < size; j++) {
			if (solution[j] == -1) {
				break;
			}
			if (solution[j] == nextIndex) {
#if _DEBUG
				cout << "Rozwiazanie juz zawiera " << nextIndex << endl;
#endif
				solutionAlreadyContains = true;
				break;
			}
		}
		if (!solutionAlreadyContains) {
			solution[i] = nextIndex;
			i++;
		}
#if _DEBUG
		cout << "Aktualna zawartosc rozwiazania: ";
		for (int j = 0; j < size; j++) {
			cout << solution[j] << ' ';
		}
		cout << endl;
#endif
	}
	return solution;
}

int * utils::ascendingSolution(int size)
{
	int *solution = new int[size];
	for (int i = 0; i < size; i++) {
		solution[i] = i;
	}
	return solution;
}

bool utils::nextPermutation(int* first, int* last) {
	if (first == last) return false;
	int* i = last;
	if (first == --i) return false;

	while (true) {
		int* i1;
		int* i2;
		i1 = i;

		if (*--i < *i1) {
			i2 = last;
			while (!(*i < *(--i2)));
			std::iter_swap(i, i2);
			std::reverse(i1, last);
			return true;
		}
		if (i == first) {
			std::reverse(first, last);
			return false;
		}
	}
}

void utils::printSolution(int* solution, int size)
{
	std::cout << solution[0];
	for (int i = 1; i < size; i++) {
		std::cout << ", " << solution[i];
	}
	std:cout << endl;
}

void utils::printSolution(vector<int> solution)
{
	cout << solution.at(0);
	for (int i = 1; i < solution.size(); i++) {
		cout << ", " << solution.at(i);
	}
	cout << endl;
}
