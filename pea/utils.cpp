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