#include "pch.h"
#include "Algorithms.h"

void BruteForce::run()
{
	int* solution = new int[instance.getSize()]();
	finalSolution = new int[instance.getSize()]();
	solution = utils::ascendingSolution(instance.getSize());
	int distance;
	do {
		distance = instance.calculateCostFunction(solution);
		if (distance < finalDistance) {
			finalDistance = distance;
			std::copy(solution, &solution[instance.getSize()], finalSolution);
		}
	} while (utils::nextPermutation(&solution[1], &solution[instance.getSize()]) == 1);
}

BruteForce::BruteForce(Instance& instance) : instance(instance) {
	finalSolution = new int[instance.getSize()]();
};

int BruteForce::getFinalDistance()
{
	return finalDistance;
}

int* BruteForce::getFinalSolution()
{
	return finalSolution;
}

void DynamicProgramming::run()
{
	finalSolution.push_back(0);
	finalDistance = dp(1, 0);
}

DynamicProgramming::DynamicProgramming(Instance &instance) : instance(instance) {
	size = instance.getSize();
	VISITED_ALL = (1 << size) - 1;
	cache = new int*[VISITED_ALL];
	for (int i = 0; i < VISITED_ALL + 1; i++) {
		cache[i] = new int[size];
	}
	for (int i = 0; i < VISITED_ALL + 1; i++) {
		for (int j = 0; j < size; j++) {
			cache[i][j] = -1;
		}
	}
	
	data = instance.getData();
};

int DynamicProgramming::getFinalDistance()
{
	return finalDistance;
}

vector<int> DynamicProgramming::getFinalSolution()
{
	return finalSolution;
}

int DynamicProgramming::dp(int mask, int position)
{
	if (mask == VISITED_ALL) {
		return data[position][0];
	}

	if (cache[mask][position] != -1) {
		return cache[mask][position];
	}

	int distance = INT_MAX;
	int city;

	for (int i = 0; i < size; i++) {
		if ((mask & (1 << i)) == 0) {
			int new_mask = mask | (1 << i);
			int newDistance = data[position][i] + dp(new_mask, i);
			if (newDistance < distance) {
				distance = newDistance;
				city = i;
			}
		}
	}

	return cache[mask][position] = distance;
}
