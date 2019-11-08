#include "pch.h"
#include "Algorithm.h"

void Algorithm::BruteForce(Instance& instance)
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

Algorithm::Algorithm()
{
	finalDistance = INT_MAX;
	int* finalSolution = nullptr;
}
