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
	int state = 1 << START_NODE;
	vector<vector<int>> memo(size, vector<int>(1 << size, -1));
	vector<vector<int>> prev(size, vector<int>(1 << size, NULL));
	finalDistance = dp(START_NODE, state, memo, prev);

	int index = START_NODE;
	while (true) {
		finalSolution.push_back(index);
		int nextIndex = prev[index][state];
		if (nextIndex == NULL) break;
		int nextState = state | (1 << nextIndex);
		state = nextState;
		index = nextIndex;
	}
	finalSolution.push_back(START_NODE);
}

DynamicProgramming::DynamicProgramming(Instance &instance, int START_NODE) : instance(instance), START_NODE(START_NODE) {
	size = instance.getSize();
	VISITED_ALL = (1 << size) - 1;
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

int DynamicProgramming::dp(int i, int state, vector<vector<int>> &memo, vector<vector<int>> &prev)
{
	if (state == VISITED_ALL) {
		return data[i][START_NODE];
	}

	if (memo[i][state] != -1) {
		return memo[i][state];
	}

	int minCost = INT_MAX;
	int index = -1;
	for (int next = 0; next < size; next++) {
		if ((state & (1 << next)) != 0) continue;

		int nextState = state | (1 << next);

		int newCost = data[i][next] + dp(next, nextState, memo, prev);
		if (newCost < minCost) {
			minCost = newCost;
			index = next;
		}
	}

	prev[i][state] = index;
	return memo[i][state] = minCost;
}
