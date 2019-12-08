#include "pch.h"
#include "Algorithms.h"

void BruteForce::run()
{
	vector<int> solution;
	for (int i = 0; i < instance.getSize(); i++) {
		solution.push_back(i);
	}
	finalSolution.resize(instance.getSize());
	int** data = instance.getData();
	int distance;
	do {
		distance = instance.calculateCostFunction(solution);
		if (distance < finalDistance) {
			finalDistance = distance;
			std::copy(solution.begin(), solution.end(), finalSolution.begin());
		}
	} while (next_permutation(solution.begin() + 1, solution.end()) == 1);
}

BruteForce::BruteForce(Instance& instance) : instance(instance) {
};

int BruteForce::getFinalDistance()
{
	return finalDistance;
}

vector<int> BruteForce::getFinalSolution()
{
	return finalSolution;
}

void DynamicProgramming::run()
{
	int state = 1 << START_NODE;
	vector<vector<int>> cache(size, vector<int>(1 << size, -1));
	vector<vector<int>> prev(size, vector<int>(1 << size, NULL));
	finalDistance = dp(START_NODE, state, cache, prev);

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

int DynamicProgramming::dp(int i, int state, vector<vector<int>> &cache, vector<vector<int>> &prev)
{
	if (state == VISITED_ALL) {
		return data[i][START_NODE];
	}

	if (cache[i][state] != -1) {
		return cache[i][state];
	}

	int minCost = INT_MAX;
	int index = -1;
	for (int next = 0; next < size; next++) {
		if ((state & (1 << next)) != 0) continue;

		int nextState = state | (1 << next);

		int newCost = data[i][next] + dp(next, nextState, cache, prev);
		if (newCost < minCost) {
			minCost = newCost;
			index = next;
		}
	}

	prev[i][state] = index;
	return cache[i][state] = minCost;
}





void BranchNBound::run()
{
	vector<int> v;

	vector<vector<int>> matrix;


	matrix.resize(size);
	for (int i = 0; i < size; i++) {
		matrix[i].resize(size);
		for (int j = 0; j < size; j++) {
			matrix[i][j] = data[i][j];
		}
	}

	Node* root = newNode(matrix, v, 0, -1, 0);

	root->cost = calculateCost(root->reducedMatrix);

	queue.push(root);

	while (!queue.empty()) {
		Node* min = queue.top();

		queue.pop();

		int i = min->vertex;

		if (min->level == size - 1) {
			min->path.push_back(0);
			savePath(min->path);
			finalDistance = min->cost;
			return;
		}

		for (int j = 0; j < size; j++) {
			if (min->reducedMatrix[i][j] != INT_MAX) {
				Node* child = newNode(min->reducedMatrix, min->path, min->level + 1, i, j);
				child->cost = min->cost + min->reducedMatrix[i][j] + calculateCost(child->reducedMatrix);
				queue.push(child);
				
			}
		}
		delete min;
	}
}

BranchNBound::BranchNBound(Instance& instance, int START_NODE) : instance(instance), START_NODE(START_NODE)
{
	size = instance.getSize();
	data = instance.getData();
}

int BranchNBound::getFinalDistance()
{
	return finalDistance;
}

vector<int> BranchNBound::getFinalSolution()
{
	return finalSolution;
}

void BranchNBound::clearQueue()
{
	while (!queue.empty()) {
		Node* min = queue.top();
		queue.pop();
		delete min;
	}
}

Node* BranchNBound::newNode(vector<vector<int>> parentMatrix, vector<int> const& path, int level, int i, int j)
{
	Node* node = new Node;
	node->path = path;
	if (level != 0) {
		node->path.push_back(j);
	}

	node->reducedMatrix = parentMatrix;

	for (int k = 0; level != 0 && k < size; k++) {
		node->reducedMatrix[i][k] = INT_MAX;
		node->reducedMatrix[k][j] = INT_MAX;
	}

	node->reducedMatrix[j][0] = INT_MAX;
	node->level = level;
	node->vertex = j;

	return node;
}

void BranchNBound::rowReduction(vector<vector<int>> &reducedMatrix, vector<int> &row)
{
	fill(row.begin(), row.end(), INT_MAX);

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (reducedMatrix[i][j] < row[i]) {
				row[i] = reducedMatrix[i][j];
			}
		}
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (reducedMatrix[i][j] != INT_MAX && row[i] != INT_MAX) {
				reducedMatrix[i][j] -= row[i];
			}
		}
	}
}

void BranchNBound::columnReduction(vector<vector<int>> &reducedMatrix, vector<int> &column)
{
	fill(column.begin(), column.end(), INT_MAX);

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (reducedMatrix[i][j] < column[j]) {
				column[j] = reducedMatrix[i][j];
			}
		}
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (reducedMatrix[i][j] != INT_MAX && column[j] != INT_MAX) {
				reducedMatrix[i][j] -= column[j];
			}
		}
	}
}

int BranchNBound::calculateCost(vector<vector<int>> &reducedMatrix)
{
	int cost = 0;

	vector<int> row(size);
	rowReduction(reducedMatrix, row);

	vector<int> column(size);
	columnReduction(reducedMatrix, column);

	for (int i = 0; i < size; i++) {
		cost += (row[i] != INT_MAX) ? row[i] : 0;
		cost += (column[i] != INT_MAX) ? column[i] : 0;
	}

	return cost;
}

void BranchNBound::savePath(vector<int> const& path)
{
	for (int i = path.size() - 1; i >= 0; i--) {
		finalSolution.push_back(path[i]);
	}
}

void SimulatedAnnealing::run()
{
	double T = INITIAL_TEMPERATURE;
	int t = 0;
	int size = instance.getSize();
	vector<int> x0;
	int i = utils::random(0, size - 1);
	x0.push_back(i);
	while (x0.size() < size) {
		int min = INT_MAX;
		int id = -1;
		for (int j = 0; j < size; j++) {
			if (i == j) continue;
			if (std::find(x0.begin(), x0.end(), j) != x0.end())	continue;
			int distance = instance.getDistance(i, j);
			if (distance < min) {
				min = distance;
				id = j;
			}
		}
		i = id;
		x0.push_back(id);
	}

	int min_distance = instance.calculateCostFunction(x0);
	vector<int> x = x0;
	while (T >= MINIMAL_TEMPTERATURE) {
		vector<int> y = x;
		int i = 0, j = 0;
		while (i >= j) {
			i = utils::random(0, size - 1);
			j = utils::random(0, size - 1);
		}
		switch(NEIGHBOURHOOD) {
			case 0:
				insert(y, i, j);
				break;
			case 1:
				swap(y, i, j);
				break;
			case 2:
				invert(y, i, j);
				break;
		}

		int distance = instance.calculateCostFunction(y);
		if (distance < min_distance) {
			x = y;
			min_distance = distance;
		}
		else {
			double minProbability = utils::random(0.0, 1.0);
			double acceptanceFactor = exp((double)((long long)min_distance - (long long)distance) / (double)T);
			if (minProbability < acceptanceFactor) {
				x = y;
				min_distance = distance;
			}
		}

		//T = cooldown(T, t);
		T = T - TEMPERATURE_FACTOR;
		t += 1;
	}

	finalDistance = min_distance;
	finalSolution = x;

}

SimulatedAnnealing::SimulatedAnnealing(Instance& instance, double initial_temperature, double minimal_temperature, double temperature_factor, int neighbourhood) :
	instance(instance), INITIAL_TEMPERATURE(initial_temperature), MINIMAL_TEMPTERATURE(minimal_temperature),
	TEMPERATURE_FACTOR(temperature_factor), NEIGHBOURHOOD(neighbourhood)
{
}

int SimulatedAnnealing::getFinalDistance()
{
	return finalDistance;
}

vector<int> SimulatedAnnealing::getFinalSolution()
{
	return finalSolution;
}

bool SimulatedAnnealing::swap(vector<int> &solution, int i, int j)
{
	if (i == j)
		return false;
	std::iter_swap(solution.begin() + i, solution.begin() + j);
	return true;
}

bool SimulatedAnnealing::insert(vector<int>& solution, int i, int j)
{
	if (i == j)
		return false;
	int to_insert = solution.at(i);
	solution.erase(solution.begin() + i);
	solution.insert(solution.begin() + j, to_insert);
	return true;
}

bool SimulatedAnnealing::invert(vector<int>& solution, int i, int j)
{
	if (i == j)
		return false;
	std::reverse(solution.begin() + i, solution.begin() + j + 1);
	return true;
}

double SimulatedAnnealing::cooldown(double T, int t)
{
	double quotient = T / INITIAL_TEMPERATURE;
	if (quotient > 0.75) {
		return T - TEMPERATURE_FACTOR * 100;
	}
	else if (quotient > 0.25) {
		return T - TEMPERATURE_FACTOR * 10;
	}
	else {
		return T - TEMPERATURE_FACTOR;
	}
}

void TabuSearch::run()
{
	int size = instance.getSize();
	vector<int> x0;
	int i = utils::random(0, size - 1);
	x0.push_back(i);
	while (x0.size() < size) {
		int min = INT_MAX;
		int id = -1;
		for (int j = 0; j < size; j++) {
			if (i == j) continue;
			if (std::find(x0.begin(), x0.end(), j) != x0.end())	continue;
			int distance = instance.getDistance(i, j);
			if (distance < min) {
				min = distance;
				id = j;
			}
		}
		i = id;
		x0.push_back(id);
	}

	vector<int> xopt = x0;
	int optimalDistance = instance.calculateCostFunction(xopt);
	tabuList = vector<vector<int>>(size, vector<int>(size, 0));

	for(int i = 0; i < NUMBER_OF_ITERATIONS; i++) {
		vector<int> bestNearbySolution = getBestNearbySolution(x0, i);
		x0 = bestNearbySolution;
		int currentDistance = instance.calculateCostFunction(x0);
		if (currentDistance < optimalDistance) {
			xopt = x0;
			optimalDistance = currentDistance;
		}
	}

	finalSolution = xopt;
	finalDistance = optimalDistance;
}

TabuSearch::TabuSearch(Instance& instance, int neighbourhood, int number_of_iterations, int tabu_length) :
	instance(instance), NEIGHBOURHOOD(neighbourhood), NUMBER_OF_ITERATIONS(number_of_iterations), TABU_LENGTH(tabu_length)
{
}

int TabuSearch::getFinalDistance()
{
	return finalDistance;
}

vector<int> TabuSearch::getFinalSolution()
{
	return finalSolution;
}

bool TabuSearch::swap(vector<int>& solution, int i, int j)
{
	if (i == j)
		return false;
	std::iter_swap(solution.begin() + i, solution.begin() + j);
	return true;
}

bool TabuSearch::insert(vector<int>& solution, int i, int j)
{
	if (i == j)
		return false;
	int to_insert = solution.at(i);
	solution.erase(solution.begin() + i);
	solution.insert(solution.begin() + j, to_insert);
	return true;
}

bool TabuSearch::invert(vector<int>& solution, int i, int j)
{
	if (i == j)
		return false;
	std::reverse(solution.begin() + i, solution.begin() + j + 1);
	return true;
}

vector<int> TabuSearch::getBestNearbySolution(vector<int> solution, int it)
{
	int bestScore = 0;
	int instanceSize = instance.getSize();
	int referenceDistance = instance.calculateCostFunction(solution);
	vector<int> candidateSolution;
	vector<int> bestNearbySolution;
	int bestI = 0, bestJ = 0;
	for (int i = 0; i < instanceSize; i++) {
		candidateSolution = solution;
		for (int j = i + 1; j < instanceSize; j++) {
			neighbouringSolution(candidateSolution, i, j);
			int currentDistance = instance.calculateCostFunction(candidateSolution);
			int currentScore = referenceDistance - currentDistance;
			if (tabuList[i][j] <= it || currentScore > bestScore) {
				bestI = i;
				bestJ = j;
				bestScore = currentScore;
				tabuList[i][j] = it + TABU_LENGTH;
				tabuList[j][i] = it + TABU_LENGTH;
			}
		}
	}

	bestNearbySolution = solution;
	neighbouringSolution(bestNearbySolution, bestI, bestJ);

	return bestNearbySolution;
}

bool TabuSearch::neighbouringSolution(vector<int> &solution, int i, int j)
{
	bool status = false;
	switch (NEIGHBOURHOOD) {
	case 0:
		status = insert(solution, i, j);
		break;
	case 1:
		status = swap(solution, i, j);
		break;
	case 2:
		status = invert(solution, i, j);
		break;
	}

	return status;
}
