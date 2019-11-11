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

struct Node {
	vector<pair<int, int>> path;
	vector<vector<int>> reducedMatrix;
	int cost;
	int vertex;
	int level;
};

struct comp {
	bool operator()(const Node* lhs, const Node* rhs) const {
		return lhs->cost > rhs->cost;
	}
};

void BranchNBound::run()
{
	priority_queue<Node*, vector<Node*>, comp> queue;

	vector<pair<int, int>> v;

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
			min->path.push_back(make_pair(i, 0));
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

Node* BranchNBound::newNode(vector<vector<int>> parentMatrix, vector<pair<int, int>> const& path, int level, int i, int j)
{
	Node* node = new Node;
	node->path = path;
	if (level != 0) {
		node->path.push_back(make_pair(i, j));
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

void BranchNBound::savePath(vector<pair<int, int>> const& path)
{
	for (int i = path.size() - 1; i >= 0; i--) {
		finalSolution.push_back(path[i].second);
		//finalSolution.push_back(path[i].second);
	}
}
