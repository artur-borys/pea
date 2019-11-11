#pragma once
#include "Instance.h"
#include "utils.h"
#include <vector>
#include <bitset>
#include <queue>

class BruteForce {
public:
	void run();
	BruteForce(Instance& instance);
	int getFinalDistance();
	int* getFinalSolution();
private:
	Instance &instance;
	int finalDistance = INT_MAX;
	int* finalSolution;
};

class DynamicProgramming {
public:
	void run();
	DynamicProgramming(Instance& instance, int START_NODE);
	int getFinalDistance();
	vector<int> getFinalSolution();
private:
	int size;
	int START_NODE;
	int VISITED_ALL;
	int** data;
	int dp(int i, int state, vector<vector<int>>& memo, vector<vector<int>>& prev);
	Instance& instance;
	int finalDistance = INT_MAX;
	vector<int> finalSolution;
};

struct Node;
struct comp;

class BranchNBound {
public:
	void run();
	BranchNBound(Instance& instance, int START_NODE);
	int getFinalDistance();
	vector<int> getFinalSolution();
private:
	int size;
	int START_NODE;
	int** data;
	Instance& instance;
	Node* newNode(vector<vector<int>> parentMatrix, vector<pair<int, int>> const& path, int level, int i, int j);
	void rowReduction(vector<vector<int>> &reducedMatrix, vector<int> &row);
	void columnReduction(vector<vector<int>> &reducedMatrix, vector<int> &column);
	int calculateCost(vector<vector<int>> &reducedMatrix);
	void savePath(vector<pair<int, int>> const& path);
	int finalDistance = INT_MAX;
	vector<int> finalSolution;
};