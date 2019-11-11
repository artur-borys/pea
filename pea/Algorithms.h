#pragma once
#include "Instance.h"
#include "utils.h"
#include <vector>
#include <bitset>

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