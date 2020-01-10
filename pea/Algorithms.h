#pragma once
#include "Instance.h"
#include "utils.h"
#include <vector>
#include <bitset>
#include <queue>
#include <algorithm>

class BruteForce {
public:
	void run();
	BruteForce(Instance& instance);
	int getFinalDistance();
	vector<int> getFinalSolution();
private:
	Instance &instance;
	int finalDistance = INT_MAX;
	vector<int> finalSolution;
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
	int dp(int i, int state, vector<vector<int>>& cache, vector<vector<int>>& prev);
	Instance& instance;
	int finalDistance = INT_MAX;
	vector<int> finalSolution;
};

struct Node {
	vector<int> path;
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

class BranchNBound {
public:
	void run();
	BranchNBound(Instance& instance, int START_NODE);
	int getFinalDistance();
	vector<int> getFinalSolution();
	void clearQueue();
private:
	priority_queue<Node*, vector<Node*>, comp> queue;
	int size;
	int START_NODE;
	int** data;
	Instance& instance;
	Node* newNode(vector<vector<int>> parentMatrix, vector<int> const& path, int level, int i, int j);
	void rowReduction(vector<vector<int>> &reducedMatrix, vector<int> &row);
	void columnReduction(vector<vector<int>> &reducedMatrix, vector<int> &column);
	int calculateCost(vector<vector<int>> &reducedMatrix);
	void savePath(vector<int> const& path);
	int finalDistance = INT_MAX;
	vector<int> finalSolution;
};

class SimulatedAnnealing {
public:
	double INITIAL_TEMPERATURE;
	double MINIMAL_TEMPTERATURE;
	double TEMPERATURE_FACTOR;
	int NEIGHBOURHOOD;
	int STARTING_SOLUTION;
	void run();
	SimulatedAnnealing(Instance& instance, double initial_temperature, double minimal_temperature, double temperature_factor, int neighbourhood, int starting_solution);
	int getFinalDistance();
	vector<int> getFinalSolution();
	bool swap(vector<int>& solution, int i, int j);
	bool insert(vector<int>& solution, int i, int j);
	bool invert(vector<int>& solution, int i, int j);
private:
	Instance& instance;
	int finalDistance = -1;
	vector<int> finalSolution;
	vector<vector<int>> solutionMemory;
	double cooldown(double T, int t);
	vector<int> generateStartingSolution();
	vector<int> generateRandomSolution();
};

class TabuSearch {
public:
	int NEIGHBOURHOOD;
	int NUMBER_OF_ITERATIONS;
	int TABU_LENGTH;
	int STARTING_SOLUTION;
	void run();
	TabuSearch(Instance& instance, int neighbourhood, int number_of_iterations, int tabu_length, int starting_solution);
	int getFinalDistance();
	vector<int> getFinalSolution();
	bool swap(vector<int>& solution, int i, int j);
	bool insert(vector<int>& solution, int i, int j);
	bool invert(vector<int>& solution, int i, int j);
private:
	Instance& instance;
	int finalDistance = -1;
	vector<int> finalSolution;
	vector<vector<int>> tabuList;
	vector<int> getBestNearbySolution(vector<int> solution, int it, int optimalDistance);
	bool neighbouringSolution(vector<int> &solution, int i, int j);
	vector<int> generateStartingSolution();
	vector<int> generateRandomSolution();
};

class Genetic {
public:
	int POPULATION_SIZE; // rozmiar populacji
	double Pc, Pm; // p-stwa krzy¿owania i mutacji
	int SELECTION_METHOD, CROSS_OPERATOR, MUTATION_OPERATOR;
	void run();
	vector<int> getFinalSolution();
	int getFinalDistance();
	Genetic(Instance& instance, int POPULATION_SIZE, double Pc, double Pm, int SELECTION_METHOD, int CROSS_OPERATOR, int MUTATION_OPERATOR);
private:
	Instance& instance;
	int finalDistance = -1;
	vector<int> finalSolution;
	vector<vector<int>> population;
	vector<vector<int>> matingPool;

	void selectMatingPool();
	double fitnessFunction(vector<int> x);

	vector<vector<int>> crossMatingPool();
	vector<vector<int>> crossPair(vector<int> p, vector<int> q);
	vector<vector<int>> _corss_PMS(vector<int> p, vector<int> q);
	vector<vector<int>> _cross_OX(vector<int> p, vector<int> q);
	vector<vector<int>> _cros_EX(vector<int> p, vector<int> q);

	void mutate(vector<int> &x);
	void mutateGeneration(vector<vector<int>> &g);
	void _mutate_inversion(vector<int>& x);
	void _mutate_insertion(vector<int>& x);
	void _mutate_displacement(vector<int>& x);
	void _mutate_transposition(vector<int>& x);
};