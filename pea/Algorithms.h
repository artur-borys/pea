#pragma once
#include "Instance.h"
#include "utils.h"
#include <vector>
#include <bitset>
#include <queue>
#include <algorithm>
#include "Tour.h"
#include <chrono>

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
	int SELECTION_METHOD, CROSS_OPERATOR, MUTATION_OPERATOR, ITERATION_COUNT, TOURNAMENT_SIZE, MAX_TIME;
	void run();
	vector<int> getFinalSolution();
	int getFinalDistance();
	Genetic(InstanceVector& instance, int ITERATION_COUNT, int POPULATION_SIZE, double Pc, double Pm, int SELECTION_METHOD, int CROSS_OPERATOR, int MUTATION_OPERATOR);
private:

	InstanceVector& instance;
	int finalDistance = -1;
	vector<int> finalSolution;
	vector<Tour> population;
	vector<Tour> matingPool;

	vector<Tour> generatePopulation();

	Tour selection();
	Tour _selection_roulette();
	Tour _selection_tournament();

	vector<Tour> crossPair(Tour p, Tour q);
	vector<Tour> _cross_PMX(Tour p, Tour q);
	vector<Tour> _cross_OX(Tour p, Tour q);

	void mutate(Tour &x);
	void _mutate_inversion(Tour& x);
	void _mutate_insertion(Tour& x);
	void _mutate_scramble(Tour& x);
	void _mutate_transposition(Tour& x);
};

class Ant {
public:
	vector<bool> tabuList;
	Ant() {};
	int starting_city;
	int size;
	Ant(int starting_city, int size);
	void resetTabu();
	int getTabu(int pos);
	void setTabu(int pos);
};

class AntColony {
public:
	int ANTS_COUNT, ITER_COUNT, STRATEGY, MAX_TIME;
	double feromone_quantity, alpha, beta, initial_feromone = 0.5;
	double vaporating_factor;
	AntColony(InstanceVector& instance) : instance(instance) {};
	vector<vector<double>> feromone;
	void run();
	int selectNextCity(Ant& ant, int city);
	double delta_feromone(int i, int j);
	vector<int> getFinalSolution();
	int getFinalDistance();
private:
	InstanceVector& instance;
	vector<int> finalSolution;
	double parametersCalculate(int i, int j);
	int finalDistance;
};