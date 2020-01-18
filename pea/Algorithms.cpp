#include "pch.h"
#include "Algorithms.h"
#include <cassert>

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
	vector<int> x0 = generateStartingSolution();

	int min_distance = instance.calculateCostFunction(x0);
	vector<int> x = x0;
	int repeat_count = 0;
	while (T >= MINIMAL_TEMPTERATURE) {
		while (T >= MINIMAL_TEMPTERATURE && repeat_count < 10) {
			vector<int> y = x;
			int i = 0, j = 0;
			while (i >= j) {
				i = utils::random(0, size - 1);
				j = utils::random(0, size - 1);
			}
			switch (NEIGHBOURHOOD) {
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
			if (distance == min_distance) {
				repeat_count++;
			}
			if (distance < min_distance) {
				x = y;
				min_distance = distance;
				repeat_count = 0;
			}
			else {
				double minProbability = utils::random(0.0, 1.0);
				double acceptanceFactor = exp((double)((long long)min_distance - (long long)distance) / (double)T);
				if (minProbability < acceptanceFactor) {
					x = y;
					min_distance = distance;
					repeat_count = 0;
				}
			}
			T = T * TEMPERATURE_FACTOR;
		}
		x0 = generateRandomSolution();
	}
	finalDistance = min_distance;
	finalSolution = x;

}

vector<int> SimulatedAnnealing::generateRandomSolution()
{
	vector<int> s;
	while (s.size() < instance.getSize()) {
		int i = utils::random(0, instance.getSize() - 1);
		if (std::find(s.begin(), s.end(), i) == s.end()) {
			s.push_back(i);
		}
	}
	return s;
}

SimulatedAnnealing::SimulatedAnnealing(Instance& instance, double initial_temperature, double minimal_temperature, double temperature_factor, int neighbourhood, int starting_solution) :
	instance(instance), INITIAL_TEMPERATURE(initial_temperature), MINIMAL_TEMPTERATURE(minimal_temperature),
	TEMPERATURE_FACTOR(temperature_factor), NEIGHBOURHOOD(neighbourhood), STARTING_SOLUTION(starting_solution)
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
	double quotient = T / (INITIAL_TEMPERATURE - MINIMAL_TEMPTERATURE);
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

vector<int> SimulatedAnnealing::generateStartingSolution()
{
	if (STARTING_SOLUTION == 0) {
		vector<int> s;
		for (int i = 0; i < instance.getSize(); i++) {
			s.push_back(i);
		}
		return s;
	}
	if (STARTING_SOLUTION == 1) {
		vector<int> s;
		while (s.size() < instance.getSize()) {
			int i = utils::random(0, instance.getSize() - 1);
			if (std::find(s.begin(), s.end(), i) == s.end()) {
				s.push_back(i);
			}
		}
		return s;
	}
	if (STARTING_SOLUTION == 2 || STARTING_SOLUTION == 3) {
		int size = instance.getSize();
		vector<int> s;
		int i = 0;
		if (STARTING_SOLUTION == 3) {
			i = utils::random(0, size - 1);
		}
		s.push_back(i);
		while (s.size() < size) {
			int min = INT_MAX;
			int id = -1;
			for (int j = 0; j < size; j++) {
				if (i == j) continue;
				if (std::find(s.begin(), s.end(), j) != s.end())	continue;
				int distance = instance.getDistance(i, j);
				if (distance < min) {
					min = distance;
					id = j;
				}
			}
			i = id;
			s.push_back(id);
		}

		return s;
	}

}

void TabuSearch::run()
{
	int size = instance.getSize();
	vector<int> x0 = generateStartingSolution();

	vector<int> xopt = x0;
	int optimalDistance = instance.calculateCostFunction(xopt);
	

	int iteration_count = 0;

	while (iteration_count < NUMBER_OF_ITERATIONS) {
		int repeat_count = 0;
		tabuList = vector<vector<int>>(size, vector<int>(size, 0));
		int i = 0;
		while (iteration_count < NUMBER_OF_ITERATIONS && repeat_count < 10) {
			vector<int> bestNearbySolution = getBestNearbySolution(x0, i, optimalDistance);
			x0 = bestNearbySolution;
			int currentDistance = instance.calculateCostFunction(x0);
			if (currentDistance == optimalDistance) {
				repeat_count++;
			} else if (currentDistance < optimalDistance) {
				xopt = x0;
				optimalDistance = currentDistance;
				repeat_count = 0;
			}
			i++;
			iteration_count++;
		}
		x0 = generateRandomSolution();
	}

	finalSolution = xopt;
	finalDistance = optimalDistance;
}

TabuSearch::TabuSearch(Instance& instance, int neighbourhood, int number_of_iterations, int tabu_length, int starting_solution) :
	instance(instance), NEIGHBOURHOOD(neighbourhood), NUMBER_OF_ITERATIONS(number_of_iterations), TABU_LENGTH(tabu_length), STARTING_SOLUTION(starting_solution)
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

vector<int> TabuSearch::getBestNearbySolution(vector<int> solution, int it, int optimalDistance)
{
	int instanceSize = instance.getSize();
	int referenceDistance = instance.calculateCostFunction(solution);
	int bestDistance = INT_MAX;
	vector<int> candidateSolution = solution;
	vector<int> bestNearbySolution;
	int bestI = 0, bestJ = 0;
	for (int i = 0; i < instanceSize; i++) {
		for (int j = i + 1; j < instanceSize; j++) {
			neighbouringSolution(candidateSolution, i, j);
			int currentDistance = instance.calculateCostFunction(candidateSolution);
			if ((tabuList[i][j] <= it && currentDistance < bestDistance) || currentDistance < optimalDistance) { // ) {
				bestI = i;
				bestJ = j;
				bestDistance = currentDistance;
				tabuList[i][j] = it + TABU_LENGTH;
				tabuList[j][i] = it + TABU_LENGTH;
			}
			candidateSolution = solution;
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

vector<int> TabuSearch::generateStartingSolution()
{
	if (STARTING_SOLUTION == 0) {
		vector<int> s;
		for (int i = 0; i < instance.getSize(); i++) {
			s.push_back(i);
		}
		return s;
	}
	if (STARTING_SOLUTION == 1) {
		vector<int> s;
		while (s.size() < instance.getSize()) {
			int i = utils::random(0, instance.getSize() - 1);
			if (std::find(s.begin(), s.end(), i) == s.end()) {
				s.push_back(i);
			}
		}
		return s;
	}
	if (STARTING_SOLUTION == 2 || STARTING_SOLUTION == 3) {
		int size = instance.getSize();
		vector<int> s;
		int i = 0;
		if (STARTING_SOLUTION == 3) {
			i = utils::random(0, size - 1);
		}
		s.push_back(i);
		while (s.size() < size) {
			int min = INT_MAX;
			int id = -1;
			for (int j = 0; j < size; j++) {
				if (i == j) continue;
				if (std::find(s.begin(), s.end(), j) != s.end())	continue;
				int distance = instance.getDistance(i, j);
				if (distance < min) {
					min = distance;
					id = j;
				}
			}
			i = id;
			s.push_back(id);
		}

		return s;
	}

}

vector<int> TabuSearch::generateRandomSolution()
{
	vector<int> s;
	while (s.size() < instance.getSize()) {
		int i = utils::random(0, instance.getSize() - 1);
		if (std::find(s.begin(), s.end(), i) == s.end()) {
			s.push_back(i);
		}
	}
	return s;
}

void Genetic::run()
{
	population = generatePopulation();

	Tour bestTour = *min_element(population.begin(), population.end(), [](Tour a, Tour b) {
		return a.length < b.length;
		});

	for (int i = 0; i < ITERATION_COUNT; i++) {
		vector<Tour> nextGeneration;
		Tour lastBest = bestTour;
		while (nextGeneration.size() < POPULATION_SIZE) {
			if (utils::random(0.0, 1.0) < Pc) {
				Tour p1 = selection();
				Tour p2 = selection();
				vector<Tour> children = crossPair(p1, p2);

				for (int j = 0; j < children.size(); j++) {
					if (utils::random(0.0, 1.0) < Pm) {
						mutate(children[j]);
					}
					if (children[j].length < bestTour.length) {
						bestTour = children[j];
					}
				}

				Tour bestParent = min({ p1, p2 }, [](Tour a, Tour b) {
					return a.length < b.length;
					});

				Tour bestChild = min({ children[0], children[1] }, [](Tour a, Tour b) {
					return a.length < b.length;
					});

				if (bestParent.length < bestChild.length) {
					nextGeneration.push_back(bestParent);
					nextGeneration.push_back(bestChild);
				}
				else {
					nextGeneration.push_back(children[0]);
					nextGeneration.push_back(children[1]);
				}
			}
			
		}

		population = nextGeneration;
	}
	
	finalSolution = bestTour.cities;
	finalDistance = bestTour.length;
}

int Genetic::getFinalDistance()
{
	return finalDistance;
}

Genetic::Genetic(InstanceVector& instance, int ITERATION_COUNT, int POPULATION_SIZE, double Pc, double Pm, int SELECTION_METHOD, int CROSS_OPERATOR, int MUTATION_OPERATOR) :
	instance(instance), ITERATION_COUNT(ITERATION_COUNT), POPULATION_SIZE(POPULATION_SIZE), Pc(Pc), Pm(Pm), SELECTION_METHOD(SELECTION_METHOD), CROSS_OPERATOR(CROSS_OPERATOR), MUTATION_OPERATOR(MUTATION_OPERATOR)
{
}

vector<Tour> Genetic::generatePopulation()
{
	vector<Tour> pop;
	while (pop.size() < POPULATION_SIZE) {
		vector<int> solution;
		int instance_size = instance.getSize();
		while (solution.size() < instance_size) {
			int i = utils::random(0, instance_size - 1);
			if (find(solution.begin(), solution.end(), i) == solution.end()) {
				solution.push_back(i);
			}
		}
		Tour t(solution, instance.calculateCostFunction(solution));
		pop.push_back(t);
	}
	return pop;
}

Tour Genetic::selection()
{
	switch (SELECTION_METHOD) {
	case 0:
		return _selection_roulette();
	case 1:
		return _selection_tournament();
	}
}

Tour Genetic::_selection_roulette()
{
	double sum_of_fittness = 0;
	for (auto t : population) {
		sum_of_fittness += t.fitness;
	}
	double p = utils::random(0.0, 1.0);
	double p_x = 0.0;
	for (auto t : population) {
		p_x += t.fitness / sum_of_fittness;
		if (p < p_x) {
			return t;
		}
	}
}

Tour Genetic::_selection_tournament()
{
	Tour* best = &population[utils::random(0, POPULATION_SIZE - 1)];
	for (int i = 1; i < TOURNAMENT_SIZE; i++) {
		int id = utils::random(0, POPULATION_SIZE - 1);
		Tour* candidate = &population[id];
		if (candidate->length < best->length) {
			best = candidate;
		}
	}
	return *best;
}

vector<Tour> Genetic::crossPair(Tour p, Tour q)
{
	switch (CROSS_OPERATOR) {
	case 0:
		return _cross_PMX(p, q);
	case 1:
		return _cross_OX(p, q);
	}
}

vector<Tour> Genetic::_cross_PMX(Tour p, Tour q)
{
	struct mapping {
		int from;
		int to;
	};

	vector<mapping> mapping1, mapping2;

	int i = utils::random(0, instance.getSize() - 2);
	int j = utils::random(0, instance.getSize() - 1);
	while (j == i) {
		j = utils::random(0, instance.getSize() - 1);
	}

	int begin = min({ i, j });
	int end = max({ i, j });

	vector<int> child1, child2;
	child1.resize(instance.getSize(), -1);
	child2.resize(instance.getSize(), -1);

	for (int i = begin; i < end; i++) {
		mapping m1, m2;
		m1.from = p.cities[i];
		m1.to = q.cities[i];
		m2.to = p.cities[i];
		m2.from = q.cities[i];
		mapping1.push_back(m1);
		mapping2.push_back(m2);
		child1[i] = m1.to;
		child2[i] = m1.from;
	}
	
	for (int i = 0, j = 0; i < instance.getSize(); i++, j++) {
		if (i < begin || i >= end) {
			int p1_city = p.cities[i];
			int p2_city = q.cities[i];
			while (find(child1.begin() + begin, child1.begin() + end, p1_city) != child1.begin() + end) {
				for (auto mapping : mapping2) {
					if (mapping.from == p1_city) {
						p1_city = mapping.to;
						break;
					}
				}
			}
			while (find(child2.begin() + begin, child2.begin() + end, p2_city) != child2.begin() + end) {
				for (auto mapping : mapping1) {
					if (mapping.from == p2_city) {
						p2_city = mapping.to;
						break;
					}
				}
			}
			child1[i] = p1_city;
			child2[i] = p2_city;
		}
	}

	Tour c1(child1, instance.calculateCostFunction(child1));
	Tour c2(child2, instance.calculateCostFunction(child2));

	return vector<Tour>({ c1, c2 });
}

vector<Tour> Genetic::_cross_OX(Tour p, Tour q)
{
	vector<int> parent1 = p.cities;
	vector<int> parent2 = q.cities;
	int size = parent1.size();

	int n1 = utils::random(0, size);
	int n2 = utils::random(0, size - 1);

	int start = min({ n1, n2 });
	int end = max({ n1, n2 });

	vector<int> child1;
	vector<int> child2;

	for (int i = start; i < end; i++) {
		child1.push_back(parent1[i]);
		child2.push_back(parent2[i]);
	}

	int geneIndex = 0;
	int geneInParent1 = 0;
	int geneInParent2 = 0;

	for (int i = 0; i < size; i++) {
		geneIndex = (end + i) % size;
		geneInParent1 = parent1[geneIndex];
		geneInParent2 = parent2[geneIndex];

		if (find(child1.begin(), child1.end(), geneInParent2) == child1.end()) {
			child1.push_back(geneInParent2);
		}

		if (find(child2.begin(), child2.end(), geneInParent1) == child2.end()) {
			child2.push_back(geneInParent1);
		}
	}

	rotate(child1.begin(), child1.begin() + start, child1.end());
	rotate(child2.begin(), child2.begin() + start, child2.end());

	Tour c1(child1, instance.calculateCostFunction(child1));
	Tour c2(child2, instance.calculateCostFunction(child2));

	return vector<Tour>({ c1, c2 });
}

void Genetic::mutate(Tour& x)
{
	switch (MUTATION_OPERATOR) {
	case 0:
		_mutate_inversion(x);
		break;
	case 1:
		_mutate_insertion(x);
	case 2:
		_mutate_transposition(x);
	case 3:
		_mutate_scramble(x);
	}

	x.length = instance.calculateCostFunction(x.cities);
}

void Genetic::_mutate_inversion(Tour& x)
{
	int i = utils::random(0, instance.getSize() - 2);
	int j = utils::random(0, instance.getSize() - 1);

	int start = min({ i, j });
	int stop = max({ i, j });

	reverse(x.cities.begin() + start, x.cities.begin() + stop);
}

void Genetic::_mutate_insertion(Tour& x)
{
	int i = utils::random(0, instance.getSize() - 2);
	int j = utils::random(0, instance.getSize() - 1);

	int begin = min({ i, j });
	int end = max({ i, j });

	int c = x.cities[end];

	x.cities.erase(x.cities.begin() + end);
	x.cities.insert(x.cities.begin() + begin, c);
}

void Genetic::_mutate_scramble(Tour& x)
{
	int i = utils::random(0, instance.getSize() - 2);
	int j = utils::random(0, instance.getSize() - 1);

	int begin = min({ i, j });
	int end = max({ i, j });

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

	shuffle(x.cities.begin() + begin, x.cities.begin() + end, default_random_engine(seed));
}

void Genetic::_mutate_transposition(Tour& x)
{
	int i = utils::random(0, instance.getSize() - 1);
	int j;
	do {
		j = utils::random(0, instance.getSize() - 1);
	} while (i == j);

	iter_swap(x.cities.begin() + i, x.cities.begin() + j);
}

void AntColony::run()
{
	vector<int> assignedCities;
	vector<Ant> ants;
	feromone = vector<vector<double>>(instance.getSize(), vector<double>(instance.getSize(), initial_feromone));
	while(ants.size() < ANTS_COUNT) {
		int city;
		do {
			city = utils::random(0, instance.getSize() - 1);
		} while (find(assignedCities.begin(), assignedCities.end(), city) != assignedCities.end());
		
		assignedCities.push_back(city);
		ants.push_back(Ant(city, instance.getSize()));
	}

	vector<int> bestPath;
	int bestLength = INT_MAX;
	chrono::high_resolution_clock::time_point startTime, endTime;
	startTime = chrono::high_resolution_clock::now();
	while(chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - startTime).count() < MAX_TIME) {
		bool uni_path = true;
		vector<int> lastPath;
		for (int i = 0; i < ants.size(); i++) {
			vector<int> path({ ants[i].starting_city });
			while (path.size() < instance.getSize()) {
				int next_city = selectNextCity(ants[i], path[path.size() - 1]);
				path.push_back(next_city);
			}
			int length = instance.calculateCostFunction(path);
			if (STRATEGY == 2) {
				double delta = feromone_quantity / (double)length;
				for (int a = 0; a < instance.getSize(); a++) {
					for (int b = 0; b < instance.getSize(); b++) {
						feromone[a][b] *= vaporating_factor;
					}
				}
		
				int previous_city = path[0];
				for (int a = 1; a < path.size(); a++) {
					int next_city = path[a];
					feromone[previous_city][next_city] += delta;
					previous_city = next_city;
				}
			}
			if (length < bestLength) {
				bestPath = path;
				bestLength = length;
			}

			if (i == 0) {
				lastPath = path;
			}
			else {
				if (path != lastPath) {
					uni_path = false;
				}
			}
			ants[i].resetTabu();
			if (chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - startTime).count() >= MAX_TIME) {
				break;
			}
		}
		if (uni_path) {
			cout << "UNI_PATH" << endl;
			break;
		}
	}

	finalSolution = bestPath;
	finalDistance = bestLength;
}

int AntColony::selectNextCity(Ant &ant, int city)
{
	double sum_of_possible = 0.0;
	for (int c = 0; c < instance.getSize(); c++) {
		if (ant.getTabu(c) == 1) continue;
		double p = parametersCalculate(city, c);
		sum_of_possible += p;
	}
	double biggestProbability = -1.0;
	int mostProbable = -1;
	// TODO: zamienic i zrobic jak w genetycznym. Bez tego nie ma sensu wiecej niz jedna iteracja
	for (int c = 0; c < instance.getSize(); c++) {
		if (ant.getTabu(c) == 1) continue;
		double p = parametersCalculate(city, c) / sum_of_possible;
		if (p >= biggestProbability) {
			mostProbable = c;
			biggestProbability = p;
		}
	}

	if (STRATEGY != 2) {
		for (int i = 0; i < instance.getSize(); i++) {
			for (int j = 0; j < instance.getSize(); j++) {
				feromone[i][j] = feromone[i][j] * vaporating_factor;
			}
		}
		feromone[city][mostProbable] += delta_feromone(city, mostProbable);
	}
	
	ant.setTabu(mostProbable);

	return mostProbable;
}

double AntColony::delta_feromone(int i, int j)
{
	switch (STRATEGY) {
	case 0:
		return feromone_quantity;
	case 1:
		return instance.getDistance(i, j) == 0 ? feromone_quantity : feromone_quantity / (double)instance.getDistance(i, j);
	}
}

vector<int> AntColony::getFinalSolution()
{
	return finalSolution;
}

int AntColony::getFinalDistance()
{
	return finalDistance;
}

double AntColony::parametersCalculate(int i, int j)
{
	int distance = instance.getDistance(i, j);
	double tau, eta;
	if (distance == 0) {
		eta = pow(2.0, beta);
	}
	else {
		eta = pow(1.0 / (double)distance, beta);
	}

	double fero = max({ feromone[i][j], 0.000001 });
	tau = pow(fero, alpha);

	return tau * eta;
}

Ant::Ant(int starting_city, int size)
{
	this->starting_city = starting_city;
	tabuList.resize(size);
	this->size = size;
	resetTabu();
}

void Ant::resetTabu()
{
	tabuList = vector<bool>(size, false);
	tabuList[starting_city] = 1;
}

int Ant::getTabu(int pos)
{
	return tabuList[pos];
}

void Ant::setTabu(int pos)
{
	tabuList[pos] = 1;
}
