#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <chrono>
#include <algorithm>

using namespace std;

namespace utils{
	int *randomSolution(int size);
	int *ascendingSolution(int size);
	bool nextPermutation(int* first, int* last);
	void printSolution(int* solution, int size);
	void printSolution(vector<int> solution);
	int random(int min, int max);
	double random(double min, double max);
}
