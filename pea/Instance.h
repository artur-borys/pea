#pragma once
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

class InstanceError : exception {

};

class Instance
{
public:
	Instance();
	Instance(string name, size_t size, int **data, bool debugging = false);
	static Instance createFromFile(string path);
	int calculateCostFunction(int *points = NULL);
	int calculateCostFunction(vector<int> points);
	int calculateCostFunction();
	void print();
	bool setDebugging(bool debugging);
	int** getData();
	size_t getSize();
	int getDistance(int i, int j);
	~Instance();
private:
	static Instance readFromFile(string path);
	string name;
	size_t size;
	int **data;
	bool debugging;
};

