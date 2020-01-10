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
	static Instance readFromFile(string path);
	string name;
	~Instance();
private:
	size_t size;
	int **data;
	bool debugging;
};

class InstanceVector {
public:
	InstanceVector();
	InstanceVector(string name, size_t size, vector<vector<int>> data);
	static InstanceVector createFromFile(string path);
	int calculateCostFunction(vector<int> points);
	void print();
	vector<vector<int>> getData();
	size_t getSize();
	int getDistance(int i, int j);
	string name;
	~InstanceVector();
private:
	size_t size;
	vector<vector<int>> data;
};
