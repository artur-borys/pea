#pragma once
#include <vector>
#include "Instance.h"

using namespace std;

class Tour
{
public:
	Tour(vector<int> cities, int length);
	vector<int> cities;
	double fitness;
	int length;
	bool operator <(const Tour& b) const {
		return length < b.length;
	}
};

