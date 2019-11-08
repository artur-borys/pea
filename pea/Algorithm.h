#pragma once
#include "Instance.h"
#include "utils.h"

class Algorithm
{
public:
	int finalDistance;
	int* finalSolution;
	void BruteForce(Instance& instance);
	Algorithm();
};

