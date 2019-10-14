// pea.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include "utils.h"
#include "Instance.h"

int main()
{
	int size = 10;
	Instance instance = Instance::createFromFile("data10.txt");
	instance.setDebugging(true);
	instance.print();
	int *points = utils::ascendingSolution(size);
	cout << "Funkcja celu (permutacja naturalna): " << instance.calculateCostFunction(points) << endl;
	delete points;
	points = utils::randomSolution(size);
	cout << "Funkcja celu (permutacja losowa): " << instance.calculateCostFunction(points) << endl;
	delete points;
	int points2[11] = { 0, 5, 3, 2, 4, 1, 6, 9, 8, 7 };
	cout << "Funkcja celu (permutacja zadana): " << instance.calculateCostFunction(points2) << endl;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
