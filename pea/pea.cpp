// pea.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include "utils.h"
#include "Instance.h"
#include <algorithm>
#include "Algorithms.h"

int main()
{
	Instance inst = Instance::createFromFile("E:\\libs\\Dokumenty\\studia\\PEA\\Projekt\\PEA\\SMALL\\data17.txt");
	//Instance inst = Instance::createFromFile("data4.txt");
	BruteForce bf(inst);
	DynamicProgramming dp(inst, 0);

	inst.print();

	//bf.run();
	dp.run();

	cout << "Distance: " << bf.getFinalDistance() << endl;
	cout << "Distance: " << dp.getFinalDistance() << endl;
	cout << "Solution: ";
	utils::printSolution(bf.getFinalSolution(), inst.getSize());
	cout << "Solution: ";
	utils::printSolution(dp.getFinalSolution());
	cout << endl;
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
