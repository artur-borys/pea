// pea.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include "utils.h"
#include "Instance.h"
#include <algorithm>
#include "Algorithms.h"
#include "measure.h"
#include <thread>
#include <fstream>

Instance readInstance() {
	string path;
	cout << "Podaj sciezke do pliku: ";
	cin >> path;
	return Instance::createFromFile(path);
}

void testBNB() {
	measure m;
	string path;
	cout << "Podaj sciezke do pliku: ";
	cin >> path;
	Instance inst = Instance::createFromFile(path);
	BranchNBound bnb(inst, 0);

	m.start();
	bnb.run();
	m.stop();

	cout << "Trasa: ";
	utils::printSolution(bnb.getFinalSolution());
	cout << "Dlugosc: " << bnb.getFinalDistance() << endl;
	cout << "Czas wykonywania: ";  m.printResult(); cout << endl;
	bnb.clearQueue();
}

void measureBF_t(string path, int repeats, double *dest) {
	measure t;
	Instance inst = Instance::createFromFile(path);
	double result = 0;

	for (int i = 0; i < repeats; i++) {
		BruteForce bf(inst);

		t.start();
		bf.run();
		t.stop();

		result += t.result();
	}

	*dest = result / (double)repeats;
}

void measureBF() {
	double results[8];
	measure t;
	vector<string> paths = {
		"SMALL\\data10.txt",
		"SMALL\\data11.txt",
		"SMALL\\data12.txt",
		"SMALL\\data13.txt",
		"SMALL\\data14.txt",
	};

	double t1_result, t2_result;
	thread t1(measureBF_t, paths[0], 100, &t1_result);
	thread t2(measureBF_t, paths[1], 100, &t2_result);
	thread t3(measureBF_t, paths[2], 100, &t1_result);
	thread t4(measureBF_t, paths[3], 100, &t2_result);
	t1.join();
	ofstream f("measurements\\bf.csv");
	f << paths[0] << ";" << t1_result << ";ms;\n";
	f.close();
	t2.join();
	f = ofstream("measurements\\bf.csv", (ios::out | ios::app));
	f << paths[1] << ";" << t2_result << ";ms;\n";
	f.close();
	t3.join();
	f = ofstream("measurements\\bf.csv", (ios::out | ios::app));
	f << paths[2] << ";" << t1_result << ";ms;\n";
	f.close();
	t4.join();
	f = ofstream("measurements\\bf.csv", (ios::out | ios::app));
	f << paths[3] << ";" << t2_result << ";ms;\n";
	f.close();

	thread t5(measureBF_t, paths[4], 50, &t1_result);
	t5.join();
	f = ofstream("measurements\\bf.csv", (ios::out | ios::app));
	f << paths[4] << ";" << t1_result << ";ms;\n";
	f.close();

	cout << "Brute force measurements finished" << endl;
}



void measureDP() {
	string paths[] = {
		"SMALL\\data10.txt",
		"SMALL\\data11.txt",
		"SMALL\\data12.txt",
		"SMALL\\data13.txt",
		"SMALL\\data14.txt",
		"SMALL\\data15.txt",
		"SMALL\\data16.txt",
		"SMALL\\data17.txt",
		"SMALL\\data18.txt",
		"TSP\\data21.txt"
	};
	measure t;

	ofstream f("measurements\\dp.csv");
	f.close();

	for (int i = 0; i < 10; i++) {
		double sum = 0.0;
		Instance inst = Instance::createFromFile(paths[i]);
		ofstream f("measurements\\dp.csv", ios::out | ios::app);
		if (!f.is_open()) {
			cout << "Cannot open" << endl;
			return;
		}
		for (int j = 0; j < 100; j++) {
			DynamicProgramming dp(inst, 0);
			t.start();
			dp.run();
			t.stop();
			sum += t.result();
		}
		f << paths[i] << ";" << (double)(sum / 100.0) << ";" << t.unit << ";\n";
		f.close();
	}
	cout << "Dynamic programming measurements finished" << endl;
}

double measureBNB(string path, int repeats) {
	Instance inst = Instance::createFromFile(path);
	measure t;
	double result = 0;

	for (int i = 0; i < repeats; i++) {
		BranchNBound bnb(inst, 0);

		t.start();
		bnb.run();
		t.stop();
		bnb.clearQueue();
		result += t.result();
	}

	return (double)(result / (double)repeats);
}

void measureBNB() {
	string paths[] = {
		"SMALL\\data10.txt",
		"SMALL\\data11.txt",
		"SMALL\\data12.txt",
		"SMALL\\data13.txt",
		"SMALL\\data14.txt",
		"SMALL\\data16.txt",
	};
	measure t;

	ofstream f("measurements\\bnb.csv");
	f.close();

	for (auto path : paths) {
		double result = measureBNB(path, 100);
		f = ofstream("measurements\\bnb.csv", (ios::out | ios::app));
		f << paths[2] << ";" << result << ";ms;\n";
		f.close();
	}
	cout << "Branch & Bound measurements completed" << endl;
}

void measurements() {
	/*measureDP();
	measureBF();*/
	measureBNB();
}

int main()
{
	int selection;
	string warning;
	Instance inst = Instance::createFromFile("TSP\\data26.txt");
	/*SimulatedAnnealing sa(inst, 10000, 1, 0.005, 0);
	sa.run();
	cout << sa.getFinalDistance() << endl;
	utils::printSolution(sa.getFinalSolution());*/
	TabuSearch ta(inst, 0, 10000, inst.getSize() * 3);
	ta.run();
	cout << ta.getFinalDistance() << endl;
	utils::printSolution(ta.getFinalSolution());
	return 0;
	while (true) {
		system("CLS");
		if (!warning.empty()) {
			cout << warning << endl;
		}

		cout << "Problem komiwojazera, algorytmy dokladne, Artur Borys 241323, W4 INF" << endl;
		cout << "1. Brute force" << endl;
		cout << "2. Dynamic Programming" << endl;
		cout << "3. Branch and Bound" << endl;
		cout << "4. Pomiary" << endl;
		cout << "0. Wyjscie" << endl;
		cout << "> ";
		cin >> selection;
		if (selection == 0) {
			break;
		}
		else if (selection == 1) {
			measure m;
			string path;
			cout << "Podaj sciezke do pliku: ";
			cin >> path;
			replace(path.begin(), path.end(), '\\', '/');
			Instance inst = Instance::createFromFile(path);
			BruteForce bf(inst);
			
			m.start();
			bf.run();
			m.stop();

			cout << "Trasa: ";
			utils::printSolution(bf.getFinalSolution());
			cout << "Dlugosc: " << bf.getFinalDistance() << endl;
			cout << "Czas wykonywania: ";  m.printResult(); cout << endl;

		}
		else if (selection == 2) {
			measure m;
			string path;
			cout << "Podaj sciezke do pliku: ";
			cin >> path;
			Instance inst = Instance::createFromFile(path);
			DynamicProgramming dp(inst, 0);

			m.start();
			dp.run();
			m.stop();

			cout << "Trasa: ";
			utils::printSolution(dp.getFinalSolution());
			cout << "Dlugosc: " << dp.getFinalDistance() << endl;
			cout << "Czas wykonywania: ";  m.printResult(); cout << endl;
		}
		else if (selection == 3) {
			testBNB();
		}
		else if (selection == 4) {
			measurements();
		}
		else {
			warning = "Musisz wybrac opcje z listy!";
			continue;
		}
		system("PAUSE");
	}
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
