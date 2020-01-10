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
#include <filesystem>
#include <iomanip>

#define INSERT 0
#define SWAP 1
#define INVERT 2

#define NATURAL_SOLUTION 0
#define RANDOM_SOLUTION 1
#define GREEDY_0 2
#define GREEDY_RANDOM 3

namespace fs = std::filesystem;

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

string getNM(int n) {
	switch (n) {
	case 0:
		return "INSERT";
	case 1:
		return "SWAP";
	case 2:
		return "INVERT";
	}
}

string getSS(int ss) {
	switch (ss) {
	case 0:
		return "NATURAL_SOLUTION";
	case 1:
		return "RANDOM_SOLUTION";
	case 2:
		return "GREEDY_0";
	case 3:
		return "GREEDY_RANDOM";
	}
}

struct results {
	double avgTime;
	double bestScore;
	double avgScore;
	string unit;
	string name;
};

results testSA(string path, int optimal, double INIT_TEMP, double MIN_TEMP, double TEMP_FACTOR, int nType, int sSolution, int repeats) {
	Instance instance = Instance::createFromFile(path);
	measure t;
	results r;
	r.bestScore = INT_MAX;
	r.name = path + "_" + to_string(INIT_TEMP) + "_" + to_string(MIN_TEMP) + "_" + to_string(TEMP_FACTOR) +
		"_" + getNM(nType) + "_" + getSS(sSolution);
	for (int i = 0; i < repeats; i++) {
		{
			SimulatedAnnealing sa(instance, INIT_TEMP, MIN_TEMP, TEMP_FACTOR, nType, sSolution);
			t.start();
			sa.run();
			t.stop();
			r.avgTime += t.result();
			double score = (double)abs(optimal - sa.getFinalDistance()) / (double)optimal * 100.0;
			if (score < r.bestScore) {
				r.bestScore = score;
			}
			r.avgScore += score;
		}
	}

	r.avgTime /= (double)repeats;
	r.avgScore /= (double)repeats;
	r.unit = t.unit;
	return r;
}

results testTS(string path, int optimal, int ITERATIONS, int tabu_length, int nType, int sSolution, int repeats) {
	Instance instance = Instance::createFromFile(path);
	measure t;
	results r;
	r.bestScore = INT_MAX;
	r.name = path + "_" + to_string(ITERATIONS) + "_" + to_string(tabu_length) + 
		"_" + getNM(nType) + "_" + getSS(sSolution);
	for (int i = 0; i < repeats; i++) {
		{
			TabuSearch sa(instance, nType, ITERATIONS, tabu_length, sSolution);
			t.start();
			sa.run();
			t.stop();
			r.avgTime += t.result();
			double score = (double)abs(optimal - sa.getFinalDistance()) / (double)optimal * 100.0;
			if (score < r.bestScore) {
				r.bestScore = score;
			}
			r.avgScore += score;
		}
	}

	r.avgTime /= (double)repeats;
	r.avgScore /= (double)repeats;
	r.unit = t.unit;
	return r;
}

results testSA(Instance instance, int optimal, double INIT_TEMP, double MIN_TEMP, double TEMP_FACTOR, int nType, int sSolution, int repeats) {
	measure t;
	results r;
	r.bestScore = INT_MAX;
	r.name = instance.name + "_" + to_string(INIT_TEMP) + "_" + to_string(MIN_TEMP) + "_" + to_string(TEMP_FACTOR) +
		"_" + getNM(nType) + "_" + getSS(sSolution);
	for (int i = 0; i < repeats; i++) {
		{
			SimulatedAnnealing sa(instance, INIT_TEMP, MIN_TEMP, TEMP_FACTOR, nType, sSolution);
			t.start();
			sa.run();
			t.stop();
			r.avgTime += t.result();
			double score = (double)abs(optimal - sa.getFinalDistance()) / (double)optimal;
			if (score < r.bestScore) {
				r.bestScore = score;
			}
			r.avgScore += score;
		}
	}

	r.avgTime /= (double)repeats;
	r.avgScore /= (double)repeats;
	r.unit = t.unit;
	return r;
}

results testTS(Instance instance, int optimal, int ITERATIONS, int tabu_length, int nType, int sSolution, int repeats) {
	measure t;
	results r;
	r.bestScore = INT_MAX;
	r.name = instance.name + "_" + to_string(ITERATIONS) + "_" + to_string(tabu_length) +
		"_" + getNM(nType) + "_" + getSS(sSolution);
	for (int i = 0; i < repeats; i++) {
		{
			TabuSearch sa(instance, nType, ITERATIONS, tabu_length, sSolution);
			t.start();
			sa.run();
			t.stop();
			r.avgTime += t.result();
			double score = (double)abs(optimal - sa.getFinalDistance()) / (double)optimal;
			if (score < r.bestScore) {
				r.bestScore = score;
			}
			r.avgScore += score;
		}
	}

	r.avgTime /= (double)repeats;
	r.avgScore /= (double)repeats;
	r.unit = t.unit;
	return r;
}

void measureSA() {
	int repeats = 50;
	double temp_factor[] = { 0.995, 0.999, 0.9995 };
	int n[] = { INSERT, SWAP, INVERT };
	int s[] = { NATURAL_SOLUTION, RANDOM_SOLUTION, GREEDY_RANDOM };
	vector<int> optimal = {
		/*212,
		202,
		264,
		269,
		125,
		291,
		156,*/
		2085,
		/*187,
		2707,
		1272,*/
		1530,
		/*1613,
		14422,*/
		25395,
		//1950,
		6942
	};
	vector<string> paths = {
		/*"data10.txt",
		"data11.txt",
		"data12.txt",
		"data13.txt",
		"data14.txt",
		"data15.txt",
		"data16.txt",*/
		"data17.txt",
		/*"data18.txt",
		"data21.txt",
		"data24.txt",*/
		"data39.txt",
		//"data45.txt",
		/*"data48.txt",*/
		"data58.txt",
		//"data71.txt",
		"data120.txt"
	};
	double min_temp = 0.995;
	double init_temp = 10.0;
	int a = 0;
	for (int a = 0; a < optimal.size(); a++) {
		Instance instance = Instance::createFromFile("instances_local/" + paths[a]);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					ofstream f("measurements/sa.csv", ios::out | ios::app);
					results r;
					r = testSA(instance, optimal[a], init_temp, min_temp, temp_factor[i], n[j], s[k], repeats);
					f << instance.getSize() << ";" << temp_factor[i] << ";"
						<< getNM(n[j]) << ";" << getSS(s[k]) << ";" << r.avgTime << ";" << r.unit << ";"
						<< r.avgScore << ";" << r.bestScore << ";" << endl;
					f.close();
				}
			}
		}
	}
	
}

void measureTS() {
	int repeats = 50;
	int ITERATIONS[] = { 100, 1000, 5000 };
	int length_div[] = { 2, 1 };
	int n[] = { INSERT, SWAP, INVERT };
	int s[] = { NATURAL_SOLUTION, RANDOM_SOLUTION, GREEDY_RANDOM };
	vector<int> optimal = {
		/*212,
		202,
		264,
		269,
		125,
		291,
		156,*/
		2085,
		/*187,
		2707,
		1272,*/
		1530,
		/*1613,
		14422,*/
		25395,
		//1950,
		6942
	};
	vector<string> paths = {
		/*"data10.txt",
		"data11.txt",
		"data12.txt",
		"data13.txt",
		"data14.txt",
		"data15.txt",
		"data16.txt",*/
		"data17.txt",
		/*"data18.txt",
		"data21.txt",
		"data24.txt",*/
		"data39.txt",
		//"data45.txt",
		/*"data48.txt",*/
		"data58.txt",
		//"data71.txt",
		"data120.txt"
	};
	int a = 0;
	for (int a = 0; a < optimal.size(); a++) {
		Instance instance = Instance::createFromFile("instances_local/" + paths[a]);
		int size = instance.getSize();
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					for (int l = 0; l < 2; l++) {
						ofstream f("measurements/ts.csv", ios::out | ios::app);
						results r;
						r = testTS(instance, optimal[a], ITERATIONS[i], size/length_div[l], n[j], s[k], repeats);
						f << size << ";" << ITERATIONS[i] << ";" << size / length_div[l] << ";" 
							<< getNM(n[j]) << ";" << getSS(s[k]) << ";" << r.avgTime << ";" << r.unit << ";"
							<< r.avgScore << ";" << r.bestScore << ";" << endl;
						f.close();
					}
				}
			}
		}
	}

}

void measureLocal() {
	thread sa(measureSA);
	thread ts(measureTS);
	sa.join();
	ts.join();
}

void measureSASmall() {
	int repeats = 50;
	vector<int> optimal = {
		212,
		202,
		264,
		269,
		125,
		291,
		156,
		2085,
		187,
		2707,
		1286,
		1530,
		1613,
		14422,
		25395,
		1950,
		6942
	};
	vector<string> paths = {
		"data10.txt",
		"data11.txt",
		"data12.txt",
		"data13.txt",
		"data14.txt",
		"data15.txt",
		"data16.txt",
		"data17.txt",
		"data18.txt",
		"data21.txt",
		"data34.txt",
		"data39.txt",
		"data45.txt",
		"data48.txt",
		"data58.txt",
		"data71.txt",
		"data120.txt"
	};
	double min_temp = 0.995;
	double init_temp = 100.0;
	int a = 0;
	for (int a = 0; a < optimal.size(); a++) {
		Instance instance = Instance::createFromFile("instances_local/" + paths[a]);
		int size = instance.getSize();
		ofstream f("measurements/sa_small_avg.csv", ios::out | ios::app);
		ofstream f1("measurements/sa_small_best.csv", ios::out | ios::app);
		ofstream f2("measurements/sa_small_time.csv", ios::out | ios::app);
		results r;
		vector<double> temp_factor = { 0.99, 0.995, 0.999, 0.9995 };
		for (int i = 0; i < temp_factor.size(); i++) {
			r = testSA(instance, optimal[a], init_temp, min_temp, temp_factor[i], INVERT, GREEDY_RANDOM, repeats);
			f << std::fixed << std::setprecision(4) << r.avgScore << " ";
			f1 << std::fixed << std::setprecision(4) << r.bestScore << " ";
			f2 << std::fixed << std::setprecision(4) << r.avgTime << " ";
		}

		f << endl;
		f1 << endl;
		f2 << endl;
	}

}

void measureTSSmall() {
	int repeats = 50;
	int ITERATIONS[] = { 100, 1000, 5000 };
	int length_div[] = { 2, 1 };
	int n[] = { INSERT, SWAP, INVERT };
	int s[] = { NATURAL_SOLUTION, RANDOM_SOLUTION, GREEDY_RANDOM };
	vector<int> optimal = {
		212,
		202,
		264,
		269,
		125,
		291,
		156,
		2085,
		187,
		2707,
		1286,
		1530,
		1613,
		14422,
		25395,
		1950,
		6942
	};
	vector<string> paths = {
		"data10.txt",
		"data11.txt",
		"data12.txt",
		"data13.txt",
		"data14.txt",
		"data15.txt",
		"data16.txt",
		"data17.txt",
		"data18.txt",
		"data21.txt",
		"data34.txt",
		"data39.txt",
		"data45.txt",
		"data48.txt",
		"data58.txt",
		"data71.txt",
		"data120.txt"
	};
	int a = 0;
	for (int a = 0; a < optimal.size(); a++) {
		Instance instance = Instance::createFromFile("instances_local/" + paths[a]);
		int size = instance.getSize();
		ofstream f("measurements/ts_small_avg.csv", ios::out | ios::app);
		ofstream f1("measurements/ts_small_best.csv", ios::out | ios::app);
		ofstream f2("measurements/ts_small_time.csv", ios::out | ios::app);
		results r;
		vector<int> iters = { 100, 500, 1000, 2000, 3000, 4000 };
		for (int i = 0; i < iters.size(); i++) {
			r = testTS(instance, optimal[a], iters[i], size / 2, INVERT, GREEDY_RANDOM, repeats);
			f << std::fixed << std::setprecision(4) << r.avgScore  << " ";
			f1 << std::fixed << std::setprecision(4) << r.bestScore << " ";
			f2 << std::fixed << std::setprecision(4) << r.avgTime << " ";
		}

		f << endl;
		f1 << endl;
		f2 << endl;
	}

}

int main()
{
	int selection;
	string warning;
	thread t1(measureTSSmall);
	thread t2(measureSASmall);
	t1.join();
	t2.join();
	return 0;
	/*{
		SimulatedAnnealing sa(inst, 5999.0, 0.995, 0.005, INVERT, GREEDY_0);
		t.start();
		sa.run();
		t.stop();
		cout << sa.getFinalDistance() << endl;
		cout << abs(optimal - sa.getFinalDistance()) * 100.0 / optimal << "% worse than optimal" << endl;
		t.printResult();
		cout << "\n\n";
	}*/
	//cout << "TS: 1000, instanze size, SWAP, GREEDY_0" << endl;
	/*{
		TabuSearch ts(inst, INVERT, 5000, inst.getSize(), GREEDY_RANDOM);
		t.start();
		ts.run();
		t.stop();
		cout << ts.getFinalDistance() << endl;
		cout << abs(optimal - ts.getFinalDistance()) * 100.0 / optimal << "% worse than optimal" << endl;
		t.printResult();
		cout << "\n\n";
	}*/
	//presentation();
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
