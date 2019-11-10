#include "pch.h"
#include "Instance.h"


Instance::Instance()
{
}

Instance::Instance(string name, size_t size, int ** data, bool debugging) : name(name), size(size), data(data), debugging(debugging)
{
}

Instance Instance::createFromFile(string path)
{
	return readFromFile(path);
}

int Instance::calculateCostFunction(int * points)
{
	int total = 0;
	int previousIndex = points[0];
	for (int i = 1; i < size; i++) {
		int nextIndex = points[i];
		int distance = data[previousIndex][points[i]];
		if (debugging) {
			cout << previousIndex << "-" << points[i] << ": " << distance << endl;
		}
		total += distance;
		previousIndex = nextIndex;
	}
	if (debugging) {
		cout << previousIndex << "-" << points[0] << ": " << data[previousIndex][points[0]] << endl;
	}
	total += data[previousIndex][points[0]];
	return total;
}

int Instance::calculateCostFunction()
{
	int i = 0, j = 1, sum = 0;
	while (j < size) {
		sum += data[i][j];
		i++, j++;
	}
	sum += data[size - 1][0];
	return sum;
}

void Instance::print()
{
	cout << "Instancja" << endl;
	cout << "Nazwa: " << name << endl;
	cout << "Rozmiar: " << size << endl;
	cout << "Dane: " << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << setw(4) << data[i][j];
		}
		cout << endl;
	}
}

bool Instance::setDebugging(bool debugging)
{
	return this->debugging = debugging;
}

int** Instance::getData()
{
	return data;
}

size_t Instance::getSize()
{
	return size;
}


Instance::~Instance()
{
	delete[] data;
}

Instance Instance::readFromFile(string path)
{
	ifstream file(path);

	if (!file.is_open()) {
		cout << "Nie mozna otworzyc pliku" << endl;
	}
	else {
		string instanceName, line;
		size_t instanceSize;
		int **data, rowNum = 0;

		file >> instanceName;
		file >> instanceSize;

		data = new int*[instanceSize];

		while (getline(file, line)) {

			// czasami linijka zaczyna siê od bia³ego znaku, wiêc nale¿y go pomin¹æ
			if (line.empty()) {
				continue;
			}

			stringstream lineStream = stringstream(line);
			
			data[rowNum] = new int[instanceSize];

			for (int i = 0; i < instanceSize; i++) {
				lineStream >> data[rowNum][i];
			}

			rowNum++;
		}

		file.close();

		return Instance(instanceName, instanceSize, data);
	}
	return Instance();
}
