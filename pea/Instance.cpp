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

int Instance::calculateCostFunction(int *points)
{
	int total = 0;
	int previousIndex = points[0];
	for (int i = 1; i < size; i++) {
		int nextIndex = points[i];
		int distance = data[previousIndex][nextIndex];
		total += distance;
		previousIndex = nextIndex;
	}
	total += data[previousIndex][points[0]];
	return total;
}

int Instance::calculateCostFunction(vector<int> points) {
	int total = 0;
	int previousIndex = points[0];
	for (int i = 1; i < size; i++) {
		int nextIndex = points[i];
		int distance = data[previousIndex][nextIndex];
		total += distance;
		previousIndex = nextIndex;
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

int Instance::getDistance(int i, int j)
{
	return data[i][j];
}


Instance::~Instance()
{

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
				int num;
				lineStream >> num;
				if (num == 0) {
					num = INT_MAX;
				}
				data[rowNum][i] = num;
			}

			rowNum++;
		}

		file.close();

		return Instance(instanceName, instanceSize, data);
	}
	return Instance();
}

InstanceVector::InstanceVector()
{
}

InstanceVector::InstanceVector(string name, size_t size, vector<vector<int>> data) : name(name), size(size), data(data)
{
}

InstanceVector InstanceVector::createFromFile(string path)
{
	string name;
	size_t size;
	vector<vector<int>> data;
	ifstream f(path);

	if (!f.is_open()) {
		cout << "Error when reading file" << endl;
	}
	else {
		string line;

		f >> name;
		f >> size;

		while (getline(f, line)) {
			stringstream ss(line);
			if (line.empty()) {
				continue;
			}

			int c;

			vector<int> row;
			for (int i = 0; i < size; i++) {
				ss >> c;
				//if(!s.empty())
				row.push_back(c);
			}

			data.push_back(row);

		}

		return InstanceVector(name, size, data);
	}
}

int InstanceVector::calculateCostFunction(vector<int> points)
{
	int total = 0;
	int previousIndex = points[0];
	for (int i = 1; i < size; i++) {
		int nextIndex = points[i];
		int distance = data[previousIndex][nextIndex];
		total += distance;
		previousIndex = nextIndex;
	}
	total += data[previousIndex][points[0]];
	return total;
}

vector<vector<int>> InstanceVector::getData()
{
	return data;
}

size_t InstanceVector::getSize()
{
	return size;
}

int InstanceVector::getDistance(int i, int j)
{
	return data[i][j];
}

InstanceVector::~InstanceVector()
{
	data.clear();
}
