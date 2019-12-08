#include "pch.h"
#include "measure.h"

void measure::start()
{
	startTime = chrono::high_resolution_clock::now();
}

void measure::stop()
{
	endTime = chrono::high_resolution_clock::now();
}

void measure::printResult()
{
	long ms = chrono::duration_cast<chrono::milliseconds>(measure::endTime - measure::startTime).count();
	if (ms > 0) {
		cout << ms << " ms";
	}
	else {
		cout << chrono::duration_cast<chrono::microseconds>(measure::endTime - measure::startTime).count() << " mikrosekund";
	}
}

long measure::result()
{
	long ms = chrono::duration_cast<chrono::milliseconds>(measure::endTime - measure::startTime).count();
	if (1) {
		unit = "ms";
		return ms;
	}
	else {
		unit = "us";
		return chrono::duration_cast<chrono::microseconds>(measure::endTime - measure::startTime).count();
	}
}