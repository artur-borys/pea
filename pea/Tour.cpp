#include "pch.h"
#include "Tour.h"

Tour::Tour(vector<int> cities, int length) : cities(cities), length(length)
{
	fitness = 1 / (double)length;
}
