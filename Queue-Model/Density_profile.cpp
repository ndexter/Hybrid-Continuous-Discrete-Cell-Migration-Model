#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <vector>
#include <iostream>
using namespace std;

void profile(int argc, const char** argv)
{
	vector<double> x;
	double x_min = DBL_MAX, x_max = 0.0;
	if (argc != 3)
	{
		cout <<	 argv[0] << " <dx> <file>" << endl;
	}
	double dx = atof(argv[1]);
	ifstream fin(argv[2]);
	while (!fin.eof())
	{
		double value;
		fin >> value;
		if (!fin.eof())
		{
			x.push_back(value);
			x_min = min(x_min,value);
			x_max = max(x_max,value);
		}
	}
	fin.close();
	int max_count = 0;
	int buckets = (int)ceil((x_max-x_min)/dx);
	int* count = new int[buckets];
	for (int i = 0; i < buckets; i++)
	{
		count[i] = 0;
	}
	for (unsigned i = 0; i < x.size(); i++)
	{
		int index = (int)((x[i]-x_min)/dx);
		count[index]++;
		max_count = max(max_count,count[index]);
	}
	for (int i = 0; i < buckets; i++)
	{
		cout << x_min+(double)i*dx << " " << 
		(double)(count[i])/max_count << endl;
	}
}


