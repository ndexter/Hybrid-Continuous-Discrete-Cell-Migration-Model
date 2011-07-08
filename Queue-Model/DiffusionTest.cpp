#include "DiffusionMath.h"
#include <iostream>
#include <fstream>
using namespace std;

void diffusionTest()
{
	cout << "t \t x_pos \t an_soln \t exp_value \t difference \n";
	int counter;
	double total_difference;
	ifstream fin("biochem/c");
	double x_pos, t, exp_value;
	while (fin >> t >> x_pos >> exp_value)
	{
		cout << x_pos << "\t" << t << "\t" ;
		double an_soln= 0;
		int n;
		for (n = -100; n<=100; n++)
		{
			an_soln += diffusion(x_pos, t ,n);
		}
		double difference = abs(an_soln-exp_value); 
		total_difference += difference;
		cout << an_soln << "\t" << exp_value << "\t" << difference << endl;
		counter ++;
	}
	fin.close();
	double average_difference = total_difference/counter;
	cout << "total difference = " << total_difference << endl;
	cout << "average distance = " << average_difference << endl;
	
}

	
