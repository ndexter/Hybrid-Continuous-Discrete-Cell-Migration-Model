#include <iostream>
#include "ThreeDHeat.h"
#include <sys/time.h>
#include <cfloat>

int main (int argc, char * const argv[]) {
	
	using namespace std;
	
	struct timeval tv;
	
	time_t begin, end, cputime;
	
	gettimeofday(&tv, NULL);
	
	begin=tv.tv_usec;
	
	//	cout << DBL_DIG << endl;
	//	cout << DBL_MANT_DIG << endl;
	//	cout << DBL_MAX_EXP << endl;
	//	cout << DBL_MAX << endl;
	//	cout << DBL_MIN_EXP << endl;
	//	cout << DBL_MIN << endl;
	
	double dt		= (double) 1.0E-2;							// seconds (was 1.0E-3)
	double dx   	= (double) 1.0;								// microns
	double W 		= (double) 10.0;							// microns
	double D 		= (double) 10.0;							// microns
	double H 		= (double) 20.0;							// microns
	double DC       = (double) 13.0 + (double)1.0/(double)3.0;	// microns^2/second
	double flux0    = (double) 0.0;								// flux at 0
	double fluxL    = (double) 0.0;								// flux at L
	double k        = (double) 1.0;								// 1/D
	double theta    = (double) 0.5;								// combined method coefficient
	
	ThreeDHeat *threedheat = new ThreeDHeat(dt, dx, W, D, H, DC, flux0, fluxL, k, theta);
	
	while (threedheat->getTime() < 40.0){
		threedheat->advanceTime(dt);
	}
	
	gettimeofday(&tv, NULL);
	
	end=tv.tv_usec;
	
	cputime = end - begin;
	
	cout << "begin of simulation: " << begin 
	<< "\nend of simulation: " << end 
	<< "\ncputime for simulation: " << cputime << endl;
	
	return 0;
}
