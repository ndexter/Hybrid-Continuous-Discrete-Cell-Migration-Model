#ifndef DIFFUSIONMATH_H_
#define DIFFUSIONMATH_H_
#include "BioChem.h"
#include "vec.h"
#include <math.h>
#include <iostream>



double func(double y)
{
	return exp(-y*y);
}
 
//define integration using a riemann sum
double integrate(double lowB, double upB)   // lowB is the lower limit, upB is the upper limit
{ 
    double dx=0.1;          //width of square for riemann sum
    double sum=0.0;           //final value
	int start = lowB/dx;
	int stop = upB/dx;
		for (int i = start; i < stop; i++)
        {
        	sum=sum+func(i*dx)*dx;
        }
	return sum;
}
 
//define the error function
double erf(double z)
{  
	 return (2/sqrt(atan(1)))*integrate(0,z);
}

// define the equation for the chemoattractant concentration
double diffusion(double pos, double time, int num)
{ 
/*	BioChem* env = BioChem::getInstance();
	vec h= env->getDomain();      //h is chamber height,
	double D = env->getDiffCoeff();     //D is diffusion constant for chemoattractant,
	double a = env->getInitialConc();      //a is concentration of cemoattractant in lower well at start*/
	double a = 1.0E-4;
	vec h( 7.0);
	double D = 5.0E-1;
	return (a/2)*(erf((h[0]*(-1+4*num)+pos)/(2*sqrt(D*time))) + erf((h[0]*(3-4*num)-pos)/(2*sqrt(D*time))));
}


#endif /*DIFFUSIONMATH_H_*/
