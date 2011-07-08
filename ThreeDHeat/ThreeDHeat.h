/*
 *  ThreeDHeat.h - A combined-method heat simulation.
 *
 *  Created by Nicholas Dexter on 2/2/09.
 *  Copyright 2009. All rights reserved.
 *
 */
#ifndef _ThreeDHeat_h_
#define _ThreeDHeat_h_
#include <fstream>
#include <cmath>
#include <sys/time.h>
//#include "adevs.h"
extern "C"{
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_permutation.h>
}

class ThreeDHeat
{
	
public:
	/**
	 * Get heat concentration at point
	 */
	double getHeat(long unsigned int *l) const;

	/**
	 * Get the Chemtoatractant Gradient at vector p = [ i, j, k ]
	 */ 
//	vec getChemtacGrad(vec p) const;
	
	/**
	 * Get the ammount of Chemtoatractant at vector p = [ i, j, k ]
	 */
//	double getChemtac(vec p) const;

	/**
	 * Get the current time
	 */
	double getTime() const;
		
	/**
	 * Get the grid spacing
	 */
	double getGridSpacing() const;
	
	/**
	 * Get index based on formula i + j*D + k*W*D
	 */
	 int getIndex(int i, int j, int k) const;
	
	/**
	 * Get coordinates x, y, and z from the index
	 */ 
	int getCoordinates(int index, int *coords) const;
	
	/**
	 * Generate the matrix equations (LHS) necessary to solve the implicit solution
	 */
	void buildLHS(gsl_matrix *LHS) const;

	/**
	 * Generate the matrix equations (RHS) necessary to solve the explicit solution
	 */
	void buildRHS(gsl_matrix *RHS) const;
	
	/**
	 * Advance the time of the simulation by dt
	 */
	void advanceTime(double dt);
	
	/**
	 * Constructor for simulation
	 */
	ThreeDHeat(double dt, double dx, double W, double D, double H, double DC, double flux0, double fluxL, double k, double theta);
	
	/**
	 * Destructor for simulation
	 */
	~ThreeDHeat();
	
private:
	
	// Length of the domain
	double W,D,H;
	
	// Diffusion coefficient
	double DC;
	
	// Current time, time step
	double t, dt;
	
	// Grid step size
	double dx;
	
	// Combined method switch
	double theta;
	
	// Flux constants at boundary
	double flux0, fluxL;
	
	// Flux scalar
	double k;
	
	// Diffusion term
	double alpha;
};

#endif
