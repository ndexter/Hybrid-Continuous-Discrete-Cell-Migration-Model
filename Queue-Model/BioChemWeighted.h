#ifndef _BioChemWeighted_h_
#define _BioChemWeighted_h_
#include "BioChemInstance.h"
#include <fstream>
#include <cmath>
#include <sys/time.h>
#include "adevs.h"
extern "C"{
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_permutation.h>
}

class BioChemWeighted:
	public BioChemInstance{
	
public:
	
	/**
	 * Get the Chemtoatractant Gradient at vector p = [ i, j, k ]
	 */ 
	vec getChemtacGrad(vec p) const;
	
	/**
	 * Get the concentration of Chemoattractant at vector p = [ i, j, k ]
	 */
	double getChemtac(vec p) const;
	
	/**
	 * Get the haptoattractant gradient at position p
	 * */
	vec getHaptacGrad(vec p) const;

	/**
	 * Get the haptoattractant concentration at position p
	 * */
	double getHaptac(vec p) const;
		
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
	 * Get the Diffusion Coefficient
	 */
	double getDiffCoeff() const { return D; }
	
	/**
	 * get the initial concentration of substance
	 */
	double getInitialConc() const { return a; }
	
	/**
	 * Get the extent of the solution domain in each dimension
	 */
	vec getDomain() const;
	/**
	 * Get the position of the bottom of the filter
	 */
	double getFilterBottomPos() const;
	/**
	 * Get the position of the top of the filter
	 */
	double getFilterTopPos() const;
	
	/**
	 * Advance the time of the simulation by dt
	 */
	void advanceTime(double dt);
	
	/**
	 * Constructor for simulation
	 */
	BioChemWeighted(double dt, double dx, double W, double D, 
					  double H, double DC, double flux0, double fluxL, 
					  double k, double theta);
	
	/**
	 * Destructor for simulation
	 */
	~BioChemWeighted();
	
private:
	
	/**
	 * No copy constructor
	 */
	BioChemWeighted(const BioChemWeighted&){}

	/**
	 * No assignment operator
	 */
	void operator=(const BioChemWeighted&){}
	
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
	
	/**
	 * Boyden Chamber (BC) specs for a single well
	 */
	double BC_hf;    // height (thickness) of filter
	double BC_hl;    // height of lower chamber
	double BC_hu;    // height of upper chamber
	double BC_d;     // diameter of well
	double BC_r;     // radius of well
	double BC_a;     // area of a single well
	/**
	 * Experimental (EX) specs
	 */
	int    EX_nhpf;  // number of high power fields (HPF) counted
	double EX_ahpf;  // area of a single HPF
	double EX_at;    // total area of all HPFs
	double EX_cc;    // concentration of cells
	double EX_chpf;  // number of cells per HPF
	double EX_c0;    // initial concentration of chemoattractant
};

#endif
