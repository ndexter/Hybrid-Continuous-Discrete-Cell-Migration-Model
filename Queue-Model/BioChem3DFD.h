#ifndef _BioChem3DFD_h_
#define _BioChem3DFD_h_
#include "BioChemInstance.h"
#include "adevs.h"
#include <cmath>
#include <sys/time.h>
extern "C"{
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_permutation.h>
}

class BioChem3DFD:
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
		double getTime() const { return t; }
		
		/**
		 * Get the grid spacing
		 */
		double getGridSpacing() const { return dx; }
		
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
		void buildLHS(gsl_matrix *LHS);
		
		/**
		 * Generate the matrix equations (RHS) necessary to solve the explicit solution
		 */
		void buildRHS(gsl_matrix *RHS);
		
		/**
		 * Get the Diffusion Coefficient
		 */
		double getDiffCoeff() const { return DC; }
		
		/**
		 * get the initial concentration of substance
		 */
		double getInitialConc() const { return EX_c0; }
		
		/**
		 * Get the extent of the solution domain in each dimension
		 */
		vec getDomain() const;
		
		/**
		 * Get the position of the bottom of the filter
		 */
		double getFilterBottomPos() const {	return fbottom; }
		
		/**
		 * Get the position of the top of the filter
		 */
		double getFilterTopPos() const { return ftop; }	
		
		/**
		 * Print out the header for the .vtk VisIt files
		 */
		std::string VTKHeader() const;
		
		/**
		 * Prints out the concentration at all points on the grid
		 */
		std::string printEnvironment() const;
		
		/**
		 * Print out the file that models the filter
		 */
		std::string printFilter() const;
		
		/**
		 * Advance the time of the simulation by dt
		 */
		void advanceTime(double dt);
		
		/**
		 * Constructor for simulation
		 */
		BioChem3DFD();
		
		/**
		 * Destructor for simulation
		 */
		~BioChem3DFD();
		
	private:
		
		/**
		 * No copy constructor
		 */
		BioChem3DFD(const BioChem3DFD&){}
		
		/**
		 * No assignment operator
		 */
		void operator=(const BioChem3DFD&){}
		
		// 3D coordinate 1D vectors
		gsl_vector *x;		// u at n+1
		gsl_vector *b;		// C*(u at n)
		gsl_vector *un;		// u at n
		
		// Central Difference matrices 
		gsl_matrix *LHS;	// Left Hand Side, Implicit solution
		gsl_matrix *RHS;	// Right Hand Side, Explicit solution
		
		// Permutation for calculating LU Decomposition
		gsl_permutation *p;
		
		/**
		 * Bottom and top of filter
		 */
		double ftop, fbottom;
		
		// Length of the domain
		double Width,Depth,Height;
		
		// Diffusion coefficient, Diffusion coefficient (filter)
		double DC, DCf;
		
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
		double alpha, alphaprev;
		
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
