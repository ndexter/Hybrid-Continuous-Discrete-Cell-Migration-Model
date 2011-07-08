#ifndef _BioChem1DFD_h_
#define _BioChem1DFD_h_
#include "BioChemInstance.h"
#include "adevs.h"
#include <cmath>
#include <sys/time.h>
extern "C"{
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_linalg.h>
}

class BioChem1DFD:
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
		double getGridSpacing() const {	return dx; }
		
		/**
		 * Generate the matrix equations (LHS) necessary to solve the implicit solution
		 */
		void buildLHS(gsl_vector *lhsdiag, gsl_vector *lhsabovediag, gsl_vector *lhsbelowdiag);
		
		/**
		 * Generate the matrix equations (RHS) necessary to solve the explicit solution
		 */
		void buildRHS(gsl_vector *rhsdiag, gsl_vector *rhsabovediag, gsl_vector *rhsbelowdiag);
		
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
		BioChem1DFD();
		
		/**
		 * Destructor for simulation
		 */
		~BioChem1DFD();
		
	private:
		
		/**
		 * No copy constructor
		 */
		BioChem1DFD(const BioChem1DFD&){}
		
		/**
		 * No assignment operator
		 */
		void operator=(const BioChem1DFD&){}
		
		// 3D coordinate 1D vectors
		gsl_vector *x;		// u at n+1
		gsl_vector *b;		// C*(u at n)
		gsl_vector *un;		// u at n
		
		// Central Difference diagonals
		gsl_vector *lhsdiag, *lhsabovediag, *lhsbelowdiag;
		gsl_vector *rhsdiag, *rhsabovediag, *rhsbelowdiag;
		
		/**
		 * Bottom and top of filter
		 */
		double ftop, fbottom;
		
		// Length of the domain
		double Width, Depth, Height;
		
		// Diffusion coefficient, Diffusion coefficient (filter)
		double DC, DCf;
		
		// Current time, time step
		double t, dt;
		
		// Grid step size
		double dx;
		
		// Initial Concentration
		double c0;
		
		// Combined method switch
		double theta;
		
		// Flux constants at boundary
		double flux0, fluxL;
		
		// Flux scalar
		double k;
		
		// Diffusion term
		double alpha, alphaprev;
		
		double EX_c0;    // initial concentration of chemoattractant
	};

#endif