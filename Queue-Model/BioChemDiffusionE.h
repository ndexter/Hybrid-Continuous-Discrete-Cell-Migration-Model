#ifndef BIOCHEMDIFFUSIONE_H_
#define BIOCHEMDIFFUSIONE_H_
#include "BioChemInstance.h"

/**
 * This class implements a diffusing gradient 
 */
class BioChemDiffusionE:
	public BioChemInstance
{
	public:
	
	/**
	 * Print out the header for the .vtk VisIt files
	 */
	std::string VTKHeader() const;
	
	/**
	 * Print the concentration at each point on the grid
	 */
	std::string printEnvironment() const;
	
	/**
	 * Print the file that models the filter
	 */
	std::string printFilter() const;
	
	/**
	  * Get the haptoattractant gradient at position p.
	  * */
	vec getHaptacGrad(vec p) const { return vec(0.0); }
	/**
	  * Get the chemotattractant gradient at position p
	  * */
	vec getChemtacGrad(vec p) const;
	/**
	  * Get the haptoattractant concentration at position p
	  * */
	double getHaptac(vec p) const { return 0.0; }
	/**
	  * Get the chemoattractant concentration at position p
	  * */
	double getChemtac(vec p) const;
	/**
	  * Advance the explicit diffusion solution
	  * */
	void advanceDiffusion(double timestep);
	/**
	  * Get the time at which the hapto- and chemo- attractant
	  * quantities are valid for (i.e., the last update time)
	  * */
	double getTime() const { return t; }
	/**
	  * Get the extent of the solution domain in each dimension
	  * */
	vec getDomain() const { return L; }
	/**
	 * Get the 'natural' grid spacing for producing output.
	 */
	double getGridSpacing() const { return dx; }
	/**
	  * Advance the solution by dt units of time.
	  * */
	void advanceTime(double dt);
	/**
	  * Virtual destructor.
	  * */
	~BioChemDiffusionE(){}
	/*
	 * Get the Diffusion Coefficient
	 */
	double getDiffCoeff() const { return DifCoef; }
	/*
	 * get the initial concentration of substance
	 */
	double getInitialConc() const { return a; }
	/**
	 * Constructor. The parameter a is the initial concentration at the far end of the chamber.
	 * The parameter dx is the grid resolution. The parameter l is the length (in meters) of the
	 * boyden chamber. The DifCoef is the diffusion coefficient.
	 */
	BioChemDiffusionE(double a, double dx, vec L, double DifCoef);

	private:

	double* c[2];
	int ArraySize;
	int current;
	int next;
	double a, dx, t, DifCoef, h; 
	vec L;
};

#endif /*BIOCHEMDIFFUSIONE_H_*/
