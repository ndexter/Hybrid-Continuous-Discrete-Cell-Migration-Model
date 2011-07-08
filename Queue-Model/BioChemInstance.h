#ifndef _BioChemInstance_h_
#define _BioChemInstance_h_
#include <sstream>
#include <string>
#include "vec.h"
#include "CellInterface.h"

/**
 * This is an interface class that must be implemented
 * by an actual biochemical process model. One particular
 * biochemical process model is used by the singleton
 * BioChem class to simulate the dynamics of a particular
 * experiment. Consequently, only one implementation
 * of the BioChemInterface will be used in any simulation run.
 * Which one is used can be changed by modifying the 
 * constructor of the BioChem class.
*/
class BioChemInstance
{
	public:
	
	/**
	  * Get the haptoattractant gradient at position p.
	  */
	virtual vec getHaptacGrad(vec p) const = 0;
	
	/**
	  * Get the chemotattractant gradient at position p
	  */
	virtual vec getChemtacGrad(vec p) const = 0;
	
	/**
	  * Get the haptoattractant concentration at position p
	  */
	virtual double getHaptac(vec p) const = 0;
	
	/**
	  * Get the chemoattractant concentration at position p
	  */
	virtual double getChemtac(vec p) const = 0;
	
	/**
	  * Get the time at which the hapto- and chemo- attractant
	  * quantities are valid for (i.e., the last update time)
	  */
	virtual double getTime() const = 0;
	
	/**
	  * Get the extent of the solution domain in each dimension
	  */
	virtual vec getDomain() const = 0;
	
	/**
	  * Get the position of the bottom of the filter
	  */
	virtual double getFilterBottomPos() const = 0;
	
	/**
	  * Get the position of the top of the filter
	  */
	virtual double getFilterTopPos() const = 0;
	
	/**
	 * Get the 'natural' grid spacing for producing output.
	 */
	virtual double getGridSpacing() const = 0;
	
	/*
	 * get the diffusion coefficient for the chemoattractant
	 */
	virtual double getDiffCoeff() const = 0;
	
	/**
	 * get the initial concentration of chemoattractant in the lower well
	 */
	virtual double getInitialConc() const = 0;
	
	/**
	  * Advance the solution by dt units of time.
	  */
	virtual void advanceTime(double dt) = 0;
	
	/**
	 * Print out the header for VisIt .vtk files
	 */ 
	virtual std::string VTKHeader() const = 0;
	
	/**
	 * Print out the concentration at each point on the grid
	 */ 
	virtual std::string printEnvironment() const = 0;
	
	/**
	 * Print out the file that models the filter
	 */
	virtual std::string printFilter() const = 0;
	
	/**
	  * Virtual destructor.
	  */
	virtual ~BioChemInstance(){}
};

#endif
