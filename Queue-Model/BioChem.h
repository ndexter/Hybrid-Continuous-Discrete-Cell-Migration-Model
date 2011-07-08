#ifndef _BioChem_h_
#define _BioChem_h_
#include "BioChemInstance.h"
#include <iostream>

/**
This class represents the evolving biochemical environment
in which the cell lives. It is a singleton class because 
any simulation run can have only one biochemical environment.
This model includes methods for accessing chemotactic and
haptotactic chemical concentrations. It uses the CellInterface
to incorporate the effect of cells on the biochemical environment
(i.e., the product and degradation of MMP2 by the cell population).
*/
class BioChem
{
    public:
	
	/**
	  * Get the unique instance of this object
	  */
	static BioChem* getInstance();
	
	/**
	  * Get the haptoattractant gradient at position p
	  */
	vec getHaptacGrad(vec p) const;
	
	/**
	  * Get the chemotattractant gradient at position p
	  */
	vec getChemtacGrad(vec p) const;
	
	/**
	  * Get the haptoattractant concentration at position p
	  */
	double getHaptac(vec p) const;
	
	/**
	  * Get the chemoattractant concentration at position p
	  */
	double getChemtac(vec p) const;
	
	/**
	  * Get the time at which the hapto- and chemo- attractant
	  * quantities are valid for (i.e., the last update time)
	  */
	double getTime() const;
	
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
	 * Get the 'natural' grid spacing for producing output.
	 */
	double getGridSpacing() const;
	
	/**
	  * Advance the solution by dt units of time.
	  */
	void advanceTime(double dt);
	
	/**
	 *  Get the Diffusion coefficient for the chemoattractant
	 */ 
	double getDiffCoeff() const;
	
	/**
	 * get the Initial Concentration in the lower? well
	 */
	double getInitialConc() const;
	
	/**
	 * Print out the header for VisIt .vtk files
	 */
	std::string VTKHeader() const;

	/**
	 * Print out the concentration at all points on the grid
	 */
	std::string printEnvironment() const;

	/**
	 * Print out the file that models the filter
	 */
	std::string printFilter() const;
	
    private:
	
	/**
	  * Object can only create instances of itself
	  */
	BioChem();
	
	/**
	  * The object can only destroy itself
	  */
	~BioChem();
	
	/**
	  * No copy constructor
	  */
	BioChem(const BioChem&){}
	
	/**
	  * No assignment operator
	  */
	void operator=(const BioChem&){}
	
	/**
	 * Instance variable
	 */
	BioChemInstance* impl;
};

#endif
