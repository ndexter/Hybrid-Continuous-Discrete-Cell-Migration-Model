#ifndef _CellInterface_h_
#define _CellInterface_h_
#include "vec.h"

/**
 * This interface must be implemented by every cell.
 * It is used to access visible cell features that
 * can impact the biochemical environment or other
 * cells.
 * */
class CellInterface
{
    public:
	/**
	 * Get the time at which the cell state last changed
	 * */
	virtual double getTime() const = 0;
    /**
     * Get the postion of the cell at getTime()
     * */
	virtual vec getPos() const = 0;
	/**
	 * Get the velocity of the cell at getTime()
	 * */
	virtual vec getVel() const = 0;
	/**
	 * Get the radius of the cell
	 * */
	virtual double getRadius() const = 0;
	/**
	  * Get the MMP production rate given the specified concentrations
	  * of chemo- and hapto- attractants at the cells location.
	  * The return value should be > 0 for production of MMPs and < 0
	  * for the destruction of MMPs
	  * */
	virtual double getMMP_Prod(double hapc, double chemc) const = 0;
	/**
	  * Virtual destructor
	  * */
	virtual ~CellInterface(){}
	/**
	  * Default constructor
	  * */
	CellInterface(){}
    protected:
	/**
	  * No public assignment operator
	  * */
	void operator=(const CellInterface&){}
	/**
	  * No public copy constructor
	  * */
	CellInterface(const CellInterface&){}
};

#endif
