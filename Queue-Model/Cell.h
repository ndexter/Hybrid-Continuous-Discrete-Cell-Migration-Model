#ifndef _Cell_h_
#define _Cell_h_
#include "CellInterface.h"
#include "BioChem.h"
#include "adevs.h"
#include "vec.h"
#include <math.h>


/**
 * This class describes the dynamics of a single cell moving
 * along a chemotactic and haptotactic gradient.
 */
class Cell:
public adevs::Atomic<vec>,
public CellInterface
{
public:
	
	/**
	 * Constructor sets the cell parameters. The
	 * initial position must be supplied.
	 */
	Cell(vec p0);
	
	Cell(vec p0, double t0);  // supply time attached
	
	/**
	 * Internal transition funtion.
	 */
	void delta_int();
	
	/**
	 * Cells are autonomous in this model. They
	 * perceive the environment only periodically.
	 */
	void delta_ext(double,const adevs::Bag<vec>&){}
	
	/**
	 * Cells are autonomous, so no confluent implementation.
	 */
	void delta_conf(const adevs::Bag<vec>&){}
	
	/**
	 * Cell output is its current position. This is used
	 * to plot the cell motion.
	 */
	void output_func(adevs::Bag<vec>& yb);
	
	double getDistance(vec p1, vec p2);
	
	/**
	 * Time advance function
	 */
	double ta();
	
	/**
	 * No garbage collection needed for output.
	 */
	void gc_output(adevs::Bag<vec>&){}
	
	/**
	 * Destructor.
	 */
	~Cell();
	
	/*
	 * Cell interface method.
	 */
	double getTime() const { return getLastEventTime(); }
	
	double getTimeAttached() const { return tstart; }
	
	vec getPos() const { return p; }
	
	vec getVel() const { return v; }
	
	double getRadius() const { return 0.0; }
	
	double getMMP_Prod(double hpac, double chemc) const { return 0.0; }
	
	double getSpd(void);
	
	double getF(void); 
	
private:
	
	/**
	 * No default constructor.
	 */
	Cell(){}
	
	/**
	 * No copy constructor
	 */
	Cell(const Cell&){}
	
	/**
	 * No assignment operator.
	 */
	void operator=(const Cell&){}
	
	/**
	 * Cell position in the space
	 */
	vec p;
	
	/**
	 * Cell sampling freq., speed, and random motility coefficient
	 */
	double f, spd, gain;
	
	/**
	 * Cell direction of travel
	 */
	vec v;
	
	/**
	 * Cell attachment time
	 */
	double tstart;
};

#endif

