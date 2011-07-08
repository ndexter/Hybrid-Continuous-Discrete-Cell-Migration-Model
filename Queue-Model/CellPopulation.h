#ifndef _CellPopulation_h_
#define _CellPopulation_h_
#include "adevs.h"
#include "Cell.h"
#include <vector>

/**
 * This class is a container for the entire cell population.
 * At present, it only holds the cells and acts as a top level
 * model for the simulation engine to act on. Extended versions
 * will manage cell-cell interactions.
 */
class CellPopulation:
public adevs::Network<vec>
{
public:
	/**
	 * Population iterator.
	 */
	typedef std::vector<CellInterface*>::iterator iterator;
	
	/**
	 * Get the singleton CellPopulation object.
	 */
	static CellPopulation* getInstance();
	
	/**
	 * Get an iterator over the population.
	 */
	iterator begin() { return pop.begin(); }
	
	/**
	 * Get an iterator over the population.
	 */
	iterator end() { return pop.end(); }
	
	/**
	 * Get the component models for this cell space. The
	 * simulator needs this method; it is inherited from the
	 * Network base class.
	 */
	void getComponents(adevs::Set<adevs::Devs<vec>*>& c);
	
	/**
	 * Event routing is not needed for this model because
	 * the cells do not interact.
	 */
	void route(const vec&, adevs::Devs<vec>*, adevs::Bag<adevs::Event<vec> >&){}

private:
	
	/**
	 * Private constructor.
	 */
	CellPopulation();
	
	/**
	 * No copy constructor
	 */
	CellPopulation(const CellPopulation&){}
	
	/**
	 * No assignment operator
	 */
	void operator=(const CellPopulation&){}
	
	/**
	 * Private destructor
	 */
	~CellPopulation();
	
	/**
	 * The population itself!
	 */
	std::vector<CellInterface*> pop;
};

#endif
