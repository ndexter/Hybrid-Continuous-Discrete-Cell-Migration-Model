/*
 *  SolutionQueue.h - a class that represents a solution full of cells
 *  CellMigration
 */

#ifndef _SolutionQueue_h_
#define _SolutionQueue_h_
#include "QModelQueues.h"
#include "adevs.h"
#include "QCell.h"
#include <list>
#include <cmath>

using namespace std;
using namespace adevs;

class SolutionQueue: 
public adevs::Atomic<QCell>,
public QModelQueues{
	
public:
	
	SolutionQueue(int initPop, double time);
	
	/**
	 * An internal event, the cell dumping cells into the 
	 * waiting queue.
	 */
	void delta_int();
	
	/**
	 * Unimplemented because the Solution has no external 
	 * influence potential.
	 */
	void delta_ext(double, const adevs::Bag<QCell>&){}
	
	/**
	 * Unimplemented because if the Solution has no external 
	 * influence potential, confluent events don't happen.
	 */
	void delta_conf(const adevs::Bag<QCell>&){}
	
	/**
	 * Output the cells moving to the WaiterQueue
	 */
	void output_func(adevs::Bag<QCell>& yb);
	
	double nextTime(int count);
	
	int getListSize() const;
	
	/**
	 * time advance function
	 */
	double ta();
	
	void gc_output(adevs::Bag<QCell>&){}
	
	~SolutionQueue();
	
private:
	
	int initialPopulation;
	
	list<QCell> solutionlist;
	
	QCell holder;
	
	double t;
	
};

#endif