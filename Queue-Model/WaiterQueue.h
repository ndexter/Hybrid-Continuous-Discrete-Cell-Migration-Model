/*
 *  WaiterQueue.h
 *  CellMigration
 */

#ifndef _WaiterQueue_h_
#define _WaiterQueue_h_
#include "QModelQueues.h"
#include "adevs.h"
#include "QCell.h"
#include <cmath>
#include <list>

using namespace std;
using namespace adevs;

class WaiterQueue:
public adevs::Atomic<QCell>,
public QModelQueues{
	
public:
	
	WaiterQueue(double pores){
		poreCount = pores;
	}
	
	/**
	 * Internal state transition
	 */
	void delta_int();
	
	/**
	 * External state transition
	 */
	void delta_ext(double, const adevs::Bag<QCell>&);

	/**
	 * Confluent state transition
	 */
	void delta_conf(const adevs::Bag<QCell>&);
	
	/**
	 * Output the cells moving to the WigglerQueue
	 */
	void output_func(adevs::Bag<QCell>& yb);
	
	/**
	 * time advance function
	 */ 
	double ta();
	
	int getListSize() const;
	
	void gc_output(adevs::Bag<QCell>&){}

private:
	
	// internal list for the QCells
	list<QCell> waiterlist;
	
	double poreCount;
	
	// current time
	double t;
	
	// Holder QCell
	QCell holder;
	
	// eject from waiter queue switch
	bool eject;

};

#endif