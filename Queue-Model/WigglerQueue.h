/*
 *  WigglerQueue.h
 *  CellMigration
 */
#ifndef _WigglerQueue_h_
#define _WigglerQueue_h_
#include "QModelQueues.h"
#include "adevs.h"
#include "QCell.h"
#include "vec.h"
#include "BioChem.h"
#include <cmath>
#include <queue>
#include <cassert>

using namespace std;
using namespace adevs;

class WigglerQueue:
public adevs::Atomic<QCell>,
public QModelQueues{
	
public:
	
	WigglerQueue(){}
	
	void delta_int();
	
	void delta_ext(double, const adevs::Bag<QCell>&);
	
	void delta_conf(const adevs::Bag<QCell>&);
	
	void output_func(adevs::Bag<QCell>&);
	
	double ta();

	int getListSize() const;
	
	double serviceTime(double z_grad);
	
	void gc_output(adevs::Bag<QCell>&){}
	
	~WigglerQueue(){}
	
private:

	// internal list for the QCells
	priority_queue<QCell> wigglerlist;
	
	// time
	double t;
	
	// # of pores, # of cells
	int poreCount;
	
	// eject from waiter queue switch
	bool eject;
};

#endif