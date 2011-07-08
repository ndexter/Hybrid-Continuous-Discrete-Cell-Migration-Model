/*
 *  BottomQueue.h
 *  CellMigration
 */
#ifndef _BottomQueue_h_
#define _BottomQueue_h_
#include "QModelQueues.h"
#include "QCell.h"
#include "adevs.h"
#include <list>

using namespace std;

class BottomQueue:
public adevs::Atomic<QCell>,
public QModelQueues{
	
public:
	
	BottomQueue(){}
	
	void delta_int(){}
	
	void delta_ext(double e, const adevs::Bag<QCell>& xb);
	
	void delta_conf(const adevs::Bag<QCell>&){}
	
	void output_func(adevs::Bag<QCell>&){}
	
	void gc_output(adevs::Bag<QCell>&){}
	
	int getListSize() const;
	
	double ta(){ return DBL_MAX; }
	
private:
	
	list<QCell> bottomlist;
	
};

#endif