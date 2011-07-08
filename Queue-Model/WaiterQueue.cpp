/*
 *  WaiterQueue.cpp
 *  CellMigration
 */

#include "WaiterQueue.h"
#include <cmath>

void WaiterQueue::delta_int(){
	// decrement #pores
//	poreCount--; 
	// pop cell off of x
	waiterlist.pop_front();
}

void WaiterQueue::delta_ext(double e, const Bag<QCell>& xb){
	
	// sample of cell coming in
	QCell sample;
	
	if(!xb.empty()){
		// add everyting in xb to x
		Bag<QCell>::iterator iter = xb.begin();
		while(iter != xb.end()){
			sample = *iter;
			if(!sample.isCellAttached()){
				sample.Attach();
				waiterlist.push_back(sample);
				if(poreCount > 0){
					eject = true;
					poreCount--;
				}
			} else {
				poreCount++;
				if(waiterlist.size() > 0){
					eject = true;
					poreCount--;
				}
			}
			iter++;
		}
	}
}

int WaiterQueue::getListSize() const{
	return waiterlist.size();
}

void WaiterQueue::delta_conf(const Bag<QCell>& xb){
	delta_int();
	delta_ext(0.0, xb);
}

void WaiterQueue::output_func(Bag<QCell>& yb){
	yb.insert(waiterlist.front());
	eject = false;
}

double WaiterQueue::ta() { 
	if (eject){ 
		return 0; 
	} else { 
		return DBL_MAX; 
	} 
}