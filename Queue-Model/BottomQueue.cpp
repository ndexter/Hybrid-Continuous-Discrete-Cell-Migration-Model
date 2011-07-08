/*
 *  BottomQueue.cpp
 *  CellMigration
 */

#include "BottomQueue.h"

void BottomQueue::delta_ext(double e, const adevs::Bag<QCell>& xb){
	QCell holder;
	
	if(!xb.empty()){
		adevs::Bag<QCell>::iterator iter = xb.begin();
		while(iter != xb.end()){
			holder = *iter;
			bottomlist.push_back(holder);
			iter++;
		}
	}
}

int BottomQueue::getListSize() const{
	return bottomlist.size();
}