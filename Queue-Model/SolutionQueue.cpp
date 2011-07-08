/*
 *  SolutionQueue.cpp
 *  CellMigration
 */

#include "SolutionQueue.h"

SolutionQueue::SolutionQueue(int initPop, double time){
	initialPopulation = initPop;
	t = time;
	QCell sample;
	for(int i = 0; i < initialPopulation; i++){
		sample = QCell(t);
		solutionlist.push_back(sample);
	}
}

double SolutionQueue::nextTime(int count){
	// Inverse of curve fit sigmoidal function
	double retval;
	if(count == 0){ 
		retval = DBL_MAX; 
	} else if(count == initialPopulation || initialPopulation - count <= 771){
		retval = 1/DBL_MAX; 
	} else { 
		retval = -900*log((double)initialPopulation/((double)initialPopulation - (double)count)-1)+5200;
	}
	

	return retval;
}

int SolutionQueue::getListSize() const{
	return solutionlist.size();
}

void SolutionQueue::delta_int(){
	t += ta();
	solutionlist.pop_front();
}

void SolutionQueue::output_func(adevs::Bag<QCell>& yb){
	yb.insert(solutionlist.front());
}

double SolutionQueue::ta(){
	double nt = nextTime(solutionlist.size()) - t;
	return nt;
}

SolutionQueue::~SolutionQueue(){}