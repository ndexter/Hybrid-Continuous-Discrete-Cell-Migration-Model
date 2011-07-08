/*
 *  QCell.h - A class that represents a cell in a queue
 *  CellMigration
 *
 *  Created by Nicholas Dexter on 3/25/09.
 *  Copyright 2009. All rights reserved.
 *
 */
#ifndef _QCell_h_
#define _QCell_h_
#include <cmath>

class QCell
{
public:
	
	QCell(){
		spawntime = 0.0;
	}
	
	// cell created with an event time
	QCell(double t){
		spawntime = t;
		isAttached = false;
		isMigrating = false;
	}
	
	// cell attaches to a surface
	void Attach(){
		isAttached = true;
	}
	
	// cell begins to migrate
	void Migrate(double t){
		starttime = t;
		isMigrating = true;
	}
	
	void stopMigratingTime(double t){
		stoptime = t;
	}
	
	// cell stops migrating
	void stopMigrating(){
		isMigrating = false;
	}
	
	QCell operator=(const QCell& src){
		isAttached = src.isAttached;
		isMigrating = src.isMigrating;
		spawntime = src.spawntime;
		starttime = src.starttime;
		stoptime = src.stoptime;
		return *this;
	}
	
	bool operator==(const QCell& src) const{
		return ((isAttached && src.isAttached) &&
				(isMigrating && src.isMigrating) &&
				(spawntime == src.spawntime) &&
				(starttime == src.starttime) &&
				(stoptime == src.stoptime));
	}
	
	bool operator!=(const QCell& src) const{
		return !(*this == src);
	}
	
	bool operator<(const QCell& src) const{
		return stoptime > src.stoptime;
	}
	
	bool isCellAttached() const{ return isAttached; }

	bool isCellMigrating() const{ return isMigrating; }
	
	double getspawntime() const{ return spawntime; }
	
	double getstarttime() const{ return starttime; }
	
	double getstoptime() const{ return stoptime; }
	
	// cell dies
	~QCell(){}
	
private:
	
	// Last time the cell was part of an event, either the start 
	// of the simulation or any change in state 
	double spawntime, starttime, stoptime;
	
	// Is the cell attached yet, is it migrating?
	bool isAttached, isMigrating;
};

#endif