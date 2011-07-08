/*
 *  QModelQueues.h
 *  CellMigration
 */
#ifndef _QModelQueues_h_
#define _QModelQueues_h_
#include "QCell.h"
#include <list>

using namespace std;

class QModelQueues{

public:
	
	/**
	 * Return the list in the queue
	 */
	virtual int getListSize() const = 0;

	/**
	 * Virtual destructor
	 */
	virtual ~QModelQueues(){}
	
	/**
	 * Default constructor
	 */
	QModelQueues(){}
	
protected:
	
//	/**
//	 * No public assignment operator
//	 */
//	void operator=(const QModelQueues&){}
	
//	/**
//	 * No public copy constructor
//	 */
//	QModelQueues(const QModelQueues&){}

};

#endif