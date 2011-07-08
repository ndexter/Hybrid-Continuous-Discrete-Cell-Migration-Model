/*
 *  QModelListener.h
 *  CellMigration
 */
#ifndef _QModelListener_h_
#define _QModelListener_h_
#include "adevs.h"
#include "QCell.h"
#include <fstream>

class QModelListener:
public adevs::EventListener<QCell>{
	
public:
	
	QModelListener(bool final_only = true);
	
	void outputEvent(adevs::Event<QCell> x, double t);
	
	~QModelListener();
	
private:
	
	std::ofstream fout;
	
	char scratch[1000];
	
	bool final_only;
	
};

#endif