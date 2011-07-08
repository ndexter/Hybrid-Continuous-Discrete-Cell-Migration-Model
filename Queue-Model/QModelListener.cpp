/*
 *  QModelListener.cpp
 *  CellMigration
 *
 *  Created by Nicholas Dexter on 3/30/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "QModelListener.h"
#include "QModelQueues.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>

using namespace std;
using namespace adevs;

static const char* cellDir = "qcells";

QModelListener::QModelListener(bool final_only):
EventListener<QCell>(),
final_only(final_only){
	errno = 0;
	mkdir(cellDir,S_IRWXU);
	if(errno != 0)
		throw errno;
}

//void QModelListener::cleanup(double t){
////	sprintf(scratch, "%s/migration_statistics",cellDir);
////	fout.open(scratch);
////	fout << "This is where migration statistics go";
////	fout.close();
//}

void QModelListener::outputEvent(Event<QCell> x, double t){
//	if(final_only) return;
//	QModelQueues* qiface = dynamic_cast<QModelQueues*> (x.model);
//	if(qiface == NULL) return;
//	sprintf(scratch, "%s/%d.vtk", cellDir, t);
//	fout.open(scratch,ios_base::app);
//	fout << t << " " << x.value.isCellAttached() << " " << x.value.isCellMigrating() << " " 
//		 << x.value.getspawntime() << " " << x.value.getattachtime() << " " 
//		 << x.value.getstarttime() << " " << x.value.getstoptime() << "\n";
//	fout.close();
}

QModelListener::~QModelListener(){
	fout.close();
	delete scratch;
}