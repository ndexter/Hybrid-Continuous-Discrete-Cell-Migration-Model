#include "BioChemOutput.h"
#include "BioChem.h"
#include "CellEventListener.h"
#include "Cell.h"
#include "CellInterface.h"
#include "CellPopulation.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
#include <iostream>

using namespace std;
using namespace adevs;

static const char* chemDir = "biochem";
static const char* cellDir = "cell";

ofstream fout;

BioChemOutput::BioChemOutput(double t_rec):
Atomic<vec>()
{
	/* Create the output directory */
	errno = 0;
	mkdir(chemDir,S_IRWXU);
	mkdir(cellDir,S_IRWXU);
	if (errno != 0)
		throw errno;
	/* Remember the output recording interval */
	dt = t_rec;
	fnum = 0;
}

double BioChemOutput::ta()
{
	return dt;
}

void BioChemOutput::output_initial(){
//	
//	// Print out Biochemical Data
//	char tmp[100];
//	sprintf(tmp, "%s/%d.vtk",chemDir, fnum);
//	fout.open(tmp);
	BioChem* biochem = BioChem::getInstance();
//	fout << biochem->VTKHeader();
//	fout << biochem->printEnvironment();
//	fout.close();
//	
//	// Print out Cell Data
//	char tmp2[100];
//	sprintf(tmp2, "%s/%d.vtk",cellDir, fnum);
//	fout.open(tmp2);	
//	CellPopulation* pop = CellPopulation::getInstance();
//	CellPopulation::iterator iter;
//	fout<< "# vtk DataFile Version 2.0\n"
//		<< "3D Cells\n"
//		<< "ASCII\n"
//		<< "DATASET UNSTRUCTURED_GRID \n"
//		<< "POINTS 36 float\n";
//	for(iter = pop->begin(); iter != pop->end(); iter++){
//		vec position = (*iter)->getPos()+
//		(*iter)->getVel()*(biochem->getTime()-(*iter)->getTime());
//		Event<vec> x(dynamic_cast<Devs<vec>*>(*iter),position);
//		fout << x.value[0] << " " << x.value[1] << " " << x.value[2] << endl;
//	}
//	fout<< "CELLS 36 72\n"
//		<< "1 35\n1 34\n1 33\n1 32\n1 31\n1 30\n1 29\n1 28\n1 27\n1 26\n1 25\n"
//		<< "1 24\n1 23\n1 22\n1 21\n1 20\n1 19\n1 18\n1 17\n1 16\n1 15\n1 14\n"
//		<< "1 13\n1 12\n1 11\n1 10\n1 9\n1 8\n1 7\n1 6\n1 5\n1 4\n1 3\n1 2\n"
//		<< "1 1\n1 0\n"
//		<< "CELL_TYPES 36\n"
//		<< "1 1 1 1 1 1\n"
//		<< "1 1 1 1 1 1\n"
//		<< "1 1 1 1 1 1\n"
//		<< "1 1 1 1 1 1\n"
//		<< "1 1 1 1 1 1\n"
//		<< "1 1 1 1 1 1\n"
//		<< "POINT_DATA 36\n"
//		<< "SCALARS scalars float 1\n"
//		<< "LOOKUP_TABLE default\n"
//		<< "1 2 3 4 5 6\n"
//		<< "7 8 9 10 11 12\n"
//		<< "13 14 15 16 17 18\n"
//		<< "19 20 21 22 23 24\n"
//		<< "25 26 27 28 29 30\n"
//		<< "31 32 33 34 35 36\n";
//	fout.close();
//	fnum++;
//	
	char tmp3[100];
	sprintf(tmp3, "filter.vtk");
	fout.open(tmp3);
	fout << biochem->printFilter();
	fout.close();
}

void BioChemOutput::output_func(Bag<vec>&)
{
	char tmp[100];
	sprintf(tmp,"%s/%d.vtk",chemDir,fnum);
	fout.open(tmp);
	BioChem* biochem = BioChem::getInstance();
	fout << biochem->VTKHeader();
	fout << biochem->printEnvironment();
	fout.close();

	//3/24/2009 note remove all below this line
	
//	// Print out Cell Data
//	char tmp2[100];
//	sprintf(tmp2, "%s/%d.vtk",cellDir, fnum);
//	fout.open(tmp2);	
//	CellPopulation* pop = CellPopulation::getInstance();
//	CellPopulation::iterator iter;
//	fout<< "# vtk DataFile Version 2.0\n"
//	<< "3D Cells\n"
//	<< "ASCII\n"
//	<< "DATASET UNSTRUCTURED_GRID \n"
//	<< "POINTS 36 float\n";
//	
//	// IS THIS RIGHT???
//	for(iter = pop->begin(); iter != pop->end(); iter++){
//		vec position = (*iter)->getPos()+
//		(*iter)->getVel()*(biochem->getTime()-(*iter)->getTime());
//		Event<vec> x(dynamic_cast<Devs<vec>*>(*iter),position);
//		fout << x.value[0] << " " << x.value[1] << " " << x.value[2] << endl;
//	}
//	fout<< "CELLS 36 72\n"
//	<< "1 35\n1 34\n1 33\n1 32\n1 31\n1 30\n1 29\n1 28\n1 27\n1 26\n1 25\n"
//	<< "1 24\n1 23\n1 22\n1 21\n1 20\n1 19\n1 18\n1 17\n1 16\n1 15\n1 14\n"
//	<< "1 13\n1 12\n1 11\n1 10\n1 9\n1 8\n1 7\n1 6\n1 5\n1 4\n1 3\n1 2\n"
//	<< "1 1\n1 0\n"
//	<< "CELL_TYPES 36\n"
//	<< "1 1 1 1 1 1\n"
//	<< "1 1 1 1 1 1\n"
//	<< "1 1 1 1 1 1\n"
//	<< "1 1 1 1 1 1\n"
//	<< "1 1 1 1 1 1\n"
//	<< "1 1 1 1 1 1\n"
//	<< "POINT_DATA 36\n"
//	<< "SCALARS scalars float 1\n"
//	<< "LOOKUP_TABLE default\n"
//	<< "1 2 3 4 5 6\n"
//	<< "7 8 9 10 11 12\n"
//	<< "13 14 15 16 17 18\n"
//	<< "19 20 21 22 23 24\n"
//	<< "25 26 27 28 29 30\n"
//	<< "31 32 33 34 35 36\n";
//	fout.close();
//	fnum++;
}


BioChemOutput::~BioChemOutput(){
}
