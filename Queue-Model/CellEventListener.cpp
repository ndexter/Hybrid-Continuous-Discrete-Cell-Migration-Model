#include "CellEventListener.h"
#include "CellPopulation.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
using namespace std;
using namespace adevs;

static const char* cellDir = "cells";

CellEventListener::CellEventListener(bool final_only):
	EventListener<vec>(),
	final_only(final_only)
{
	errno = 0;
	mkdir(cellDir,S_IRWXU);
	if (errno != 0)
		throw errno;
}

void CellEventListener::cleanup(double t)
{
	 /* Create a gnuplot file for the cell trajectories 
	 * This is done at the end to make sure any cells created
	 * during the simulation are captured in the plot list 
	 */
	sprintf(scratch,"%s/plot.gnu",cellDir);
	fout.open(scratch);
	fout << "set ylabel 'x'" << endl;
	fout << "set xlabel 't'" << endl;
	fout << "unset key" << endl;
	CellPopulation* pop = CellPopulation::getInstance();
	CellPopulation::iterator iter = pop->begin();
	sprintf(scratch,"'%lx'",(unsigned long)(*iter));
	fout << "plot " << scratch << " with lines";
	iter++;
	for (; iter != pop->end(); iter++)
	{
		sprintf(scratch,"'%lx'",(unsigned long)(*iter));
		fout << ", " << scratch << " with lines";
	}
	fout << endl;
	fout.close();
	/*
	 * Dump the final cell positions
	 */
	sprintf(scratch,"%s/final",cellDir);
	ofstream final(scratch);
	for (iter = pop->begin(); iter != pop->end(); iter++)
	{
		vec position = (*iter)->getPos()+
			(*iter)->getVel()*(t-(*iter)->getTime());
		Event<vec> x(dynamic_cast<Devs<vec>*>(*iter),position);
		outputEvent(x,t);
                // KLK - write y & z position also
		final << position[0] << " " << position[1] << " " << position[2] << endl;
	}
	final.close();
}

CellEventListener::~CellEventListener()
{
}

void CellEventListener::outputEvent(Event<vec> x, double t)
{
	if (final_only) return;
	CellInterface* iface = dynamic_cast<CellInterface*>(x.model);
	if (iface == NULL) return;
	sprintf(scratch,"%s/%lx",cellDir,(unsigned long)(iface));
	fout.open(scratch,ios_base::app);
        // KLK - write y & z position also
	fout << t << " " << x.value[0] << " " << x.value[1] << " " << x.value[2] << endl;
	fout.close();
}
