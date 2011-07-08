#include "CellPopulation.h"
using namespace std;
using namespace adevs;

/* Get the singleton instance */
CellPopulation* CellPopulation::getInstance()
{
	/* Create the singleton instance */
	static CellPopulation inst;
	return &inst;
}

CellPopulation::CellPopulation():
Network<vec>()
{
	// Start near the far edge
	vec dims = BioChem::getInstance()->getDomain();
	double startPos = BioChem::getInstance()->getFilterTopPos();
	// This loop creates the cells and adds them to
	// the population set
	vec pstart;
	double tstart;
	int istart = 1;
	int inum = 6;
	int jnum = 6;
	for (int j = 0; j < jnum; j++){
		for (int i = 0; i < inum; i++){
			pstart[0] = dims[0]*(i+1)/(inum+1);
			pstart[1] = dims[1]*(j+1)/(jnum+1);
			pstart[2] = startPos;
			if (istart <= 1 ) tstart = (double)0.5;
			else if(istart <= 8) tstart = (double)1.0;
			else if(istart <= 20) tstart = (double)1.5;
			else if(istart <= 32) tstart = (double)2.0;
			else if(istart <= 34) tstart = (double)3.0;
			else tstart = (double)4.0;
			istart = istart + 1;
//			static rv r(543534543);
//			pstart[0] = r.uniform(startPos,dims[0]);
//			pstart[1] = r.uniform(startPos,dims[1]);
//			pstart[2] = r.uniform(startPos,dims[2]);
			pop.push_back(new Cell(pstart,tstart));
		}
	}
}

void CellPopulation::getComponents(adevs::Set<adevs::Devs<vec>*>& c)
{
	for (unsigned i = 0; i < pop.size(); i++)
	{
		c.insert(dynamic_cast<Devs<vec>*>(pop[i]));
	}
}

CellPopulation::~CellPopulation()
{
	for (unsigned i = 0; i < pop.size(); i++)
	{
		delete pop[i];
	}
}
