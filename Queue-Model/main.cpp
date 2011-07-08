//#include "CellPopulation.h"
//#include "CellInterface.h"
//#include "Cell.h"
//#include "CellEventListener.h"
#include "BioChem.h"
#include "BioChemOutput.h"
#include "adevs.h"
#include "QModelQueues.h"
#include "QModelListener.h"
#include "SolutionQueue.h"
#include "WaiterQueue.h"
#include "WigglerQueue.h"
#include "BottomQueue.h"
#include <iostream>
#include <fstream>
#include <cerrno>
using namespace std;
using namespace adevs;

/**
 * Very top level simulation model that contains
 * the BioChemOutput process and the
 * CellPopulation model.
 */
/*
class TopLevelModel:public Network<vec>{
public:
	
	TopLevelModel(double biochem_t_rec);
	
	void route(const vec&, Devs<vec>*, Bag<Event<vec> >&){}
	
	void getComponents(Set<Devs<vec>*>& c);
	
	BioChemOutput* getBioChemOutput(){
		return bio_chem_output;
	}
	
	~TopLevelModel();
	
private:
	
	BioChemOutput* bio_chem_output;
	
};

TopLevelModel::TopLevelModel(double biochem_t_rec):Network<vec>(){
	bio_chem_output = new BioChemOutput(biochem_t_rec);
}

void TopLevelModel::getComponents(Set<Devs<vec>*>& c){
	c.insert(bio_chem_output);
	c.insert(CellPopulation::getInstance());
}

TopLevelModel::~TopLevelModel(){
	delete bio_chem_output;
}
*/

class QueueModel:public SimpleDigraph<QCell>{
public:
	
	QueueModel(double biochem_t_rec);
	
	BioChemOutput* getBioChemOutput(){
		return bio_chem_output;
	}
	
	string printSizes();
	
	string getServiceTime(double gradient);
	
	~QueueModel();
	
private:
	
	BioChemOutput* bio_chem_output;
	
	SolutionQueue* solution;
	
	WaiterQueue* waiters;
	
	WigglerQueue* wigglers;
	
	BottomQueue* bottom;
};

QueueModel::QueueModel(double biochem_t_rec):SimpleDigraph<QCell>(){
	bio_chem_output = new BioChemOutput(biochem_t_rec);
	
	solution = new SolutionQueue(250000, 0.0);
	waiters = new WaiterQueue(30000);
	wigglers = new WigglerQueue();
	bottom = new BottomQueue();
	
	add(solution);
	add(waiters);
	add(wigglers);
	add(bottom);
	
	couple(solution, waiters);
	couple(waiters, wigglers);
	couple(wigglers, waiters);
	couple(wigglers, bottom);
}

string QueueModel::printSizes(){
	ostringstream outstr;
	outstr << solution->getListSize()  << " " 
			<< waiters->getListSize() << " " 
			<< wigglers->getListSize() << " " 
			<< bottom->getListSize();
	return outstr.str();
}

string QueueModel::getServiceTime(double gradient){
	ostringstream outstr;
	outstr << wigglers->serviceTime(gradient);
	return outstr.str();
}

QueueModel::~QueueModel(){}

/**
 * Main simulation loop
 */
int main(){
//	CellEventListener* listener = NULL;
	QModelListener* listener = NULL;
//	Simulator<vec>* sim = NULL;
	Simulator<QCell>* sim = NULL;
//	TopLevelModel* model = NULL;
	QueueModel* qmodel = NULL;
	try
	{
		// Create the simulation objects
		
		// HVSMC total time in seconds
		double tend = (double) 14400;
		
		// 200 snapshots of biochem process
//		model		= new TopLevelModel(tend/(double)1440);
		qmodel		= new QueueModel(tend/(double)1440);
		
//		listener	= new CellEventListener(false);
		listener	= new QModelListener(false);
		
//		sim			= new Simulator<vec>(model);
		sim			= new Simulator<QCell>(qmodel);
		
		sim->addEventListener(listener);
		
		// Run the simulation
		int output_count = 0;
		BioChem* biochem = BioChem::getInstance();
		double tL = (double)0.0;
		double tLast = (double)-1.0;
		int hours;
		int minutes;
		int seconds;
		
		// Output initial conditions
		qmodel->getBioChemOutput()->output_initial();
		
		tL = sim->nextEventTime();
		
		while (tL < tend){
			
			biochem->advanceTime(tL-biochem->getTime());
			
			sim->execNextEvent();
			
			if (output_count == 0 && tL != tLast){
				tLast = tL;
				hours = tLast/3600;
				minutes = (tLast - hours*3600)/60;
				seconds = tLast - hours*3600 - minutes*60;
				double f_middle = BioChem::getInstance()->getFilterTopPos();
				vec grad_c = BioChem::getInstance()->getChemtacGrad(vec(0.0, 0.0, f_middle));
				double z_grad = grad_c[2];
				cout.precision(3);
				cout << "\r" << 100.0*tL/tend << " %";
				cout << "\r" << hours << ":" << minutes << ":" << seconds 
					 << " " << qmodel->printSizes() << " "
					 << qmodel->getServiceTime(z_grad);
				cout.flush();
			}
			
			output_count = (output_count+1)%1000;
			//exit(-1);
			
			tL = sim->nextEventTime();
		}
		
		cout << "\r" << tL << " " << qmodel->printSizes();
		cout.flush();
		
//		// Output final cell positions
//		listener->cleanup(tL);
		
	} catch(adevs::exception err) {
		cout << endl;
		cerr << "Simulator exception: ";
		cerr << err.what() << endl;
	} catch(...) {
		cout << endl;
		if (errno != 0)
			perror("System error: ");
	}
	cout << endl;
	
	// Cleanup
	if (sim != NULL) delete sim;
//	if (listener != NULL) delete listener;
//	if (model != NULL) delete model;
	if (qmodel != NULL) delete qmodel;
	return 0;
}
