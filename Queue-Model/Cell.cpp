#include "Cell.h"
#include "CellInterface.h"
#include "CellPopulation.h"
#include "BioChem.h"
#include <iostream>

using namespace std;
using namespace adevs;


/* Create a randomly oriented unit vector */
static vec rand_unit_vec(){
	static rv r(543534543);
	vec result;
	result[0] = r.uniform(-1.0,1.0);
	result[1] = r.uniform(-1.0,1.0);
	result[2] = r.uniform(-1.0,1.0);
	return result/result.norm2();
}

Cell::Cell(vec p0):Atomic<vec>(),CellInterface(){
	gain = (double)1.66666E8;
	p = p0;
	tstart = 0.0;
	/* Not moving initially */
	v = (double)0.0;
}

Cell::Cell(vec p0, double t0):Atomic<vec>(),CellInterface(){
	gain = (double)1.66666E8;
	p = p0;
	tstart = t0;
	/* Not moving initially */
	v = (double)0.0;
}

double Cell::getSpd(){

//	BioChem* env = BioChem::getInstance();
//	double conc = env->getChemtac(p);
/* This is the variable speed model developed by Matt et al over the summer */
/*	spd = 14.1285* (conc) / (16.5833-0.739225*conc+ 0.0264594*conc*conc);
	spd = spd * 1000;  //convert from micrometers to millimeters
*/
//	spd = 1E-3;  // millimeters per hour
//	spd = (double)17.5;  // microns/hour, interpolated from DiMilla 1993
	spd = (double) 4.86111111E-3; // microns/second
	return spd;
}

double Cell::getF(){
//	BioChem* env = BioChem::getInstance();
//	double conc = env->getChemtac(p);
/* This is the variable persistance time model developed by Matt et al over the summer */
///	double Persis = 10.3323*exp( -0.0029333-30.062*conc+conc*conc);
//	double Persis = (double)0.1;
//	double Persis = (double)0.5;   // hours, min P from DiMilla 1992
//	double Persis = (double)5.77;  // hours, interpolated from DiMilla 1993
	// Time is in seconds, frequency in 1/seconds
	double Persis = (double)10;
//	double Persis = (double)1800;
	f = 1/ Persis;
	return f;
}

double Cell::getDistance(vec p1, vec p2){
	return sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) +
				(p1[1]-p2[1])*(p1[1]-p2[1]) +
				(p1[2]-p2[2])*(p1[2]-p2[2]));
}

void Cell::delta_int(){
	/*
	 * Position was updated in the output function.
	 * 
	 * Updating the model state in the output function
	 * is, in general, not a good idea (it confuses large models).
	 * But it sure simplifies the implementation of this model.
	 */
	/* Pick new direction based on current position */
	/* EXPERIMENTAL BEHAVIOR CODE */
	//	if(p[2] > BioChem::getInstance()->getFilterTopPos()){
	//		vec grad_c;
	//		grad_c[0] = 0.0;
	//		grad_c[1] = 0.0;
	//		grad_c[2] = 1.0;
	//		vec rand_dir = rand_unit_vec();
	//		vec v_total = grad_c+rand_dir;
	//		v = getSpd()*v_total/v_total.norm2();
	//        return;
	//	} else {
	vec grad_c = BioChem::getInstance()->getChemtacGrad(p);
	vec rand_dir = rand_unit_vec();
	vec v_total = gain*grad_c+rand_dir;
	v = getSpd()*v_total/v_total.norm2();
	return;
	//	}
}

void Cell::output_func(adevs::Bag<vec>& yb){
/* EXPERIMENTAL BEHAVIOR CODE */
//	// If the cell hasn't attached
//	if(getTimeAttached() == 0.0){
//		// check to see if the cell could attach
//		vec dist = p + ta()*v;
//		
//		if(dist[2] < BioChem::getInstance()->getFilterTopPos()){
//			// see what would happen if the cell did attach
//			dist[2] = BioChem::getInstance()->getFilterTopPos();
//			
//			// look at all the other cells
//			CellPopulation* pop = CellPopulation::getInstance();
//			CellPopulation::iterator iter;
//			bool settled = false;
//			for(iter = pop->begin(); iter != pop->end(); iter++){
//				vec position = (*iter)->getPos();
//				do{
//					// could this cell collide?
//					if(getDistance(position, dist) < 1.0){
//						// avoid collision by stacking
//						vec stackem;
//						stackem[0] = r.uniform(-1.0,1.0);
//						stackem[1] = r.uniform(-1.0,1.0);
//						stackem[2] = 0.0;
//						dist = position + stackem;
//					}
//				} while (!settled);
//			}
//			
//			p += ta()*v;
//			p[2] = BioChem::getInstance()->getFilterTopPos();
//		} else {
//			// set the position to the new position
//			p = dist;
//		}
//	}
	
	
	/* Move in current direction */
	if (getTime() >= getTimeAttached()) p += ta()*v;
//	if (getTime() >= 3.0) p += ta()*v;
	
	// Don't go below 0.0
	if (p[0] < (double)0.3)
		p[0] = (double)0.3;
	if (p[1] < (double)0.3)
		p[1] = (double)0.3;
	
	// Don't go above Width or Depth
	if (p[0] > BioChem::getInstance()->getDomain()[0]-(double)0.3)
		p[0] = BioChem::getInstance()->getDomain()[0]-(double)0.3;
	if (p[1] > BioChem::getInstance()->getDomain()[1]-(double)0.3)
		p[1] = BioChem::getInstance()->getDomain()[1]-(double)0.3;
	
	// Don't go above or below the filter
	if (p[2] < BioChem::getInstance()->getFilterBottomPos())
		p[2] = BioChem::getInstance()->getFilterBottomPos();
	if (p[2] > BioChem::getInstance()->getFilterTopPos())
		p[2] = BioChem::getInstance()->getFilterTopPos();
	/* Output new position */
	yb.insert(p);
}

double Cell::ta(){
	return 1.0/getF();
}

Cell::~Cell(){ }
