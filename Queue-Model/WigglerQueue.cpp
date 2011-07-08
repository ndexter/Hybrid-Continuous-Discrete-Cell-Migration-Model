/*
 *  WigglerQueue.cpp
 *  CellMigration
 */

#include "WigglerQueue.h"

void WigglerQueue::delta_int(){
	t += ta();
	wigglerlist.pop();
}

void WigglerQueue::delta_ext(double e, const adevs::Bag<QCell>& xb){
	t += e;
	QCell holder;
	if(!xb.empty()){
		Bag<QCell>::iterator iter = xb.begin();
		while(iter != xb.end()){
			holder = *iter;
			holder.Migrate(t);
			double f_middle = BioChem::getInstance()->getFilterTopPos();
			vec grad_c = BioChem::getInstance()->getChemtacGrad(vec(0.0, 0.0, f_middle));
			double z_grad = grad_c[2];
			double servicetime = serviceTime(z_grad);
			holder.stopMigratingTime(servicetime + t);
			wigglerlist.push(holder);
			iter++;
		}
	}
}

int WigglerQueue::getListSize() const{
	return wigglerlist.size();
}

void WigglerQueue::delta_conf(const adevs::Bag<QCell>& xb){
	delta_int();
	delta_ext(0.0, xb);
}

void WigglerQueue::output_func(adevs::Bag<QCell>& yb){
		yb.insert(wigglerlist.top());
}

double WigglerQueue::ta(){
	if(!wigglerlist.empty()){
		QCell sample = wigglerlist.top();
		double timeadv = wigglerlist.top().getstoptime() - t;
		return timeadv;
	} else {
		return DBL_MAX;
	}
}

double WigglerQueue::serviceTime(double z_grad){
	static rv r(543534543);
	vec randvec;
	randvec[0] = r.uniform(0.0, 1.0);
	randvec[1] = r.uniform(0.0, 1.0);
	randvec[2] = r.uniform(0.0, 1.0);
	randvec = randvec/randvec.norm2();
	
	double retval =  (double)10.0 / 
					((double)4.861E-3*((double)1.6E8 * z_grad + randvec[2]));
	
	return retval;
}