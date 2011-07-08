#ifndef _BioChemOutput_h_
#define _BioChemOutput_h_
#include <fstream>
#include "adevs.h"
#include "vec.h"
/**
 * Model for periodically generating biochemistry
 * output.
 */
class BioChemOutput:
public adevs::Atomic<vec>
{
public:
	/**
	 * Create an atomic model to record the biochemistry
	 * output at regular intervals. This makes plotting
	 * the data much easier. This model must be part of the
	 * top level simulation model in order to do its work.
	 * The object creates a directory "biochem" to store
	 * its output and throws an errno code of there is
	 * a problem creating the directory.
	 */
	BioChemOutput(double t_rec);
	/*
	 * DEVS methods. Do nothing - all data is recorded in the
	 * output function where the biochemistry is updated as well.
	 */
	void delta_int(){}
	void delta_ext(double, const adevs::Bag<vec>&){}
	void delta_conf(const adevs::Bag<vec>&){}
	/* No garbage collection needed */
	void gc_output(adevs::Bag<vec>&){}
	
	int getFnum() { return fnum; }
	
	/**
	 * Output initial conditions in the domain
	 */
	void output_initial();
	/**
	 * Output function that updates and records the bio chemistry data
	 */
	void output_func(adevs::Bag<vec>&);
	/**
	 * Time advance schedules periodic data dumps
	 */
	double ta();
	/**
	 * Destructor
	 */
	~BioChemOutput();
private:
	/**
	 * Output file for recording the biochemistry process
	 */
	std::ofstream fout;
	std::ofstream fout_grad;
	
	/**
	 * file number
	 */
	int fnum;
	
	/**
	 * Output recording interval
	 */
	double dt;
};

#endif
