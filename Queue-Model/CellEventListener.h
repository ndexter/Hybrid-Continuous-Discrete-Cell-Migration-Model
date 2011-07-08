#ifndef _CellEventListener_h_
#define _CellEventListener_h_
#include "adevs.h"
#include "vec.h"
#include <fstream>

/**
 * This class listens for and records the trajectories
 * of individual cells. Cell trajectories are stored in
 * individual files within the directory "cells"
 */
class CellEventListener:
	public adevs::EventListener<vec>
{
	public:
		/**
		 * Create a listener for cell events. It will
		 * need to be registered with the Simulator to
		 * actually record trajectories. This will throw
		 * an errno code if the directory could not be created.
		 */
		CellEventListener(bool final_only = true);
		/**
		 * This method is activated whenever a model in the
		 * simulation produces an output event.
		 */
		void outputEvent(adevs::Event<vec> x, double t);
		/**
		 * Force an output at time t
		 */
		void cleanup(double t);
		/**
		 * Destructor
		 */
		~CellEventListener();
	private:
		/**
		 * Output stream
		 */
		std::ofstream fout;
		/**
		 * Scratch space for creating file names
		 */
		char scratch[1000];
		/**
		 * Output final positions only?
		 */
		bool final_only;
};

#endif
