#ifndef _BioChemConstGrad_h_
#define _BioChemConstGrad_h_
#include "BioChemInstance.h"

/**
 * This class implements a constant gradient
 * model derived from c(x)yeah = a x, with 'a' constant.
*/
class BioChemConstGrad:
	public BioChemInstance
{
	public:
	/**
	  * Get the haptoattractant gradient at position p.
	  * */
	vec getHaptacGrad(vec p) const { return vec(0.0); }
	/**
	  * Get the chemotattractant gradient at position p
	  * */
	vec getChemtacGrad(vec p) const { return vec(a); }
	/**
	  * Get the haptoattractant concentration at position p
	  * */
	double getHaptac(vec p) const { return 0.0; }
	/**
	  * Get the chemoattractant concentration at position p
	  * */
	double getChemtac(vec p) const { return p[0]*a; }
	/**
	  * Get the time at which the hapto- and chemo- attractant
	  * quantities are valid for (i.e., the last update time)
	  * */
	double getTime() const { return t; }
	/**
	  * Get the extent of the solution domain in each dimension
	  * */
	vec getDomain() const { return L; }
	/**
	 * Get the 'natural' grid spacing for producing output.
	 */
	double getGridSpacing() const { return dx; }
	/**
	  * Advance the solution by dt units of time.
	  * */
	void advanceTime(double dt) { t += dt; }
	
	/**
	 * Print out the header for VisIt .vtk files
	 */
	std::string VTKHeader() const{
		std::ostringstream outstr;
		outstr << "# vtk DataFile Version 2.0\n" 
		<< "3D Biochemical Diffusion\n" 
		<< "ASCII\n" 
		<< "DATASET STRUCTURED_POINTS\n" 
//		<< "DIMENSIONS " << W << " " << D << " " << H << "\n" 
//		<< "SPACING " << dx << " " << dx << " " << dx << "\n" 
		<< "ORIGIN 0 0 0\n" 
//		<< "POINT_DATA " << W*D*H << "\n" 
		<< "SCALARS volume_scalars double\n" 
		<< "LOOKUP_TABLE default\n";
		return outstr.str();
	}
	
	/**
	 * Print out the concentration at each point on the grid
	 */
	std::string printEnvironment() const{
		std::ostringstream outstr;
		for(int i = 0; i*dx < L[0]; i++){
			outstr << getChemtac(*(new vec(i, 0.0, 0.0))) << "\n";
		}
		return outstr.str();
	}
	
	std::string printFilter() const{
		std::ostringstream outstr;
		outstr << "temp holder";
		return outstr.str();
	}

	/**
	 * not used in the math for this class, but is in the instance 
	 */
	double getDiffCoeff() const{
		return 0;
	}
	
	/**
	 * not used in the math for this class, but is in the instance 
	 */
	double getInitialConc() const{
		return 0;
	}
	
	/**
	  * Virtual destructor.
	  */
	~BioChemConstGrad(){}
	
	/**
	 * Constructor.
	 */
	BioChemConstGrad(double a = -5E-2, double dx = 1.0E-5, double l = 0.02):
		BioChemInstance(), a(a), dx(dx), t(0.0), L(l) {}
	
	private:
	
	double a, dx, t;
	
	vec L;
};

#endif
