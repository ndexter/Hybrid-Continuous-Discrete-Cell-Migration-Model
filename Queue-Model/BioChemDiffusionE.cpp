#include "BioChemDiffusionE.h"

vec BioChemDiffusionE::getChemtacGrad(vec p) const 
{
	vec grad;
	p[0] = floor(p[0]/dx)*dx ;
	double nextPos = p[0]+dx;
	double prevPos = p[0];
	nextPos = nextPos/dx;
	prevPos = prevPos/dx; 
	int x_2 = int(nextPos);
	int x_1 = int(prevPos);
	grad[0] = (c[current][x_2]-c[current][x_1]) / (dx);
	return grad;
}

double BioChemDiffusionE::getChemtac(vec p) const 
{
	p = p / getGridSpacing();
	return c[current][(int)p[0]];
}

void BioChemDiffusionE::advanceDiffusion(double timestep)
{
	// h is change in time, DifCoef is diffusion coefficient of chemoattractant
	double r=(DifCoef*timestep/(getGridSpacing()*getGridSpacing()));//<=.5 for system to be stable
	// Set the boundary conditions
	c[current][0] = c[current][1];
	c[current][ArraySize-1] = c[current][ArraySize-2];
	// Compute value at each point in space
	for (int j=1; j < ArraySize-1; j++)
	{
		c[next][j]=(1.0-2.0*r)*c[current][j] + r*c[current][j-1] + r*c[current][j+1];
	}
	current = next;
	next = (next+1)%2;
}

void BioChemDiffusionE::advanceTime(double dt)
{
	int steps = (int)(dt/h);
	for (int i = 0; i < steps; i++)
		advanceDiffusion(h);
	advanceDiffusion(dt-steps*h);
	t += dt;
}

std::string BioChemDiffusionE::VTKHeader() const{
	std::ostringstream outstr;
	outstr << "# vtk DataFile Version 2.0\n" 
	<< "3D Biochemical Diffusion\n" 
	<< "ASCII\n" 
	<< "DATASET STRUCTURED_POINTS\n" 
//	<< "DIMENSIONS " << W << " " << D << " " << H << "\n" 
//	<< "SPACING " << dx << " " << dx << " " << dx << "\n" 
	<< "ORIGIN 0 0 0\n" 
//	<< "POINT_DATA " << W*D*H << "\n" 
	<< "SCALARS volume_scalars double\n" 
	<< "LOOKUP_TABLE default\n";
	return outstr.str();
}

std::string BioChemDiffusionE::printEnvironment() const{
	std::ostringstream outstr;
	outstr << " ";
	return outstr.str();
}

std::string BioChemDiffusionE::printFilter() const{
	std::ostringstream outstr;
	outstr << " ";
	return outstr.str();
}

/*
 * BioChemDiffusionE::~BioChemDiffusionE(){}
 * */
BioChemDiffusionE::BioChemDiffusionE(double a, double dx, vec L, double DifCoef):
	BioChemInstance(),
	a(a),
	dx(dx),
	t(0.0),
	DifCoef(DifCoef),
	L(L) 
	{
		// makes an r value a maximum of .45 which is slightly less than the value 
		// it will take to become unstable = 0.5
		h = 4.0*getGridSpacing()*getGridSpacing()/(10.0*DifCoef);
		ArraySize = int (L[0] / getGridSpacing());
		current = 0;
		next = 1;
		c[0] = new double[ArraySize];
		c[1] = new double[ArraySize];
		c[current][1] = c[current][ArraySize-2] = a;
	}

