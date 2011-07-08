/*
 *  BioChem1DFD.cpp - A combined-method diffusion simulation.
 *
 *  Created by Nicholas Dexter on 2/2/09.
 *  Copyright 2009. All rights reserved.
 *
 */
#include "BioChem1DFD.h"
#include <iostream>
#include <sys/time.h>
#include <time.h>

using namespace std;

// No haptoattractant in this version of the model
vec BioChem1DFD::getHaptacGrad(vec p) const{
	return vec(0.0);
}

// No haptoattractant in this version of the model
double BioChem1DFD::getHaptac(vec p) const{
	return 0.0;
}

// The chemoattractant solution at the point is obtained by linear
// interpolation on the computational grid.
double BioChem1DFD::getChemtac(vec p) const{
	
	//	cout << p[0] << " " << p[1] << " " << p[2] << endl;
	
	// Make sure the point is in the solution domain
	assert(p[0] >= 0.0 && p[0] <= Width);
	assert(p[1] >= 0.0 && p[1] <= Depth);
	assert(p[2] >= 0.0 && p[2] <= Height);
	
	double flxd = floor(p[0]/dx)*dx;
	double flyd = floor(p[1]/dx)*dx;
	double flzd = floor(p[2]/dx)*dx;
	
	int flz = floor(p[2]/dx);	
	
	int clz = ceil(p[2]/dx);
	
	double xd = p[0] - flxd;
	double yd = p[1] - flyd;
	double zd = p[2] - flzd;
	
	double un_flz = gsl_vector_get(un, flz);
	double un_clz = gsl_vector_get(un, clz);
	
	double i1 = un_flz*((double)1.0 - zd) + un_clz*zd;
	double i2 = un_flz*((double)1.0 - zd) + un_clz*zd;	
	double j1 = un_flz*((double)1.0 - zd) + un_clz*zd;	
	double j2 = un_flz*((double)1.0 - zd) + un_clz*zd;	
	
	double w1 = i1*(1-yd) + i2*yd;
	double w2 = j1*(1-yd) + j2*yd;
	
	double cons = w1*(1-xd) + w2*xd;
	
	return cons;
}

// Approximate the chemoattractant gradient
vec BioChem1DFD::getChemtacGrad(vec p) const{
	// Make sure the point is in the domain
	//	cout << p[0] << " " << p[1] << " " << p[2] << endl;
	assert(p[0] >= 0.0 && p[0] <= Width);
	assert(p[1] >= 0.0 && p[1] <= Depth);
	assert(p[2] >= 0.0 && p[2] <= Height);
	
	double cdat[3], tdat[3], xder, yder, zder;
	
	double flxd = floor(p[0]/dx)*dx;
	double flyd = floor(p[1]/dx)*dx;
	double flzd = floor(p[2]/dx)*dx;
	
	double clxd = ceil(p[0]/dx)*dx;
	double clyd = ceil(p[1]/dx)*dx;
	double clzd = ceil(p[2]/dx)*dx;
	
	//Function value at point vec p, is used in all InterPoly calls
	cdat[0] = getChemtac(vec(flxd, p[1], p[2]));
	cdat[1] = getChemtac(vec(p[0], p[1], p[2]));
	cdat[2] = getChemtac(vec(clxd, p[1], p[2]));
	
	tdat[0] = flxd;
	tdat[1] = p[0];
	tdat[2] = clxd;
	
	//fout << "getChemtacGrad 0 tdat: " 
	//		<< tdat[0] << " " << tdat[1] << " " << tdat[2] << endl;
	//fout << "getChemtacGrad 0 cdat: " 
	//		<< cdat[0] << " " << cdat[1] << " " << cdat[2] << endl;
	
	if (tdat[0] == tdat[1] ||  tdat[1] == tdat[2]) {
		tdat[0] = p[0] - dx;
		tdat[2] = p[0] + dx;
	} 
	
	adevs::InterPoly px(cdat, tdat, 3);
	xder =  px.derivative(p[0]);
	
	cdat[0] = getChemtac (vec(p[0], flyd, p[2]));
	cdat[2] = getChemtac (vec(p[0], clyd, p[2]));
	
	tdat[0] = flyd;
	tdat[1] = p[1];
	tdat[2] = clyd;
	
	if (tdat[0] == tdat[1] || tdat[1] == tdat[2]) {
		tdat[0] = p[1] - dx;
		tdat[2] = p[1] + dx;
	}
	
	adevs::InterPoly py(cdat, tdat, 3);
	yder =  py.derivative(p[1]);
	
	cdat[0] = getChemtac(vec(p[0],p[1], flzd));
	cdat[2] = getChemtac(vec(p[0],p[1], clzd));
	
	tdat[0] = flzd;
	tdat[1] = p[2];
	tdat[2] = clzd;
	
	if (tdat[0] == tdat[1] || tdat[1] == tdat[2]) {
		tdat[0] = p[2] - dx;
		tdat[2] = p[2] + dx;
	} 
	
	adevs::InterPoly pz(cdat, tdat, 3);
	zder =  pz.derivative(p[2]);
	
	return vec(xder, yder, zder);
}

void BioChem1DFD::buildLHS(gsl_vector *lhsdiag, gsl_vector *lhsabovediag, gsl_vector *lhsbelowdiag){

	// Hold the value of alpha incase it needs to change
	double holder = alpha;
	
	// Set the diagonals
	for(int i = 0; i < Height; i++){
		if(i >= getFilterBottomPos() && i <= getFilterBottomPos()){
			alpha = DCf*dt/(dx*dx);
		}
		gsl_vector_set(lhsdiag, i, (double)1.0 + (double)2.0*alpha*theta);
		// Restore alpha
		alpha = holder;
	}
	
	// No flux conditions du/dt = 0
	gsl_vector_set(lhsabovediag, 0, -(double)2.0*alpha*theta);
	//					+ ((double)4.0*alpha*dx*flux0/k)*((double)1.0 - theta)));
	gsl_vector_set(lhsbelowdiag, Height-2, -(double)2.0*alpha*theta);
	//					+ ((double)4.0*alpha*dx*flux0/k)*((double)1.0 - theta)));
	
	// Set the upper and lower diagonals -- Matches book
	for(int i = 1; i < Height-1; i++){
		if(i >= getFilterBottomPos() && i <= getFilterBottomPos()){
			alpha = DCf*dt/(dx*dx);
		}
		gsl_vector_set(lhsabovediag, i, -(double)1.0*alpha*theta);
		// Restore alpha
		alpha = holder;
	}
	for(int i = 0; i < Height-2; i++){
		if(i+1 >= getFilterBottomPos() && i+1 <= getFilterBottomPos()){
			alpha = DCf*dt/(dx*dx);
		}
		gsl_vector_set(lhsbelowdiag, i, -(double)1.0*alpha*theta);
		// Restore alpha
		alpha = holder;
	}
	
}

void BioChem1DFD::buildRHS(gsl_vector *rhsdiag, gsl_vector *rhsabovediag, gsl_vector *rhsbelowdiag) {
	
	// Hold the value of alpha incase it needs to change
	double holder = alpha;
	
	// Set the diagonals
	for(int i = 0; i < Height; i++){
		if(i >= getFilterBottomPos() && i <= getFilterBottomPos()){
			alpha = DCf*dt/(dx*dx);
		}
		gsl_vector_set(rhsdiag, i, (double)1.0 - (double)2.0*alpha*((double)1.0 - theta));
		// Restore alpha
		alpha = holder;
	}
	
	// No flux conditions du/dt = 0
	gsl_vector_set(rhsabovediag, 0, ((double)2.0*alpha*((double)1.0 - theta)));
	//					+ ((double)4.0*alpha*dx*flux0/k)*((double)1.0 - theta)));
	gsl_vector_set(rhsbelowdiag, Height-2, ((double)2.0*alpha*((double)1.0 - theta)));
	//					+ ((double)4.0*alpha*dx*flux0/k)*((double)1.0 - theta)));
	
	// Set the upper and lower diagonals
	for(int i = 1; i < Height-1; i++){
		if(i >= getFilterBottomPos() && i <= getFilterBottomPos()){
			alpha = DCf*dt/(dx*dx);
		}
		gsl_vector_set(rhsabovediag, i, alpha*((double)1.0 - theta));
		// Restore alpha
		alpha = holder;
	}
	for(int i = 0; i < Height-2; i++){
		if(i+1 >= getFilterBottomPos() && i+1 <= getFilterBottomPos()){
			alpha = DCf*dt/(dx*dx);
		}
		gsl_vector_set(rhsbelowdiag, i, alpha*((double)1.0 - theta));
		// Restore alpha
		alpha = holder;
	}
	
	// Set alphaprev so we don't have to do this every time step
	alphaprev = alpha;
}

vec BioChem1DFD::getDomain() const{
	return vec(Width-1,Depth-1,Height-1);
}

std::string BioChem1DFD::VTKHeader() const{
	ostringstream outstr;
	outstr << "# vtk DataFile Version 2.0\n" 
	<< "3D Biochemical Diffusion @ time = " << t << " and alpha = " << alpha << "\n" 
	<< "ASCII\n" 
	<< "DATASET STRUCTURED_POINTS\n" 
	<< "DIMENSIONS " << Width << " " << Depth << " " << Height << "\n" 
	<< "SPACING " << dx << " " << dx << " " << dx << "\n" 
	<< "ORIGIN 0 0 0\n" 
	<< "POINT_DATA " << Width*Depth*Height << "\n" 
	<< "SCALARS volume_scalars double\n" 
	<< "LOOKUP_TABLE default\n";
	return outstr.str();
}

std::string BioChem1DFD::printEnvironment() const{
	ostringstream outstr;
	// Print out the diffusion data
	for(int k = 0; k < Height; k++){
		for(int j = 0; j < Depth; j++) {
			for(int i = 0; i < Width; i++){
				outstr << gsl_vector_get(un, k) << "\n";
			}
		}
	}
	return outstr.str();
}

std::string BioChem1DFD::printFilter() const{
	ostringstream outstr;
	outstr << "# vtk DataFile Version 3.0\n"
	<< "10 micron filter\n"
	<< "ASCII\n"
	<< "DATASET POLYDATA\n"
	<< "POINTS 8 float\n"
	<< "0 0 " << getFilterBottomPos() << "\n"
	<< Width-1 << " 0 " << getFilterBottomPos() << "\n"
	<< Width-1 << " " << Depth-1 << " " << getFilterBottomPos() << "\n"
	<< "0 " << Depth-1 << " " << getFilterBottomPos() << "\n"
	<< "0 0 " << getFilterTopPos() << "\n"
	<< Width-1 << " 0 " << getFilterTopPos() << "\n"
	<< Width-1 << " " << Depth-1 << " " << getFilterTopPos() << "\n"
	<< "0 " << Depth-1 << " " << getFilterTopPos() << "\n\n"
	
	<< "POLYGONS 6 30\n"
	<< "4 0 1 2 3\n"
	<< "4 4 5 6 7\n"
	<< "4 0 1 5 4\n"
	<< "4 2 3 7 6\n"
	<< "4 0 4 7 3\n"
	<< "4 1 2 6 5\n"
	
	<< "CELL_DATA 6\n"
	<< "SCALARS volume_scalars int 1\n"
	<< "LOOKUP_TABLE default\n"
	<< "0\n1\n2\n3\n4\n5\n"
	<< "NORMALS cell_normals float\n"
	<< "0 0 -1\n0 0 1\n0 -1 0\n0 1 0\n-1 0 0\n1 0 0\n"
	<< "FIELD FieldData 2\n"
	<< "cellIds 1 6 int\n"
	<< "0 1 2 3 4 5\n"
	<< "faceAttributes 2 6 float\n"
	<< "0.0 1.0 1.0 2.0 2.0 3.0 3.0 4.0 4.0 5.0 5.0 6.0\n"
	<< "POINT_DATA 8\n"
	<< "SCALARS sample_scalars float 1\n"
	<< "LOOKUP_TABLE my_table\n"
	<< "0.0\n1.0\n2.0\n3.0\n4.0\n5.0\n6.0\n7.0\n"
	<< "LOOKUP_TABLE my_table 8\n"
	<< "0.0 0.0 0.0 1.0\n"
	<< "1.0 0.0 0.0 1.0\n"
	<< "0.0 1.0 0.0 1.0\n"
	<< "1.0 1.0 0.0 1.0\n"
	<< "0.0 0.0 1.0 1.0\n"
	<< "1.0 0.0 1.0 1.0\n"
	<< "0.0 1.0 1.0 1.0\n"
	<< "1.0 1.0 1.0 1.0\n";
	return outstr.str();
}


void BioChem1DFD::advanceTime(double dt_req){
	
	double dt_max;
	
	double tend = t + dt_req;
	
	dt = (double) dt_req; // 1 second time step default
	
	alpha = DC*dt/(dx*dx);	
	
	if(tend != t){
		while(t < tend){
			
			if(theta < 0.5){
				
				dt_max = dx*dx/(DC*(double)6.01); // stability requirement plus safety margin
				
				dt = min(dt_max, tend-t);
				
				if(alpha != DC*dt/(dx*dx)){
					alpha = DC*dt/(dx*dx);
					
					// Build the LHS matrix
					buildLHS(lhsdiag, lhsabovediag, lhsbelowdiag);
					
					// Build the RHS matrix
					buildRHS(rhsdiag, rhsabovediag, rhsbelowdiag);
				}
				
			} else if(alpha != alphaprev){
				
				// Build the LHS matrix
				buildLHS(lhsdiag, lhsabovediag, lhsbelowdiag);
				
				// Build the RHS matrix
				buildRHS(rhsdiag, rhsabovediag, rhsbelowdiag);
				
			} else if(tend <= t + dt){
				
				dt = tend - t;
				
				alpha = DC*dt/(dx*dx);
				
				// Build the LHS matrix
				buildLHS(lhsdiag, lhsabovediag, lhsbelowdiag);
				
				// Build the RHS matrix
				buildRHS(rhsdiag, rhsabovediag, rhsbelowdiag);
			}
			
			t += dt;
			
			// Tridiagonal RHS matrix by vector un 
			
			// boundary condition
			gsl_vector_set(b, 0, 
						   gsl_vector_get(rhsdiag, 0)*gsl_vector_get(un, 0) + 
						   gsl_vector_get(rhsabovediag, 0)*gsl_vector_get(un, 1));
			// central difference
			for(int i = 1; i < Height-1; i++){
				gsl_vector_set(b, i, 
							   gsl_vector_get(rhsbelowdiag, i - 1)*gsl_vector_get(un, i - 1) +
							   gsl_vector_get(rhsdiag, i)*gsl_vector_get(un, i) + 
							   gsl_vector_get(rhsabovediag, i)*gsl_vector_get(un, i + 1));
			}	 
			// boundary condition
			gsl_vector_set(b, Height-1, 
						   gsl_vector_get(rhsbelowdiag, Height-2)*gsl_vector_get(un, Height-2) + 
						   gsl_vector_get(rhsdiag, Height-1)*gsl_vector_get(un, Height-1));
			
			// Solve tridiagonal system A*x = b where x = un+1 and b = RHS*un
			gsl_linalg_solve_tridiag(lhsdiag, 
									 lhsabovediag, 
									 lhsbelowdiag, 
									 b, 
									 x);		
		}
		*un = *x;
	}
}

BioChem1DFD::BioChem1DFD(){
	
	t	    = (double) 0.0;						// seconds
	
	dx	    = (double) 1.0;						// microns
	Width   = (double) 3000.0;					// microns
	Depth   = (double) 3000.0;					// microns
	Height  = (double) 3510.0;					// microns
	ftop	= (double) Height/2 + (double)5.0;	// microns
	fbottom = (double) Height/2 - (double)5.0;	// microns
	DC	    = (double) 2.7;						// microns^2/second
	DCf	    = (double) 0.0148;					// microns^2/second
	theta   = (double) 0.5;						// combined method coefficient
	flux0   = (double) 0.0;						// flux at 0
	fluxL   = (double) 0.0;						// flux at L
	k	    = (double) 1.0;						// 1/DC
	c0	    = (double) 2.0E-5;					// femtograms/micron^3
	
	alpha = -1;
	
	alphaprev = 0;
	
	x =  gsl_vector_alloc(Height);
	b =  gsl_vector_alloc(Height);
	un = gsl_vector_alloc(Height);
	
	lhsdiag = gsl_vector_alloc(Height);
	rhsdiag = gsl_vector_alloc(Height);
	lhsabovediag = gsl_vector_alloc(Height-1);
	rhsabovediag = gsl_vector_alloc(Height-1);
	lhsbelowdiag = gsl_vector_alloc(Height-1);
	rhsbelowdiag = gsl_vector_alloc(Height-1);
	
	for(int i = 0; i < fbottom; i++) {
			gsl_vector_set(un, i, (double) c0);
	}
	for(int i = fbottom; i < Height; i++) {
			gsl_vector_set(un, i, (double) 0.0);
	}
	
}

BioChem1DFD::~BioChem1DFD(){
	gsl_vector_free(x);
	gsl_vector_free(b);
	gsl_vector_free(lhsdiag);
	gsl_vector_free(lhsabovediag);
	gsl_vector_free(lhsbelowdiag);
	gsl_vector_free(rhsdiag);
	gsl_vector_free(rhsabovediag);
	gsl_vector_free(rhsbelowdiag);
}