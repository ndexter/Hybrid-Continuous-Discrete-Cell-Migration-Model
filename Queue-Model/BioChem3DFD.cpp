/*
 *  BioChem3DFD.cpp - A combined-method diffusion simulation.
 *
 *  Created by Nicholas Dexter on 2/2/09.
 *  Copyright 2009. All rights reserved.
 *
 */
#include "BioChem3DFD.h"
#include <iostream>
#include <sys/time.h>
#include <time.h>

using namespace std;

//static const char* outDir = "";
//ofstream fout;

// No haptoattractant in this version of the model
vec BioChem3DFD::getHaptacGrad(vec p) const{
	return vec(0.0);
}

// No haptoattractant in this version of the model
double BioChem3DFD::getHaptac(vec p) const{
	return 0.0;
}

// The chemoattractant solution at the point is obtained by linear
// interpolation on the computational grid.
double BioChem3DFD::getChemtac(vec p) const{
	
//	cout << p[0] << " " << p[1] << " " << p[2] << endl;

	// Make sure the point is in the solution domain
	assert(p[0] >= 0.0 && p[0] <= Width);
	assert(p[1] >= 0.0 && p[1] <= Depth);
	assert(p[2] >= 0.0 && p[2] <= Height);
	
	double flxd = floor(p[0]/dx)*dx;
	double flyd = floor(p[1]/dx)*dx;
	double flzd = floor(p[2]/dx)*dx;
	
	int flx = floor(p[0]/dx);
	int fly = floor(p[1]/dx);
	int flz = floor(p[2]/dx);	
	
	int clx = ceil(p[0]/dx);
	int cly = ceil(p[1]/dx);
	int clz = ceil(p[2]/dx);
	
	double xd = p[0] - flxd;
	double yd = p[1] - flyd;
	double zd = p[2] - flzd;
	
	double i1 = gsl_vector_get(un, getIndex(flx, fly, flz))*((double)1.0 - zd) + 
		gsl_vector_get(un, getIndex(flx, fly, clz))*zd;
	double i2 = gsl_vector_get(un, getIndex(flx, cly, flz))*((double)1.0 - zd) + 
		gsl_vector_get(un, getIndex(flx, cly, clz))*zd;	
	double j1 = gsl_vector_get(un, getIndex(clx, fly, flz))*((double)1.0 - zd) + 
		gsl_vector_get(un, getIndex(clx, fly, clz))*zd;	
	double j2 = gsl_vector_get(un, getIndex(clx, cly, flz))*((double)1.0 - zd) + 
		gsl_vector_get(un, getIndex(clx, cly, clz))*zd;	
	
	double w1 = i1*(1-yd) + i2*yd;
	double w2 = j1*(1-yd) + j2*yd;
	
	double cons = w1*(1-xd) + w2*xd;
	
	return cons;
}

// Approximate the chemoattractant gradient
vec BioChem3DFD::getChemtacGrad(vec p) const{
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

int BioChem3DFD::getIndex(int i, int j, int k) const{
	return (int)(i + j*Width + k*Width*Depth);
}

int BioChem3DFD::getCoordinates(int index, int *coords) const{
	int z, y, x, product = 0;
	if(index >= Width*Depth){
		z = index / (Width*Depth);
		product = z * (int)(Width*Depth);
		index -= product;
		product = 0;
	} else {
		z = 0;
	}
	if(index >= Width){
		y = index / Width;
		product = y * (int)Width;
	} else {
		y = 0;
	}
	x = index - product;
	coords[0] = x;
	coords[1] = y;
	coords[2] = z;
	return *coords;
}

void BioChem3DFD::buildLHS(gsl_matrix *LHS){
	double holder = alpha;
	for(int r = 0; r < Width*Depth*Height; r++){
		int coords[3];
		getCoordinates(r, coords);
		int i = coords[0];
		int j = coords[1];
		int k = coords[2];
		
		if(k < getFilterTopPos() && k > getFilterBottomPos()){
			alpha = DCf*dt/(dx*dx);
		}
		
		gsl_matrix_set(LHS, r, getIndex(i,j,k), (double)1.0 + (double)6.0*alpha*theta);
		if(i == 0){
			gsl_matrix_set(LHS, r, getIndex(i+1,j,k), -(double)2.0*alpha*theta);
		} else if(i == Width-1){
			gsl_matrix_set(LHS, r, getIndex(i-1,j,k), -(double)2.0*alpha*theta);
		} else {
			gsl_matrix_set(LHS, r, getIndex(i-1,j,k), -(double)1.0*alpha*theta);
			gsl_matrix_set(LHS, r, getIndex(i+1,j,k), -(double)1.0*alpha*theta);
		}
		if(j == 0){
			gsl_matrix_set(LHS, r, getIndex(i,j+1,k), -(double)2.0*alpha*theta);
		} else if(j == Depth-1){
			gsl_matrix_set(LHS, r, getIndex(i,j-1,k), -(double)2.0*alpha*theta);
		} else {
			gsl_matrix_set(LHS, r, getIndex(i,j-1,k), -(double)1.0*alpha*theta);
			gsl_matrix_set(LHS, r, getIndex(i,j+1,k), -(double)1.0*alpha*theta);
		}
		if(k == 0){
			gsl_matrix_set(LHS, r, getIndex(i,j,k+1), -(double)2.0*alpha*theta);
		} else if(k == Height-1){
			gsl_matrix_set(LHS, r, getIndex(i,j,k-1), -(double)2.0*alpha*theta);
		} else {
			gsl_matrix_set(LHS, r, getIndex(i,j,k-1), -(double)1.0*alpha*theta);
			gsl_matrix_set(LHS, r, getIndex(i,j,k+1), -(double)1.0*alpha*theta);
		}
		alpha = holder;
	}
	
	// Want a look at the LHS matrix? Uncomment these lines
	/*for(int i = 0; i < Width*Depth*Height; i++){
		int coords[3];
		getCoordinates(i, coords);
		int x = coords[0];
		int y = coords[1];
		int z = coords[2];
		
		fout << "( " << x <<  ", " << y << ", " << z << ")   ";
		
		for(int j = 0; j < Width*Depth*Height; j++){
			if(gsl_matrix_get(LHS, i, j) == 0){
				fout << " 0        ";
			} else if (i == j){
				fout << " " << gsl_matrix_get(LHS, i, j) << "      "; 
			} else {
				fout << gsl_matrix_get(LHS, i, j) << " ";
			}
			if(j == Width*Depth*Height-1){
				fout << endl << endl << endl << endl;
			}
		}
	}*/
	
	int s;
	
	gsl_linalg_LU_decomp (LHS, p, &s);
	
}

void BioChem3DFD::buildRHS(gsl_matrix *RHS){
	double holder = alpha;
	for(int r = 0; r < Width*Depth*Height; r++){
		int coords[3];
		getCoordinates(r, coords);
		int i = coords[0];
		int j = coords[1];
		int k = coords[2];
		
		if(k < getFilterTopPos() && k > getFilterBottomPos()){
			alpha = DCf*dt/(dx*dx);
		}
		
		gsl_matrix_set(RHS, r, getIndex(i,j,k), (double)1.0 - (double)6.0*alpha*((double)1.0 - theta));
		if(i == 0){
			gsl_matrix_set(RHS, r, getIndex(i+1,j,k), ((double)2.0*alpha*((double)1.0 - theta)));
			//															+ ((double)4.0*alpha*dx*flux0*((double)1.0 - theta)/k)));
		} else if(i == Width-1){
			gsl_matrix_set(RHS, r, getIndex(i-1,j,k), ((double)2.0*alpha*((double)1.0 - theta)));
			//															+ ((double)4.0*alpha*dx*fluxL*((double)1.0 - theta)/k)));
		} else {
			gsl_matrix_set(RHS, r, getIndex(i-1,j,k), alpha*((double)1.0 - theta));
			gsl_matrix_set(RHS, r, getIndex(i+1,j,k), alpha*((double)1.0 - theta));
		}
		if(j == 0){
			gsl_matrix_set(RHS, r, getIndex(i,j+1,k), ((double)2.0*alpha*((double)1.0 - theta)));
			//															+ ((double)4.0*alpha*dx*flux0*((double)1.0 - theta)/k)));
		} else if(j == Depth-1){
			gsl_matrix_set(RHS, r, getIndex(i,j-1,k), ((double)2.0*alpha*((double)1.0 - theta)));
			//															+ ((double)4.0*alpha*dx*fluxL*((double)1.0 - theta)/k)));
		} else {
			gsl_matrix_set(RHS, r, getIndex(i,j-1,k), alpha*((double)1.0 - theta));
			gsl_matrix_set(RHS, r, getIndex(i,j+1,k), alpha*((double)1.0 - theta));
		}
		if(k == 0){
			gsl_matrix_set(RHS, r, getIndex(i,j,k+1), ((double)2.0*alpha*((double)1.0 - theta)));
			//															+ ((double)4.0*alpha*dx*flux0*((double)1.0 - theta)/k)));
		} else if(k == Height-1){
			gsl_matrix_set(RHS, r, getIndex(i,j,k-1), ((double)2.0*alpha*((double)1.0 - theta)));
			//															+ ((double)4.0*alpha*dx*fluxL*((double)1.0 - theta)/k)));
		} else {
			gsl_matrix_set(RHS, r, getIndex(i,j,k-1), alpha*((double)1.0 - theta));
			gsl_matrix_set(RHS, r, getIndex(i,j,k+1), alpha*((double)1.0 - theta));
		}
		alpha = holder;
	}
	
	// Want a look at the RHS matrix? Uncomment these lines
	/*for(int i = 0; i < Width*Depth*Height; i++){
		int coords[3];
		getCoordinates(i, coords);
		int x = coords[0];
		int y = coords[1];
		int z = coords[2];
		
		fout << "( " << x <<  ", " << y << ", " << z << ")   ";
		
		for(int j = 0; j < Width*Depth*Height; j++){
			if(gsl_matrix_get(RHS, i, j) == 0){
				fout << " 0        ";
			} else if (i == j){
				fout << " " << gsl_matrix_get(RHS, i, j) << "      "; 
			} else {
				fout << gsl_matrix_get(RHS, i, j) << " ";
			}
			if(j == Width*Depth*Height-1){
				fout << endl << endl << endl << endl;
			}
		}
	}*/
	
	alphaprev = alpha;
}

void BioChem3DFD::advanceTime(double dt_req){
	
	double dt_max;
	
	double tend = t + dt_req;
	
	dt = (double) 1.0; // 1 second time step default
	
	alpha = DC*dt/(dx*dx);
	
	while(t < tend){
		
		if(theta < 0.5){
			
			dt_max = dx*dx/(DC*(double)6.01); // stability requirement plus safety margin
			
			dt = min(dt_max, tend-t);
			
			if(alpha != DC*dt/(dx*dx)){
				alpha = DC*dt/(dx*dx);
				
				// Build the LHS matrix
				buildLHS(LHS);
				
				// Build the RHS matrix
				buildRHS(RHS);
			}
			
		} else if(alpha != alphaprev){
		
			// Build the LHS matrix
			buildLHS(LHS);
			
			// Build the RHS matrix
			buildRHS(RHS);
			
		} else if(dt+t > tend){

			dt = tend - t;
			
			alpha = DC*dt/(dx*dx);
			
			// Build the LHS matrix
			buildLHS(LHS);
			
			// Build the RHS matrix
			buildRHS(RHS);
		}

		
		t += dt;
		
		// Multiply the RHS matrix by un
		double sum = 0.0;
		
		for(int row = 0; row < Width*Depth*Height; row++){
			for(int col = 0; col < Width*Depth*Height; col++){ 
				sum += (gsl_matrix_get(RHS, row, col)*gsl_vector_get(un, col));
			}
			gsl_vector_set(b, row, sum);
			sum = 0.0;
		}
		
		// Want a look at b? Uncomment these lines
		/*for(int row = 0; row < W*D*H; row++){
			fout << gsl_vector_get(b, row) << endl;
		}*/
		 
		
		gsl_linalg_LU_solve (LHS, p, b, x);		
		
	}
	
	//	gsl_vector_set(un, getIndex(1, 3, 4), 1.0);
	//	gsl_vector_set(un, getIndex(2, 7, 1), 1.0);
	//	gsl_vector_set(un, getIndex(5, 6, 6), 1.0);
	//	gsl_vector_set(un, getIndex(0, 4, 2), 1.0);
	//	gsl_vector_set(un, getIndex(7, 7, 3), 1.0);
	//	gsl_vector_set(un, getIndex(4, 1, 9), 1.0);
	
	*un = *x;
	
	//	gettimeofday(&tv, NULL);
	//	end=tv.tv_usec;
	//	cputime = end - begin;
	
	//	cout << "begin: " << begin 
	//			<< " end: " << end 
	//			<< " cputime: " << cputime << endl;
}


vec BioChem3DFD::getDomain() const{
	return vec(Width-1,Depth-1,Height-1);
}

std::string BioChem3DFD::VTKHeader() const{
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

std::string BioChem3DFD::printEnvironment() const{
	ostringstream outstr;
	// Print out the diffusion data
	for(int k = 0; k < Height; k++){
		for(int j = 0; j < Depth; j++) {
			for(int i = 0; i < Width; i++){
				outstr << gsl_vector_get(un, getIndex(i,j,k)) << "\n";
			}
		}
	}
	return outstr.str();
}

std::string BioChem3DFD::printFilter() const{
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

BioChem3DFD::BioChem3DFD(){	
	
	t = (double) 0.0;
	
	dx   	= (double) 1.0;						// microns
	Width 	= (double) 5.0;						// microns
	Depth 	= (double) 5.0;						// microns
	Height	= (double) 50.0;					// microns
	ftop	= (double) Height/2 + (double)5.0;	// microns
	fbottom = (double) Height/2 - (double)5.0;	// microns
	DC      = (double) 2.7;						// microns^2/second
	DCf		= (double) 0.0148;					// microns^2/second
	flux0   = (double) 0.0;						// flux at 0
	fluxL   = (double) 0.0;						// flux at L
	k       = (double) 1.0;						// 1/DC
	theta   = (double) 0.5;						// combined method coefficient
	EX_c0	= (double) 2.0E-5;					// femtograms/micron^3

	alpha = -1;
	
	alphaprev = 0;
	
	x  = gsl_vector_alloc(Width*Depth*Height);
	b  = gsl_vector_alloc(Width*Depth*Height);
	un = gsl_vector_alloc(Width*Depth*Height);
	
	LHS = gsl_matrix_calloc(Width*Depth*Height, Width*Depth*Height);
	RHS = gsl_matrix_calloc(Width*Depth*Height, Width*Depth*Height);
	
	p = gsl_permutation_alloc (Width*Depth*Height);
	
//	for(int k = 0; k < Height/2; k++){
//		for(int i = 0; i < Width; i++) {
//			for(int j = 0; j < Depth; j++){
//				gsl_vector_set(un, getIndex(i,j,k), EX_c0);
//			}
//		}
//	}
//	for(int k = Height/2; k < Height; k++){
//		for(int i = 0; i < Width; i++) {
//			for(int j = 0; j < Depth; j++){
//				gsl_vector_set(un, getIndex(i,j,k), (double)0.0);
//			}
//		}
//	}
	
	for(int k = 0; k < getFilterBottomPos(); k++){
		for(int i = 0; i < Width; i++) {
			for(int j = 0; j < Depth; j++){
				gsl_vector_set(un, getIndex(i,j,k), EX_c0);
			}
		}
	}
	for(int k = getFilterBottomPos(); k < Height; k++){
		for(int i = 0; i < Width; i++) {
			for(int j = 0; j < Depth; j++){
				gsl_vector_set(un, getIndex(i,j,k), (double)0.0);
			}
		}
	}
	
	//	for(int k = 0; k < Height; k++){
	//		for(int i = 0; i < Width; i++) {
	//			for(int j = 0; j < Depth; j++){
	//				gsl_vector_set(un, getIndex(i,j,k), (double)0.0);
	//			}
	//		}
	//	}
	
//	gsl_vector_set(un, getIndex(1, 3, 4), EX_c0);
//	gsl_vector_set(un, getIndex(2, 7, 1), EX_c0);
//	gsl_vector_set(un, getIndex(5, 6, 6), EX_c0);
//	gsl_vector_set(un, getIndex(0, 4, 2), EX_c0);
//	gsl_vector_set(un, getIndex(7, 7, 3), EX_c0);
//	gsl_vector_set(un, getIndex(4, 1, 9), EX_c0);
	
}

BioChem3DFD::~BioChem3DFD(){
	gsl_vector_free(x);
	gsl_vector_free(b);
	gsl_permutation_free(p);
	gsl_matrix_free(LHS);
	gsl_matrix_free(RHS);
}
