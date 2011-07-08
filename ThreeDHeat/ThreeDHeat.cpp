/*
 *  ThreeDHeat.cpp - A combined-method heat simulation.
 *
 *  Created by Nicholas Dexter on 2/2/09.
 *  Copyright 2009. All rights reserved.
 *
 */
#include "ThreeDHeat.h"
#include <iostream>
#include <sys/time.h>
#include <time.h>

using namespace std;

//static const char* outDir = "";
ofstream fout;

int f = 0;

gsl_vector *x;
gsl_vector *b;
gsl_vector *un;

gsl_matrix *LHS;
gsl_matrix *RHS;

gsl_permutation *p;

// No haptoattractant in this version of the model
/*vec ThreeDHeat::getHaptacGrad(vec p) const{
    return vec(0.0);
}*/

// No haptoattractant in this version of the model
/*double ThreeDHeat::getHaptac(vec p) const{
    return 0.0;
}*/

// The chemoattractant solution at the point is obtained by linear
// interpolation on the computational grid.
/*double ThreeDHeat::getChemtac(vec p) const{
    // Make sure the point is in the solution domain
   	assert(p[0] >= 0.0 && p[0] <= W);
	assert(p[1] >= 0.0 && p[1] <= D);
	assert(p[2] >= 0.0 && p[2] <= H);
	
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
}*/

// Approximate the chemoattractant gradient
/*vec ThreeDHeat::getChemtacGrad(vec p) const{
    // Make sure the point is in the domain
	assert(p[0] >= 0.0 && p[0] <= W);
	assert(p[1] >= 0.0 && p[1] <= D);
	assert(p[2] >= 0.0 && p[2] <= H);
	
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
	
	InterPoly px(cdat, tdat, 3);
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
	
	InterPoly py(cdat, tdat, 3);
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
	
	InterPoly pz(cdat, tdat, 3);
	zder =  pz.derivative(p[2]);
	
	return vec(xder, yder, zder);
}*/


double ThreeDHeat::getTime() const{
	return t;
}

double ThreeDHeat::getGridSpacing() const{
	return dx;
}

int ThreeDHeat::getIndex(int i, int j, int k) const{
	return (int)(i + j*W + k*W*D);
}

int ThreeDHeat::getCoordinates(int index, int *coords) const{
	int z, y, x, product = 0;
	if(index >= W*D){
		z = index / (W*D);
		product = z * (int)(W*D);
		index -= product;
		product = 0;
	} else {
		z = 0;
	}
	if(index >= W){
		y = index / W;
		product = y * (int)W;
	} else {
		y = 0;
	}
	x = index - product;
	coords[0] = x;
	coords[1] = y;
	coords[2] = z;
	return *coords;
}

void ThreeDHeat::buildLHS(gsl_matrix *LHS) const{
	for(int r = 0; r < W*D*H; r++){
		int coords[3];
		getCoordinates(r, coords);
		int i = coords[0];
		int j = coords[1];
		int k = coords[2];
		
		gsl_matrix_set(LHS, r, getIndex(i,j,k), (double)1.0 + (double)6.0*alpha*theta);
		if(i == 0){
			gsl_matrix_set(LHS, r, getIndex(i+1,j,k), -(double)2.0*alpha*theta);
		} else if(i == W-1){
			gsl_matrix_set(LHS, r, getIndex(i-1,j,k), -(double)2.0*alpha*theta);
		} else {
			gsl_matrix_set(LHS, r, getIndex(i-1,j,k), -(double)1.0*alpha*theta);
			gsl_matrix_set(LHS, r, getIndex(i+1,j,k), -(double)1.0*alpha*theta);
		}
		if(j == 0){
			gsl_matrix_set(LHS, r, getIndex(i,j+1,k), -(double)2.0*alpha*theta);
		} else if(j == D-1){
			gsl_matrix_set(LHS, r, getIndex(i,j-1,k), -(double)2.0*alpha*theta);
		} else {
			gsl_matrix_set(LHS, r, getIndex(i,j-1,k), -(double)1.0*alpha*theta);
			gsl_matrix_set(LHS, r, getIndex(i,j+1,k), -(double)1.0*alpha*theta);
		}
		if(k == 0){
			gsl_matrix_set(LHS, r, getIndex(i,j,k+1), -(double)2.0*alpha*theta);
		} else if(k == H-1){
			gsl_matrix_set(LHS, r, getIndex(i,j,k-1), -(double)2.0*alpha*theta);
		} else {
			gsl_matrix_set(LHS, r, getIndex(i,j,k-1), -(double)1.0*alpha*theta);
			gsl_matrix_set(LHS, r, getIndex(i,j,k+1), -(double)1.0*alpha*theta);
		}
	}
	
	// Want a look at the LHS matrix? Uncomment these lines
	/*
	for(int i = 0; i < W*D*H; i++){
		int coords[3];
		getCoordinates(i, coords);
		int x = coords[0];
		int y = coords[1];
		int z = coords[2];
		
		fout << "( " << x <<  ", " << y << ", " << z << ")   ";
		
		for(int j = 0; j < W*D*H; j++){
			if(gsl_matrix_get(LHS, i, j) == 0){
				fout << " 0        ";
			} else if (i == j){
				fout << " " << gsl_matrix_get(LHS, i, j) << "      "; 
			} else {
				fout << gsl_matrix_get(LHS, i, j) << " ";
			}
			if(j == W*D*H-1){
				fout << endl << endl << endl << endl;
			}
		}
	}
	*/

	int s;

	gsl_linalg_LU_decomp (LHS, p, &s);
	
}

void ThreeDHeat::buildRHS(gsl_matrix *RHS) const{
	for(int r = 0; r < W*D*H; r++){
		int coords[3];
		getCoordinates(r, coords);
		int i = coords[0];
		int j = coords[1];
		int k = coords[2];
		
		gsl_matrix_set(RHS, r, getIndex(i,j,k), (double)1.0 - (double)6.0*alpha*((double)1.0 - theta));
		if(i == 0){
			gsl_matrix_set(RHS, r, getIndex(i+1,j,k), ((double)2.0*alpha*((double)1.0 - theta)));
			//															+ ((double)4.0*alpha*dx*flux0*((double)1.0 - theta)/k)));
		} else if(i == W-1){
			gsl_matrix_set(RHS, r, getIndex(i-1,j,k), ((double)2.0*alpha*((double)1.0 - theta)));
			//															+ ((double)4.0*alpha*dx*fluxL*((double)1.0 - theta)/k)));
		} else {
			gsl_matrix_set(RHS, r, getIndex(i-1,j,k), alpha*((double)1.0 - theta));
			gsl_matrix_set(RHS, r, getIndex(i+1,j,k), alpha*((double)1.0 - theta));
		}
		if(j == 0){
			gsl_matrix_set(RHS, r, getIndex(i,j+1,k), ((double)2.0*alpha*((double)1.0 - theta)));
			//															+ ((double)4.0*alpha*dx*flux0*((double)1.0 - theta)/k)));
		} else if(j == D-1){
			gsl_matrix_set(RHS, r, getIndex(i,j-1,k), ((double)2.0*alpha*((double)1.0 - theta)));
			//															+ ((double)4.0*alpha*dx*fluxL*((double)1.0 - theta)/k)));
		} else {
			gsl_matrix_set(RHS, r, getIndex(i,j-1,k), alpha*((double)1.0 - theta));
			gsl_matrix_set(RHS, r, getIndex(i,j+1,k), alpha*((double)1.0 - theta));
		}
		if(k == 0){
			gsl_matrix_set(RHS, r, getIndex(i,j,k+1), ((double)2.0*alpha*((double)1.0 - theta)));
			//															+ ((double)4.0*alpha*dx*flux0*((double)1.0 - theta)/k)));
		} else if(k == H-1){
			gsl_matrix_set(RHS, r, getIndex(i,j,k-1), ((double)2.0*alpha*((double)1.0 - theta)));
			//															+ ((double)4.0*alpha*dx*fluxL*((double)1.0 - theta)/k)));
		} else {
			gsl_matrix_set(RHS, r, getIndex(i,j,k-1), alpha*((double)1.0 - theta));
			gsl_matrix_set(RHS, r, getIndex(i,j,k+1), alpha*((double)1.0 - theta));
		}
	}
	
	// Want a look at the RHS matrix? Uncomment these lines
	/*
	for(int i = 0; i < W*D*H; i++){
		int coords[3];
		getCoordinates(i, coords);
		int x = coords[0];
		int y = coords[1];
		int z = coords[2];
	
		fout << "( " << x <<  ", " << y << ", " << z << ")   ";
		
		for(int j = 0; j < W*D*H; j++){
			if(gsl_matrix_get(RHS, i, j) == 0){
				fout << " 0        ";
			} else if (i == j){
				fout << " " << gsl_matrix_get(RHS, i, j) << "      "; 
			} else {
				fout << gsl_matrix_get(RHS, i, j) << " ";
			}
			if(j == W*D*H-1){
				fout << endl << endl << endl << endl;
			}
		}
	}
	*/
}

void ThreeDHeat::advanceTime(double dt_req){
	
//	struct timeval tv;
//	
//	time_t begin, end, cputime;
//		
//	gettimeofday(&tv, NULL); 
//		
//	begin=tv.tv_usec;

	char tmp[100];
	sprintf(tmp,"%d.vtk",f);
	fout.open(tmp);

	fout << "# vtk DataFile Version 2.0\n"
		 << "3D Biochemical Diffusion alpha = " << alpha << "\n" 
		 << "ASCII\n"
		 << "DATASET STRUCTURED_POINTS\n"
		 << "DIMENSIONS " << W << " " << D << " " << H << "\n"
		 << "SPACING " << dx << " " << dx << " " << dx << "\n"
		 << "ORIGIN 0 0 0\n"
		 << "POINT_DATA " << W*D*H << "\n"
		 << "SCALARS volume_scalars double\n"
		 << "LOOKUP_TABLE default\n";
	
	double dt_max;
	
	double tend = t + dt_req;
	
	while(t < tend){
		
		if(theta < 0.5){
			
			dt_max = dx*dx/(DC*(double)6.01); // stability requirement plus safety margin
			
			dt = min(dt_max, tend-t);
			
			t += dt;
			
			alpha = DC*dt/(dx*dx);
			
			// Build the LHS matrix
			buildLHS(LHS);
			
			// Build the RHS matrix
			buildRHS(RHS);
		
		} else if (alpha != DC*dt_req/(dx*dx)){
			
			t += dt_req;
			
			alpha = DC*dt_req/(dx*dx);
			
			// Build the LHS matrix
			buildLHS(LHS);
			
			// Build the RHS matrix
			buildRHS(RHS);			
			
		} else {
			
			t += dt_req;
			
		}
				
		// Multiply the RHS matrix by un
		double sum = 0.0;
		
		for(int row = 0; row < W*D*H; row++){
			for(int col = 0; col < W*D*H; col++){ 
				sum += (gsl_matrix_get(RHS, row, col)*gsl_vector_get(un, col));
			}
			gsl_vector_set(b, row, sum);
			sum = 0.0;
		}
		
		// Want a look at b? Uncomment these lines
		/*
		for(int row = 0; row < W*D*H; row++){
			fout << gsl_vector_get(b, row) << endl;
		}
		*/
		
		gsl_linalg_LU_solve (LHS, p, b, x);		

	}
	
//	gsl_vector_set(un, getIndex(1, 3, 4), 1.0);
//	gsl_vector_set(un, getIndex(2, 7, 1), 1.0);
//	gsl_vector_set(un, getIndex(5, 6, 6), 1.0);
//	gsl_vector_set(un, getIndex(0, 4, 2), 1.0);
//	gsl_vector_set(un, getIndex(7, 7, 3), 1.0);
//	gsl_vector_set(un, getIndex(4, 1, 9), 1.0);
	
	// Print out the diffusion data
	for(int k = 0; k < H; k++){
		for(int j = 0; j < D; j++) {
			for(int i = 0; i < W; i++){
//				fout << t << " " 
//					 << (i) << " " << (j) << " " << (k) << " "
//					 << gsl_vector_get(x, getIndex(i,j,k)) << endl;
				fout << gsl_vector_get(x, getIndex(i,j,k)) << "\n";
			}
		}
	}
	
	*un = *x;
	
	f++;
	
	fout.close();
	
//	gettimeofday(&tv, NULL);
//	end=tv.tv_usec;
//	cputime = end - begin;
	
//	cout << "begin: " << begin 
//			<< " end: " << end 
//			<< " cputime: " << cputime << endl;
}

ThreeDHeat::ThreeDHeat(double dt, double dx, double W, double D, double H, 
					   double DC, double flux0, double fluxL, double k, double theta){
	
	char tmp[100];
	sprintf(tmp,"%d.vtk",f);
	fout.open(tmp);

	
	t = (double) 0.0;
	
	this->dt		= dt;
	this->W 		= W;
	this->D 		= D;
	this->H 		= H;
	this->DC        = DC;
	this->dx   	    = dx;
	this->flux0     = flux0;
	this->fluxL     = fluxL;
	this->k         = k;
	this->theta     = theta;
	
	alpha = DC*dt/(dx*dx);
	
	x  = gsl_vector_alloc(W*D*H);
	b  = gsl_vector_alloc(W*D*H);
	un = gsl_vector_alloc(W*D*H);
	
	LHS = gsl_matrix_calloc(W*D*H, W*D*H);
	RHS = gsl_matrix_calloc(W*D*H, W*D*H);
	
	p = gsl_permutation_alloc (W*D*H);
	
	if(!(theta < 0.5)){
		// Build the LHS matrix
		buildLHS(LHS);
		
		// Build the RHS matrix
		buildRHS(RHS);
	}
	
	fout << "# vtk DataFile Version 2.0\n"
 		 << "3D Biochemical Diffusion alpha = " << alpha << "\n" 
		 << "ASCII\n"
		 << "DATASET STRUCTURED_POINTS\n"
		 << "DIMENSIONS " << W << " " << D << " " << H << "\n"
		 << "SPACING " << dx << " " << dx << " " << dx << "\n"
		 << "ORIGIN 0 0 0\n"
		 << "POINT_DATA " << W*D*H << "\n"
		 << "SCALARS volume_scalars double\n"
		 << "LOOKUP_TABLE default\n";
	
	for(int k = 0; k < H/2; k++){
		for(int i = 0; i < W; i++) {
			for(int j = 0; j < D; j++){
				gsl_vector_set(un, getIndex(i,j,k), (double)1.0);
			}
		}
	}
	for(int k = H/2; k < H; k++){
		for(int i = 0; i < W; i++) {
			for(int j = 0; j < D; j++){
				gsl_vector_set(un, getIndex(i,j,k), (double)0.0);
			}
		}
	}
	
//	for(int k = 0; k < H; k++){
//		for(int i = 0; i < W; i++) {
//			for(int j = 0; j < D; j++){
//				gsl_vector_set(un, getIndex(i,j,k), (double)0.0);
//			}
//		}
//	}
	
//	gsl_vector_set(un, getIndex(1, 3, 4), 1.0);
//	gsl_vector_set(un, getIndex(2, 7, 1), 1.0);
//	gsl_vector_set(un, getIndex(5, 6, 6), 1.0);
//	gsl_vector_set(un, getIndex(0, 4, 2), 1.0);
//	gsl_vector_set(un, getIndex(7, 7, 3), 1.0);
//	gsl_vector_set(un, getIndex(4, 1, 9), 1.0);
	
	// Print out the diffusion data
	for(int k = 0; k < H; k++){
		for(int j = 0; j < D; j++) {
			for(int i = 0; i < W; i++){
				fout << gsl_vector_get(un, getIndex(i,j,k)) << endl;
			}
		}
	}
	
	f++;
	
	fout.close();
}

ThreeDHeat::~ThreeDHeat(){
	gsl_vector_free(x);
	gsl_vector_free(b);
	gsl_permutation_free(p);
	gsl_matrix_free(LHS);
	gsl_matrix_free(RHS);
}
