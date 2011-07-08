#include "BioChem.h"
#include "BioChem1DFD.h"
#include "BioChem3DFD.h"
#include "BioChemConstGrad.h"
#include "BioChemDiffusionE.h"
#include "adevs.h"
#include <cmath>
using namespace std;
using namespace adevs;

BioChem* BioChem::getInstance(){
	static BioChem inst;
	return &inst;
}

vec BioChem::getHaptacGrad(vec p) const{
	return impl->getHaptacGrad(p);
}

double BioChem::getHaptac(vec p) const{
	return impl->getHaptac(p);
}

double BioChem::getChemtac(vec p) const{
	return impl->getChemtac(p);
}

vec BioChem::getChemtacGrad(vec p) const{
	return impl->getChemtacGrad(p);
}

double BioChem::getTime() const{
	return impl->getTime();
}

double BioChem::getGridSpacing() const{
	return impl->getGridSpacing();
}

double BioChem::getDiffCoeff() const{
	return impl-> getDiffCoeff();
}
double BioChem::getInitialConc() const{
	return impl-> getInitialConc();
}

void BioChem::advanceTime(double dt_req){
	impl->advanceTime(dt_req);
}

BioChem::BioChem(){
	impl = new BioChem1DFD();
}

vec BioChem::getDomain() const{
	return impl->getDomain();
}

double BioChem::getFilterBottomPos() const{
	return impl->getFilterBottomPos();
}

double BioChem::getFilterTopPos() const{
	return impl->getFilterTopPos();
}

std::string BioChem::VTKHeader() const{
	return impl->VTKHeader();
}

std::string BioChem::printEnvironment() const{
	return impl->printEnvironment();
}

std::string BioChem::printFilter() const{
	return impl->printFilter();
}

BioChem::~BioChem(){
	delete impl;
}
