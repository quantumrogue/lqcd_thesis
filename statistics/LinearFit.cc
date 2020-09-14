#include "LinearFit.h"
#include <iostream> //allows to write files
#include <fstream> //allows to read files
#include <vector> //Using vectors instead of normal arrays - range is not know a priori
#include <cmath> //use mathematical operators
using namespace std;
using namespace fit;




//CONSTRUCTOR
//resize vector of points 
LinearFit::LinearFit (int datasize){
	datapoints.resize(datasize);
	for(int i = 0; i < datasize; i++){
		datapoints[i].resize(3);
	}
}

//returns size of vector - number of points to fit
int LinearFit::datasize(){
	return datapoints.size();
}

//returns fit parameters y = a + bx
double* LinearFit::fitparameters(){
	double a,b,delta_a,delta_b;
	double weight, xi, yi;
	double sumweights = 0;
	double wxy = 0;
	double wx = 0;
	double wy = 0;
	double wxx = 0;
	for(int i = 0; i < datapoints.size(); i++){
		xi = datapoints[i][0];
		yi = datapoints[i][1];
		weight = 1.0/(datapoints[i][2]*datapoints[i][2]);
		sumweights += weight;
		wxy += weight*xi*yi;
		wx += weight*xi;
		wy += weight*yi;
		wxx += weight*xi*xi;
	}
	double delta = sumweights*wxx - wx*wx;
	a = (wxx*wy - wx*wxy)/delta;
	delta_a = sqrt(wxx/delta);
	b = (sumweights*wxy - wx*wy)/delta;
	delta_b = sqrt(sumweights/delta);
	double *parameters_point;
	double parameters[4] = {a,delta_a,b,delta_b};
	parameters_point = parameters;
	return parameters_point;
}

//Returns Chisquared of the Linear Fit
double LinearFit::chisquared(){
	double a,b;
	a = fitparameters()[0];
	b = fitparameters()[2];
	double chisquared = 0;
	double weight, xi, yi;
	for(int i = 0; i < datapoints.size(); i++){
		xi = datapoints[i][0];
		yi = datapoints[i][1];
		weight = 1.0/(datapoints[i][2]*datapoints[i][2]);
		chisquared += weight*(yi - a - b*xi)*(yi - a - b*xi);
	}
	return chisquared;
}	