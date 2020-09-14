#include <vector>
#include <iostream> //allows to write files
#include <fstream> //allows to read files
using namespace std;
#ifndef LinearFit_h
#define LinearFit_h

namespace fit
{
	class LinearFit
	{
	public:
		//y = a + bx
		//class Members;
		//datapoints = < ... <xi,yi,delta yi> ... >
		vector<vector<double> > datapoints;
		LinearFit (int datasize);

		int datasize();

		double* fitparameters();
		double chisquared();
	};
}

#endif