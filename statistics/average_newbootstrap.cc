#include <iostream> //allows to write files
#include <fstream> //allows to read files
#include <vector> //Using vectors instead of normal arrays - range is not know a priori
#include <string> //allows to use strings - useful for filenames
#include <sstream>
#include <random> //to generate random numbers
#include <ctime> //using time to generate 'random' seed
#include <algorithm> //to use 'find' function to find an element in a vector
#include <map> //using maps
#include <cmath> //use mathematical operators
#include <iomanip>
using namespace std;

/*
This program reads all raw files and averages its values by a bootstrap procedure;

NECESSARY CHANGES:
nr of files
nr of samples
nr of columns
nr of independent columns
file names template
*/


void srand(unsigned int seed);
//User defined functions
int getfilesize();
void getfilescontent(vector<vector<vector<double> > > &v_filescontent, vector<vector<double> > &v_average, const int Nconfig, const int Ncolumns, const int Nindependent);
void bootstrap(vector<vector<vector<double> > > &v_filescontent , vector<vector<double> > &v_errors, const int Nconfig, const int Nsamples, const int line_number, const int Ncolumns, const int Nindependent);

int main(){

	int Nconfig = 50; //number of configurations - files
	int Nsamples = 100; //Number of samples for the Bootstrap
	int Ncolumns = 22; //Number of columns of information (including momentum!!!)
	int Nindependent = 6; //Number of independent variables (SHOULD BE THE FIRST COLUMNS ALWAYS)

	int line_number; //number of lines in the first file - should be the same for each file
	line_number = getfilesize();

	//Creating vector(line_number) of vectors(Ncolumns)
	//each collumn representes a given physical quantity
	vector<vector<double> > v_average(line_number, vector<double>(Ncolumns));

	//Creating a vector that will save the content of all files (except the independent variables)
	vector<vector<vector<double> > > v_filescontent(Nconfig, vector<vector<double> >(line_number, vector<double>(Ncolumns - Nindependent)));

	//filling the vector v_filescontent with the files columns, and v_average with the corresponding average
	getfilescontent(v_filescontent, v_average, Nconfig, Ncolumns, Nindependent);
	
	//Creating vector(line_number) of vectors(Ncolumns) - same as v_average but with errors
	vector<vector<double> > v_errors(line_number, vector<double>(Ncolumns - Nindependent));

	//bootstrat adds the errors to the vector v_errors
	//the momentum column will be made up of zeros
	bootstrap(v_filescontent, v_errors, Nconfig, Nsamples, line_number, Ncolumns, Nindependent);
	
	
	//Printing information
	for (int line = 0; line < line_number; line++){
		for (int column = 0; column < Ncolumns; column++){
			if (column < Nindependent){
				cout << setprecision(6) << v_average[line][column] << "	";
			}
			else{
				cout << setprecision(15) << v_average[line][column] << "	" << v_errors[line][column - Nindependent] << "	";
			} 
		}
		cout << endl;
	}

	

	return 0;
} 

//counts the number of lines of the first document/configuration (vertex1.dat)
int getfilesize(){
	
	//string firstfilename = "raw_datadiag/gluon_b6.0_80.4-landau-2.dat";
	//string firstfilename = "raw_datagen_imp/gluon-genimp_b6.0_80.4-landau-1.dat";
	string firstfilename = "raw_datagen_naive/gluon_b6.0_80.4-landau-1.dat";
	//Open the first file to count number of lines
	ifstream firstfile;
	int line_number = 0;
	firstfile.open(firstfilename);
	string line;

	if (firstfile.is_open()){
		while (getline(firstfile, line)){
			line_number++;
		}
	}
	else {
		cout << "Could not open first file! - Ending program." << endl;
		return 0;
	}
	firstfile.close();
	//cout << "Number of momentum entries - lines: " << line_number << endl;
	return line_number;
}


//opens each file and adds each column to each entry of the vector v_average, averaging it (except the independent variables)
//saves all the information from the files to v_filesscontent
void getfilescontent(vector<vector<vector<double> > > &v_filescontent, vector<vector<double> > &v_average, const int Nconfig, const int Ncolumns, const int Nindependent){

	//string filename_base = "raw_datadiag/gluon_b6.0_80.4-landau-"; 					//change directory if needed
	//string filename_base = "raw_datagen_imp/gluon-genimp_b6.0_80.4-landau-";
	string filename_base = "raw_datagen_naive/gluon_b6.0_80.4-landau-";
	string extension = ".dat";
	
	//Iterating through each file/configuration
	for (int i = 0; i < Nconfig; i++){

		//converting variable i to string
		string file_number = to_string(i+1);
		//filename of current configuration
		string filename = filename_base + file_number + extension;

		//Open file to read
		ifstream file;
		file.open(filename);
		string line;
		int line_count = 0;
		
		if (file.is_open()){
			while (getline(file, line)){
				stringstream ss(line);
				for (int column = 0; column < Ncolumns; column++){ 			//cycle through each column
					double value;
					ss >> value;
					if (column < Nindependent){						//only adds the independent variables once 
						v_average[line_count][column] = value;
					}
					else{
						v_average[line_count][column] += value/Nconfig; 		//To get the average of all configurations
						v_filescontent[i][line_count][column - Nindependent] = value;
					}	
				}
				line_count++;
			}
			file.close();
		}
	}
}

//adds a third number to each interior vector of v_average with the standard error
//calculated using the bootstrap method
void bootstrap(vector<vector<vector<double> > > &v_filescontent, vector<vector<double> > &v_errors, const int Nconfig, const int Nsamples, const int line_number, const int Ncolumns, const int Nindependent){

	//using a for loop in the columns avoids using higher dimensional vectors (sample_vector)
	for (int column = Nindependent; column < Ncolumns; column++){

		std::mt19937_64 generator(0);
   		std::uniform_int_distribution<int> distribution(0, Nconfig-1);

		//define the vector that will have all samples information -<samples...<files to open>...>
		//each element of this vector will start at zero, for each collumn
		vector<vector<int> > sample_vector(Nsamples, vector<int>(Nconfig));

		//iterating/creating through each sample
		for (int sample = 0; sample < Nsamples; sample++){
			//generate an array of Nconfig random numbers between 0 and Nconfig-1 (correct index in the vector v_filescontent)
			for (int i = 0; i < Nconfig; i++){
				sample_vector[sample][i] = distribution(generator);
			}
		}

		double C = 0.675;

		//Iterating through each line (momentum) to get the corresponding error
		for (int line = 0; line < line_number; line++){
			//Create a vector containing the average of each sample for each line (momentum)
			vector<double> sample_averages(Nsamples);
			double global_sample_average = 0;
			for (int sample = 0; sample < Nsamples; sample++){
				double sample_average = 0;
				for (int entry = 0; entry < Nconfig; entry++){
					int file;
					file = sample_vector[sample][entry]; 	//this is the number of the file for each entry of the sample
					sample_average += v_filescontent[file][line][column - Nindependent]/Nconfig;
				}
				sample_averages[sample] = sample_average;
				global_sample_average += sample_average/Nsamples;
			}
			//Having the average of all samples for each line, we now sort these values by its value
			sort(sample_averages.begin(), sample_averages.end());
			double upper_fraction =  Nsamples*(1 + C)/2;
			double down_fraction = Nsamples*(1 - C)/2;
			int upper_index = (int)(upper_fraction + 0.5);
			int down_index = (int)(down_fraction + 0.5);

			double down_error, upper_error;

			down_error = global_sample_average - sample_averages[down_index - 2];
			upper_error = sample_averages[upper_index] - global_sample_average;

			
			v_errors[line][column - Nindependent] = max(down_error,upper_error);
		}
	}
}


