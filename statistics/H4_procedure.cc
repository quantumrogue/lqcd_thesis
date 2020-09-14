#include "LinearFit.h" //Class defined for linear regressions
#include <iostream> //allows to write files
#include <fstream> //allows to read files
#include <vector> //Using vectors instead of normal arrays - range is not know a priori
#include <cmath> //use mathematical operators
#include <iomanip>
using namespace std;
using namespace fit;

int* getfilesize(string filename);
void getfilescontent(vector<vector<double> >  &v_filescontent, const int Ncolumns, string filename);

int main(int argc, char** argv){

	int Nd = 4;
	double Pi = 3.14159265358979323846264338;

	//READ FILE
	//Obtaining the name of the file to open
    string filename;
    if (argc >= 2){
    	filename = argv[1];
    }
    else{
    	cout << "First command line argument should be the filename" << endl;
    	return 0;
    }
    //column of the ordered momentum
    int Pcolumn;
	if (argc >= 3){
        std::istringstream iss( argv[2] );
        if (iss >> Pcolumn){
        }
    }
    else{
    	cout << "Second command line argument should be the column where the normal/sine momentum is." << endl;
    	return 0;
    }
    int Nindependent;
    if (argc >= 4){
        std::istringstream iss( argv[3] );
        if (iss >> Nindependent){
        }
    }
    else{
    	cout << "Third command line argument should be the number of independent variables." << endl;
    	return 0;
    }
    int latticesize;
    if (argc >= 5){
        std::istringstream iss( argv[4] );
        if (iss >> latticesize){
        }
    }
    else{
    	cout << "Fourth command line argument should be the lattice size." << endl;
    	return 0;
    }

    int line_number; //number of lines 
	int Ncolumns; //number of columns
	int *dimensionspointer;
	dimensionspointer = getfilesize(filename);
	line_number = dimensionspointer[0];
	Ncolumns = dimensionspointer[1];
	//define vector that will have all file information
	vector<vector<double> > v_filescontent(line_number, vector<double>(Ncolumns));
	getfilescontent(v_filescontent, Ncolumns, filename);

	//LINEAR REGRESSION PROCEDURE
	//iterate through same momentum points
	double pmin,pmax;
	pmin = 0;
	pmax = 3;
	double momentum = v_filescontent[0][Pcolumn];
	
	//plus 1 because will add p^4
	int finaldimension = Ncolumns - Nindependent + 1;
	//vectors to store the info and groups of points
	vector<vector<double> > v_pointstofit(finaldimension, vector<double>());
	int multiplicity = 0;

	//create new file
	string newfile = filename;
	if(!newfile.empty()){
		newfile.resize(newfile.size() - 4);
	}
	newfile += "_H4.dat";
	//write to new file
	ofstream myfile;
	myfile.open(newfile);

	//cycle through content vector
	for(int line = 0; line < line_number; line++){

		double p4;	
		//change momentum
		double newmomentum = v_filescontent[line][Pcolumn];
		if(newmomentum > momentum){

			int nr_of_values = (finaldimension-1)/2;

			if(multiplicity > 2 && (momentum > pmin && momentum < pmax)){

				myfile << setprecision(15) << momentum << "	";

				//only perform linear regression if multiplicity > 1 and momentum falls inside considered range [p_min, p_max]
				for(int j = 1; j < nr_of_values+1; j++){
				
					LinearFit points(multiplicity);

					//vector to store infor in the case of multiplicity = 2
					vector<vector<double> > v_mult2(2, vector<double>(3));

					for(int point = 0; point < multiplicity; point++){
						//compute p^4
						p4 = 0;
						for(int i = 0; i < Nd; i++){
							double p_i = v_filescontent[line-multiplicity+point][i];
							p4 += p_i*p_i*p_i*p_i;
						}
						p4 *= pow((2*Pi/latticesize),4);
						points.datapoints[point][0] = p4;
						points.datapoints[point][1] = v_filescontent[line - multiplicity + point][2*j-1 + Nindependent-1];
						points.datapoints[point][2] = v_filescontent[line - multiplicity + point][2*j + Nindependent-1];
						if(multiplicity == 2){
							v_mult2[point][0] = p4;
							v_mult2[point][1] = v_filescontent[line - multiplicity + point][2*j-1 + Nindependent-1];
							v_mult2[point][2] = v_filescontent[line - multiplicity + point][2*j + Nindependent-1];
 						}
					}

					//linear regression - y= a + bx
					double *parameters = points.fitparameters();
					double a = parameters[0];
					double delta_a;
					if(multiplicity == 2){
						//treatment for 2 point fit
						double dy1, dy2, y1, y2, x1, x2, delta_a1, delta_a2;
						y1= v_mult2[0][1];
						dy1 = v_mult2[0][2];
						y2 = v_mult2[1][1];
						dy2 = v_mult2[1][2];
						x1 = v_mult2[0][0];
						x2 = v_mult2[1][0];
						double dx = x2 - x1;
						delta_a1 = sqrt((dy1*(1+x1/dx))*(dy1*(1+x1/dx)) + (x1*dy2/dx)*(x1*dy2/dx));
						delta_a2 = sqrt((x2*dy1/dx)*(x2*dy1/dx) + (dy2*(1-x2/dx))*(dy2*(1-x2/dx)));
						delta_a = max(delta_a1,delta_a2);
					}
					else{
						delta_a = parameters[1];	
					}
					myfile << setprecision(15) << a << " " << delta_a << "	";
				}//for column j
				myfile << endl;
			}//if multiplicity 
			else{
				//in the case of the momentum considered being out of the range
				//will disregard completely cases when there are only 2 points to perform regression
				if(multiplicity > 2 || momentum <	 pmin){
				for(int point = 0; point < multiplicity; point++){
					myfile << setprecision(7) << momentum << "	";
					for(int j = 1; j < nr_of_values+1; j++){
						myfile << setprecision(15) << v_filescontent[line - multiplicity + point][2*j-1 + Nindependent-1] << "   "; 
						myfile << setprecision(15) << v_filescontent[line - multiplicity + point][2*j + Nindependent-1] << "   ";
					}
					myfile << endl;
				}
				}
			}
			multiplicity = 0;
		}//reset momentum
		++multiplicity;
		momentum = newmomentum;
	}//line cycle
	myfile.close();
}




//This function counts the number of lines of the first document/configuration (vertex1.dat)
int* getfilesize(string filename){
	
	//string firstfilename = "vertex_data/vertex1.dat";
	string firstfilename = filename;
	//Open the first file to count number of lines
	ifstream firstfile;
	int line_number = 0;
	int Ncolumns = 0;
	firstfile.open(firstfilename);
	string line;
	int i;

	if (firstfile.is_open()){
		while (getline(firstfile, line)){
			if(line_number == 0){
				stringstream ss(line);
				double value;
				while(ss >> value){
		            ++Ncolumns; 
	        	}
			}
			line_number++;
		}
	}
	else {
		cout << "Could not open file! - Ending program." << endl;
		return 0;
	}
	firstfile.close();


	int dimensions[2] = {line_number, Ncolumns};
	int *dimensionspointer;
	dimensionspointer = dimensions;
	return dimensionspointer;
}

void getfilescontent(vector<vector<double> >  &v_filescontent, const int Ncolumns, string filename){

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
					v_filescontent[line_count][column] = value;
			}
			line_count++;
		}
		file.close();
	}
}