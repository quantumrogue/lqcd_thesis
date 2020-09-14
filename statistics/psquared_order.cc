#include <iostream> //allows to write files
#include <fstream> //allows to read files
#include <vector> //Using vectors instead of normal arrays - range is not know a priori
#include <map>
#include <string> //allows to use strings - useful for filenames
#include <sstream>
#include <cmath> //use mathematical operators
#include <algorithm>
#include <iomanip>
using namespace std;

//DECLARE THIS VARIABLE HERE SO IT IS A GLOBAL VARIABLE - in order to be recognized by the lambda function below
int Pcolumn;

int* getfilesize(string filename);
void getfilescontent(vector<vector<double> >  &v_filescontent, const int Ncolumns, string filename);

int main(int argc, char** argv){

	int Nindependent = 6;

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
    //ORDER INFO BY P^2 - NUMBER OF COLUMN WHERE P^2 SHOULD BE THE SECOND ARGUMENT OF THE COMMAND LINE
	if (argc >= 3){
        std::istringstream iss( argv[2] );
        if (iss >> Pcolumn){
        }
    }
    else{
    	cout << "Second command line argument should be the column where the normal/sine momentum is." << endl;
    	return 0;
    }

	int line_number; //number of lines 
	int Ncolumns; //number of columns
	int *dimensionspointer;
	dimensionspointer = getfilesize(filename);
	line_number = dimensionspointer[0];
	Ncolumns = dimensionspointer[1];
	//printing stuff to command line
	cout << "Lines:	" << line_number << endl;
	cout << "Columns:	" << Ncolumns << endl;
	//define vector that will have all file information
	vector<vector<double> > v_filescontent(line_number, vector<double>(Ncolumns));

	getfilescontent(v_filescontent, Ncolumns, filename);

	//sorting
	sort(v_filescontent.begin(),v_filescontent.end(),
		[](const std::vector<double>& a, const std::vector<double>& b) {
		  return a[Pcolumn] < b[Pcolumn];
		});

	//create new file
	string newfile = filename;
	if(!newfile.empty()){
		newfile.resize(newfile.size() - 4);
	}
	newfile += "_sorted.dat";
	//write to new file
	ofstream myfile;
	myfile.open (newfile);

	for(int l = 0; l < line_number; l++){
		for(int c = 0; c < Ncolumns; c++){
			//only printing full precision for the dependent variables
			if(c < Nindependent){
				myfile << v_filescontent[l][c] << "	";
			}
			else{
				myfile << setprecision(15) << v_filescontent[l][c] << "	";
			}
		}
		myfile << endl;
	}
	myfile.close();

	return 0;
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