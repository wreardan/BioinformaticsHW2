#include "Sequence.h"
#include <iostream>
#include <fstream>

using namespace std;

//Create a Blank Sequence
Sequence::Sequence(void)
{
	sequence_data = "";
}

//Create a sequence from a file
Sequence::Sequence(char * filename)
{
	ifstream file(filename);
	if(!file.is_open()) {
		cerr << "file \"" << filename << "\" could not be opened" << endl;
		sequence_data = "";
	} else {
		getline(file, sequence_data);
	}

}

Sequence::~Sequence(void)
{
}
