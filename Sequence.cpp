#include "Sequence.h"
#include <iostream>
#include <fstream>

using namespace std;

//Create a Blank Sequence
Sequence::Sequence(void)
{
	sequence_data = "";
	score = 0;
	used = false;
}

vector<Sequence> static ReadMultipleSequences(int argc, char *argv[])
{
	return vector<Sequence>();
}

//Read in two aligned profiles
void Sequence::ReadProfile(std::vector<Sequence> & result, char * filename)
{
	string line;

	ifstream file(filename);
	if(!file.is_open()) {
		cerr << "file \"" << filename << "\" could not be opened" << endl;
		result.empty();
		return;
	}

	while(getline(file, line)) {
		result.push_back(Sequence(line));
	}

	return;
}

void Sequence::ReadFASTASet(std::vector<Sequence> & result, char * filename)
{
	string line;

	ifstream file(filename);
	if(!file.is_open()) {
		cerr << "file \"" << filename << "\" could not be opened" << endl;
		result.empty();
		return;
	}

	while(getline(file, line)) {
		Sequence temp;
		temp.names.push_back(line.substr(1));
		getline(file, line);
		temp.sequence_data = line;
		result.push_back(temp);
	}

	return;
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

Sequence::Sequence(std::string data)
{
	this->sequence_data = data;
}


Sequence::~Sequence(void)
{
}
