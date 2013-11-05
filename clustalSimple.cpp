#include "Sequence.h"
#include "Matrix.h"

#include <iostream>
#include <vector>

using namespace std;

int main4(int argc, char *argv[])
{
	if(argc < 2) {
		cout << "usage: clustalSimple sequence_set.txt" << endl;
		return 1;
	}
	//Read in sequences
	vector<Sequence> sequence_set;
	Sequence::ReadFASTASet(sequence_set, argv[1]);
	//Perform ClustalW
	Matrix m;
	vector<Sequence> sequences = m.UPGMA_Sequence(sequence_set);
	for(unsigned int i = 0; i < sequences.size(); i++)
		cout << sequences[i].sequence_data << endl;
	//Output the result
#ifdef _DEBUG
	cout << "press enter to close" << endl;
	int nothing = getchar();
#endif
	return 0;
}
