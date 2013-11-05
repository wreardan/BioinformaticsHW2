/*
Problem 3. UPGMA (15)
Write a program makeTree, which takes as input a set of sequences and outputs a tree inferred using the
UPGMA algorithm. You can use the code already written for Problem 1 to generate pairwise alignments
for these sequences. Estimate distance as the fractional mismatch for a pair of sequences. The makeTree
program should have the following usage:
makeTree sequence set.txt
Here sequence set.txt has the set of sequences in fasta format as shown below.
>s1
ATAATAA
>s2
ATCGATT
The format of your output should have one line per edge in the tree, each line has two columns, the
first column corresponding to the parent internal node, and the second column corresponding to the child
node. The leaf nodes are named using the name in the fasta file (after “>”). Name the internal nodes as a
concatenation of the names of the children it is merging, concatenated in lexicographical order. For example,
if our input file had sequences named, s1, s2 and s3, the output containing a possible tree would be
s1-s2 s1
s1-s2 s2
s1-s2-s3 s1-s2
s1-s2-s3 s3
*/

#include "Sequence.h"
#include "Matrix.h"

#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
	if(argc < 2) {
		cout << "usage: makeTree sequence_set.txt" << endl;
		return 1;
	}
	//Read in sequences
	vector<Sequence> sequence_set;
	Sequence::ReadFASTASet(sequence_set, argv[1]);
	//Perform ClustalW
	Matrix m;
	cout << m.UPGMA(sequence_set);
	//Output the result
#ifdef _DEBUG
	cout << "press enter to close" << endl;
	int nothing = getchar();
#endif
	return 0;
}
