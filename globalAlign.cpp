/* 
BMI/CS 576 Introduction to Bioinformatics
Fall 2013 - Wesley Reardan

Problem 1. Global alignment (10)
Write a program globalAlign that implements the Needleman-Wunsch global alignment algorithm for
a pair of sequences. Your program should take two sequences as two input files and print out the alignment
on the command line. Assume a score matrix for DNA sequence with matches scored at 2, mismatches at 1
and gap penalty of 0.

You can use either C++, Java, Perl or Python to implement this program. We should be able to run it
using one of the following commands depending upon the programming language you selected.

C++: ./globalAlign seq1.txt seq2.txt
Java: java globalAlign seq1.txt seq2.txt
Perl: perl globalAlign seq1.txt seq2.txt
Python: python globalAlign seq1.txt seq2.txt

Here seq1.txt and seq2.txt each have one sequence. Assume that the entire sequence is on one
line, that is, there are no line breaks. The program should output an alignment of the two sequences. For
example, if seq1.txt has AAGG and seq1.txt has AA, the alignment printed out should be

AAGG
AA--

*/

#include "Sequence.h"
#include "Matrix.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
	if(argc < 3) {
		cout << "usage: ./globalAlign seq1.txt seq2.txt" << endl;
		return 1;
	}
	Sequence seq1(argv[1]);
	Sequence seq2(argv[2]);

#ifdef _DEBUG
	cout << "Original Alignments: " << endl;
	cout << "\t" << seq1.sequence_data << endl;
	cout << "\t" << seq2.sequence_data << endl;
#endif

	Matrix matrix(seq1, seq2);
	Matrix matrix2(seq2, seq1);
	cout << matrix.globalAlign() << endl;
	cout << matrix2.globalAlign() << endl;

#ifdef _DEBUG
	cout << "press enter to close" << endl;
	int nothing = getchar();
#endif
	return 0;
}