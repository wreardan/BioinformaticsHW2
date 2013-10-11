/*
Problem 2. Profile alignment (10)
Write a program alignProfile which aligns two profiles using average of scores from pairwise comparisons
of sequences, one from each alignment. This is similar to the unweighted scoring procedure in
ClustalW. As in Problem 1, a match has a score of 2, mismatch a score of 1 and a gap penalty 0. The
program should have the following usage (shown only for C++),
1
alignProfile seq1.txt seq2.txt
Here seq1.txt and seq2.txt correspond to files containing sequence profiles with separate lines for
each sequence. The program should output an alignment of the sequence profiles.
For example if seq1.txt’s contents are
AAAC
AGAC
and seq2.txt’s contents are
AGC
ACC
The output should be
AAAC
AGAC
AG-C
AC-C
*/

#include "Sequence.h"
#include "Matrix.h"
#include <iostream>

using namespace std;
/*
int main(int argc, char *argv[])
{
	return 0;
}
*/