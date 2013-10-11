#include "Matrix.h"
#include <assert.h>
#include <algorithm>

using namespace std;

const int gap_penalty = 0;
const int matches = 2;
const int mismatches = 1;

//Create a Blank Matrix Node
MatrixNode::MatrixNode()
{
	this->left = 0;
	this->letter = 0;
	this->up = 0;
	this->upleft = 0;
	this->score = 0;
}

//Create a Matrix Node in the grid
MatrixNode::MatrixNode(unsigned int x, unsigned int y)
{
	this->upleft = -1000000;
	//Left Column, Up is the only choice
	if(x == 0) {
		this->left = -1000000;
		this->up = this->score = gap_penalty * y;
		this->letter = '-';
	//Top Row, Left is the only choice
	} else if(y == 0) {
		this->left = this->score = gap_penalty * x;
		this->up = -1000000;
		this->letter = '-';
	} else {
		assert(false);
	}
}

void Matrix::SetupMatrix()
{
	width = seq1->sequence_data.length() + 1;
	height = seq2->sequence_data.length() + 1;

	//Create 2d array - the easy way (using vector)
	data.resize(height);
	for(unsigned int i = 0; i < height; i++) {
		data[i].resize(width);
	}
}

/* Perform Global Alignment on two sequences in Matrix
TODO: modify to accept a 4x4 scoring matrix instead of using fixed values

Calculate matrix starting at top left (x=0,y=0)

Recurrence Relation:
F(j,i) = max | F(j-1,i-1) + S(xi,yj)
             | F(j-1, i) - d
             | F(j, i-1) - d


      +--- First Column, all gaps
      |
+---+---+---+---+---+---+---+
|   |   | A | G | A | T | T |
+---+---+---+---+---+---+---+
|   | 0 |-3 |-6 |-9 |-12|-15|---First Row, all gaps
+---+---+---+---+---+---+---+
| A |-3 | 1 |-2 |-5 |-8 |-11|
+---+---+---+---+---+---+---+
| G |-6 |-2 | 2 |-1 |-4 |-7 |
+---+---+---+---+---+---+---+
| T |-9 |-5 |-1 | 1 | 0 |-3 |
+---+---+---+---+---+---+---+
| T |-12|-8 |-4 |-2 | 2 | 1 |
+---+---+---+---+---+---+---+---End Node

Reverse from End Node to find sequence

*/
string Matrix::globalAlign()
{
	unsigned int x, y;
	string result = "";

	//Top left box is 0,0,0,0
	data[0][0] = MatrixNode();

	//fill in first column
	for(y = 1; y < height; y++) {
		data[y][0] = MatrixNode(0,y);
	}

	//fill in first row
	for(x = 1; x < width; x++) {
		data[0][x] = MatrixNode(x,0);
	}

	//fill in the rest of the matrix
	for(y = 1; y < height; y++) {
		for(x = 1; x < width; x++) {
			MatrixNode & node = data[y][x];
			node = MatrixNode();

			node.left = data[y][x-1].score - gap_penalty;
			node.up = data[y-1][x].score - gap_penalty;

			unsigned int cur_score = (seq1->sequence_data[x-1] == seq2->sequence_data[y-1]) ? matches : mismatches;
			node.upleft = data[y-1][x-1].score + cur_score;

			if(node.left > node.up && node.left > node.upleft)
			{
				node.letter = '-';
				node.score = node.left;
			}
			else if(node.up > node.left && node.up > node.upleft)
			{
				node.letter = '-';
				node.score = node.up;
			}
			else
			{
				node.letter = seq1->sequence_data[x];
				node.score = node.upleft;
			}
		}
	}

	//build alignment result
	x = width - 1;
	y = height - 1;
	while(x > 0 && y > 0) {
		MatrixNode & node = data[y][x];
		char buf[4] = { seq1->sequence_data[x-1], 0, 0, 0 };
		if(node.left > node.up && node.left > node.upleft)
		{	
			result.append(buf);
			x--;
		}
		else if(node.up > node.left && node.up > node.upleft)
		{
			result.append("-");
			y--;
		}
		else
		{
			result.append(buf);
			x--; y--;
		}
	}
	//reverse resulting string because it's backwards
	reverse(result.begin(), result.end());

	return result;
}

Matrix::Matrix(Sequence & s1, Sequence & s2)
{
	if(false && s1.sequence_data.length() < s2.sequence_data.length()) {
		seq1 = &s2; //turned off at the moment
		seq2 = &s1;
	}
	else {
		seq1 = &s1;
		seq2 = &s2;
	}
	SetupMatrix();
}


Matrix::~Matrix(void)
{
}
