#include "Matrix.h"
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <sstream>

using namespace std;

static int gap_penalty = 0;
static int matches = 2;
static int mismatches = 1;

//Create a Blank Matrix Node
MatrixNode::MatrixNode()
{
	this->left = 0;
	this->up = 0;
	this->upleft = 0;
	this->score = 0;
	child1 = child2 = -1;
}

//Create a Matrix Node in the grid
MatrixNode::MatrixNode(unsigned int x, unsigned int y)
{
	this->upleft = -1000000;
	//Left Column, Up is the only choice
	if(x == 0) {
		this->left = -1000000;
		this->up = this->score = -gap_penalty * y;
	//Top Row, Left is the only choice
	} else if(y == 0) {
		this->left = this->score = -gap_penalty * x;
		this->up = -1000000;
	} else {
		assert(false);
	}
}

void Matrix::SetupMatrix(unsigned int width, unsigned int height)
{
	this->width = ++width;
	this->height = ++height;
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
Sequence Matrix::GlobalAlign(Sequence & seq1, Sequence & seq2)
{
	int x, y;
	Sequence result;

	SetupMatrix(seq1.sequence_data.size(), seq2.sequence_data.size());

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

			unsigned int cur_score = (seq1.sequence_data[x-1] == seq2.sequence_data[y-1]) ? matches : mismatches;
			node.upleft = data[y-1][x-1].score + cur_score;

			if(node.upleft >= node.left && node.upleft >= node.up)
				node.score = node.upleft;
			else if(node.left > node.up && node.left > node.upleft)
				node.score = node.left;
			else //if(node.up > node.left && node.up > node.upleft)
				node.score = node.up;
		}
	}

	//build alignment result
	x = width - 1;
	y = height - 1;
	while(x > 0 || y > 0) {
		MatrixNode & node = data[y][x];
		char buf[4] = {'-', 0, 0, 0 };
		if(x > 0)
			buf[0] = seq1.sequence_data[x-1];
		if(node.upleft >= node.left && node.upleft >= node.up) {
			result.sequence_data.append(buf);
			x--; y--;
		}
		else if(node.up > node.left && node.up > node.upleft)
		{
			result.sequence_data.append("-");
			y--;
		}
		else //if(node.left > node.up && node.left > node.upleft)
		{	
			//result.sequence_data.append("-");
			result.sequence_data.append(buf);
			x--;
		}
	}
	//reverse resulting string because it's backwards
	reverse(result.sequence_data.begin(), result.sequence_data.end());

#ifdef _DEBUG
	stringstream ss;
	ss << result.sequence_data << " score: " << data[height-1][width-1].score;

	//cout << ss.str() << endl;
#endif
	result.score = data[height-1][width-1].score;
	return result;
}

string Matrix::MultipleAlignment(std::vector<Sequence> profile1, std::vector<Sequence> profile2)
{
	int x, y;
	stringstream resultString;
	string result1, result2, result3, result4;

	SetupMatrix(profile1[0].sequence_data.size(), profile2[0].sequence_data.size());

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

			int score1=0, score2=0, score3=0, score4=0;
			char a = profile1[0].sequence_data[x-1];
			char b = profile1[1].sequence_data[x-1];
			char c = profile2[0].sequence_data[y-1];
			char d = profile2[1].sequence_data[y-1];
			if(a != '-') {
				if(c != '-')
					score1 = (a == c) ? matches : mismatches;
				if(d != '-')
					score2 = (a == d) ? matches : mismatches;
			}
			if(b != '-') {
				if(c != '-')
					score3 = (b == c) ? matches : mismatches;
				if(d != '-')
					score4 = (b == d) ? matches : mismatches;
			}
			int cur_score = (score1 + score2 + score3 + score4);
			node.upleft = data[y-1][x-1].score + cur_score;

			if(node.upleft >= node.left && node.upleft >= node.up)
				node.score = node.upleft;
			else if(node.left > node.up && node.left > node.upleft)
				node.score = node.left;
			else //if(node.up > node.left && node.up > node.upleft)
				node.score = node.up;
		}
	}

	//build alignment result
	x = width - 1;
	y = height - 1;
	while(x > 0 || y > 0) {
		MatrixNode & node = data[y][x];
		char buf1[4] = {'-', 0, 0, 0 };
		char buf2[4] = {'-', 0, 0, 0 };
		char buf3[4] = {'-', 0, 0, 0 };
		char buf4[4] = {'-', 0, 0, 0 };
		if(x > 0) {
			buf1[0] = profile1[0].sequence_data[x-1];
			buf2[0] = profile1[1].sequence_data[x-1];
		}
		if(y > 0) {
			buf3[0] = profile2[0].sequence_data[y-1];
			buf4[0] = profile2[1].sequence_data[y-1];
		}
		if(node.upleft > node.left && node.upleft > node.up) {
			result1.append(buf1);
			result2.append(buf2);
			result3.append(buf3);
			result4.append(buf4);
			x--; y--;
		}
		else if(node.left > node.up && node.left > node.upleft)
		{	
			result1.append(buf1);
			result2.append(buf2);
			result3.append("-");
			result4.append("-");
			x--;
		}
		else //if(node.up > node.left && node.up > node.upleft)
		{
			result1.append("-");
			result2.append("-");
			result3.append(buf3);
			result4.append(buf4);
			y--;
		}
	}
	//reverse resulting strings because they're backwards
	reverse(result1.begin(), result1.end());
	reverse(result2.begin(), result2.end());
	reverse(result3.begin(), result3.end());
	reverse(result4.begin(), result4.end());

	resultString << result1 << endl << result2 << endl << result3 << endl << result4 << endl;
	return resultString.str();
}

/* HW2-Problem2
ClustalW Overview
-Generate pairwise alignments for the two profiles
*/
string Matrix::ClustalW(std::vector<Sequence> profile1, std::vector<Sequence> profile2)
{
	stringstream ss;
	Matrix m;
	Sequence aligned1, aligned2;
	SetupMatrix(profile1.size()-1, profile2.size()-1);
	int x, y;
	for(y = 0; y < height; y++) {
		for(x = 0; x < width; x++) {
			aligned1 = m.GlobalAlign(profile2[x], profile1[y]);
			aligned2 = m.GlobalAlign(profile1[y], profile2[x]);
			int similarity = 0;
			for(unsigned int i = 0; i < aligned1.sequence_data.size(); i++)
				if(aligned1.sequence_data[i] == aligned2.sequence_data[i])
					similarity++;
			float distance = 1.0f - ((float)similarity / aligned1.sequence_data.size());
			data[y][x].distance = distance;
			ss << distance << "\t";
		}
		ss << endl;
	}
	return ss.str();
}

void Matrix::DistanceMatrix(std::vector<Sequence> & sequence_set)
{
	Matrix m;
	Sequence aligned1, aligned2;
	SetupMatrix(sequence_set.size()-1, sequence_set.size()-1);
	int x, y;
	for(y = 0; y < height; y++) {
		for(x = y+1; x < width; x++) {
			aligned1 = m.GlobalAlign(sequence_set[x], sequence_set[y]);
			aligned2 = m.GlobalAlign(sequence_set[y], sequence_set[x]);
			int similarity = 0;
			for(unsigned int i = 0; i < aligned1.sequence_data.size(); i++)
				if(aligned1.sequence_data[i] == aligned2.sequence_data[i])
					similarity++;
			float distance = 1.0f - ((float)similarity / aligned1.sequence_data.size());
			data[y][x].distance = distance;
			data[x][y].distance = distance;
		}
	}
}

void print_names(vector<string> & names, stringstream & ss) {
	sort(names.begin(), names.end());
	unsigned int i;
	for(i = 0; i < names.size()-1; i++)
		ss << names[i] << "-";
	ss << names[i];
}

void combine_names(vector<string> & names1, vector<string> & names2) {
	names1.insert(names1.end(), names2.begin(), names2.end());
}



string Matrix::UPGMA(std::vector<Sequence> & sequence_set)
{
	stringstream ss;
	vector<string> names;

	DistanceMatrix(sequence_set);

	for(int i = 0; i < max(width, height); i++) {
		//Find minimum distance
		int min_y = 0, min_x = 1;
		float min_distance = 1000000;
		bool all_used = true;
		for(int y = 0; y < height; y++) {
			for(int x = y+1; x < width; x++) {
				if(sequence_set[y].used || sequence_set[x].used)
					continue;
				if(data[y][x].distance < min_distance) {
					min_distance = data[y][x].distance;
					min_y = y;
					min_x = x;
					all_used = false;
				}
			}
		}
		if(all_used)
			break;

		SetupMatrix(height, width);
		for(int y = 0; y < height-1; y++) {
			int x = width - 1;
			MatrixNode & node = data[y][x];
			if(x == min_x || y == min_y)
				node.distance = 10000;
			else
				node.distance = (data[y][min_x].distance + data[y][min_y].distance) / 2.0f;
			node.child1 = min_x;
			node.child2 = min_y;
		}

		sequence_set[min_y].used = true;
		sequence_set[min_x].used = true;
		/*if(min_x < sequence_set.size())
			names.push_back(sequence_set[min_x].name);
		if(min_y < sequence_set.size())
			names.push_back(sequence_set[min_y].name);*/
		names.clear();
		combine_names(names, sequence_set[min_x].names);
		combine_names(names, sequence_set[min_y].names);
		print_names(names, ss);
		ss << " ";
		print_names(sequence_set[min_y].names, ss);
		ss << endl;
		print_names(names, ss);
		ss << " ";
		print_names(sequence_set[min_x].names, ss);
		ss << endl;

		sequence_set.resize(sequence_set.size() + 1);
		combine_names(sequence_set[sequence_set.size()-1].names, sequence_set[min_x].names);
		combine_names(sequence_set[sequence_set.size()-1].names, sequence_set[min_y].names);
	}


	return ss.str();
}

Matrix::Matrix()
{
	this->width = 0;
	this->height = 0;
}

Matrix::~Matrix(void)
{
}
