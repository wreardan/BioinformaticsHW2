#include "Matrix.h"
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <sstream>

using namespace std;

static float gap_penalty = 0;
static float matches = 2;
static float mismatches = 1;

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

			float cur_score = (seq1.sequence_data[x-1] == seq2.sequence_data[y-1]) ? matches : mismatches;
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

vector<Sequence> Matrix::MultipleAlignment(std::vector<Sequence> profile1, std::vector<Sequence> profile2)
{
	int x, y;
	vector<Sequence> results;
	vector<string> result1, result2, buf1, buf2;
	vector<int> scores;
	unsigned int i, j;

	SetupMatrix(profile1[0].sequence_data.size(), profile2[0].sequence_data.size());

	//Resize result/temp vectors
	result1.resize(profile1.size());
	result2.resize(profile2.size());
	buf1.resize(profile1.size());
	buf2.resize(profile2.size());

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

			float cur_score = 0;
			for(i = 0; i < profile1.size(); i++) {
				for(j = 0; j < profile2.size(); j++) {
					char a = profile1[i].sequence_data[x-1];
					char b = profile2[j].sequence_data[y-1];
					if(a == '-' || b == '-')
						continue;
					cur_score += (a == b) ? matches : mismatches;
				}
			}
			cur_score = cur_score / (profile1.size() * profile2.size());

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

		for(i = 0; i < profile1.size(); i++)
			if(x > 0)
				buf1[i] =  profile1[i].sequence_data[x-1];
			else
				buf1[i] = "-";

		for(i = 0; i < profile2.size(); i++)
			if(y > 0)
				buf2[i] =  profile2[i].sequence_data[y-1];
			else
				buf2[i] = "-";

		if(node.upleft >= node.left && node.upleft >= node.up) {
			for(i = 0; i < profile1.size(); i++)
				result1[i].append(buf1[i]);
			for(i = 0; i < profile2.size(); i++)
				result2[i].append(buf2[i]);
			x--; y--;
		}
		else if(node.left > node.up && node.left > node.upleft) {
			for(i = 0; i < profile1.size(); i++)
				result1[i].append(buf1[i]);
			for(i = 0; i < profile2.size(); i++)
				result2[i].append("-");
			x--;
		}
		else //if(node.up > node.left && node.up > node.upleft)
		{
			for(i = 0; i < profile1.size(); i++)
				result1[i].append("-");
			for(i = 0; i < profile2.size(); i++)
				result2[i].append(buf2[i]);
			y--;
		}
	}

	//Reverse Strings
	for(i = 0; i < profile1.size(); i++)
		reverse(result1[i].begin(), result1[i].end());
	for(i = 0; i < profile2.size(); i++)
		reverse(result2[i].begin(), result2[i].end());

	Sequence s;
	//Append Strings
	for(i = 0; i < profile1.size(); i++) {
		s.sequence_data = result1[i];
		results.push_back(s);
	}
	for(i = 0; i < profile2.size(); i++) {
		s.sequence_data = result2[i];
		results.push_back(s);
	}

	return results;
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


std::vector<Sequence> Matrix::UPGMA_Sequence(std::vector<Sequence> & sequence_set)
{
	std::vector<Sequence> alignment;
	Matrix m;
	Sequence seq;

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
		
		//combine_names(sequence_set[sequence_set.size()-1].names, sequence_set[min_x].names);
		//combine_names(sequence_set[sequence_set.size()-1].names, sequence_set[min_y].names);
		if(sequence_set[min_x].profile.size() == 0) {
			Sequence s;
			s.sequence_data = sequence_set[min_x].sequence_data;
			sequence_set[min_x].profile.push_back(s);
		}
		if(sequence_set[min_y].profile.size() == 0) {
			Sequence s;
			s.sequence_data = sequence_set[min_y].sequence_data;
			sequence_set[min_y].profile.push_back(s);
		}

		seq = Sequence();
		//seq.profile.insert(seq.profile.begin(), sequence_set[min_x].profile.begin(), sequence_set[min_x].profile.end());
		//seq.profile.insert(seq.profile.begin(), sequence_set[min_y].profile.begin(), sequence_set[min_y].profile.end());
		seq.profile = m.MultipleAlignment(sequence_set[min_y].profile, sequence_set[min_x].profile);
		
		sequence_set.push_back(seq);
	}



	return seq.profile;
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
