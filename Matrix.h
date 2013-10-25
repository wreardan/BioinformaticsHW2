#pragma once
#include <vector>
#include <string>
#include "Sequence.h"

class MatrixNode
{
public:
	float left, upleft, up, score, distance;
	int child1, child2;

	MatrixNode();
	MatrixNode(unsigned int x, unsigned int y);
};

class Matrix
{
protected:
	void SetupMatrix(unsigned int width, unsigned int height);
	int substitutionMatrix[4][4];

public:
//Variables
	std::vector<std::vector<MatrixNode> > data;
	int width, height;

//Constructors
	Matrix();
	~Matrix(void);

//Methods
	Sequence GlobalAlign(Sequence & seq1, Sequence & seq2);
	std::string ClustalW(std::vector<Sequence> profile1, std::vector<Sequence> profile2);
	std::vector<Sequence> MultipleAlignment(std::vector<Sequence> profile1, std::vector<Sequence> profile2);
	std::string UPGMA(std::vector<Sequence> & sequence_set);
	std::vector<Sequence> UPGMA_Sequence(std::vector<Sequence> & sequence_set);
	void DistanceMatrix(std::vector<Sequence> & sequence_set);
};

