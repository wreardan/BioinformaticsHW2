#pragma once
#include <vector>
#include <string>
#include "Sequence.h"
class MatrixNode
{
public:
	int left, upleft, up, score, letter;

	MatrixNode();
	MatrixNode(unsigned int x, unsigned int y);
};
class Matrix
{
protected:

	unsigned int width, height;

	void SetupMatrix();
public:
	std::vector<std::vector<MatrixNode>> data;
	Sequence * seq1;
	Sequence * seq2;

	Matrix(Sequence & s1, Sequence & s2);
	~Matrix(void);

	std::string globalAlign();
};

