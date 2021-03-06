#pragma once
#include <string>
#include <vector>

class Sequence
{
protected:
public:
	std::string sequence_data;
	float score;
	std::vector<std::string> names;
	std::vector<Sequence> profile;
	bool used;

	Sequence(void);
	Sequence(char * filename);
	Sequence(std::string data);
	~Sequence(void);

	std::vector<Sequence> static ReadMultipleSequences(int argc, char *argv[]);
	static void ReadProfile(std::vector<Sequence> & result, char * filename);
	static void ReadFASTASet(std::vector<Sequence> & result, char * filename);
};

