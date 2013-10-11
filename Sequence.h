#pragma once
#include <string>
class Sequence
{
protected:
public:
	std::string sequence_data;

	Sequence(void);
	Sequence(char * filename);
	~Sequence(void);
};

