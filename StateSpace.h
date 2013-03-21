

#ifndef STATESPACE_H
#define STATESPACE_H

#include "BiologicalSequences.h"

// pure interface
//
class StateSpace	{

	public:

	virtual ~StateSpace() {}

	virtual string GetState(int state) = 0;
	virtual int GetNstate() = 0;

	virtual int GetState(string from) = 0;

};

// simple state space: assumes that states are referred to using a one-letter code
//
class SimpleStateSpace : public StateSpace	{

	public:


	int GetState(string from);

	int GetNstate() {
		return Nstate;
	}
	
	string GetState(int state);

	protected:
	int Nstate;
	char* Alphabet;
	int NAlphabetSet;
	char* AlphabetSet;
};

class DNAStateSpace : public SimpleStateSpace	{

	public:

	DNAStateSpace();
	~DNAStateSpace();
};

class RNAStateSpace : public SimpleStateSpace	{

	public:

	RNAStateSpace();
	~RNAStateSpace();
};

class ProteinStateSpace : public SimpleStateSpace	{

	public:

	ProteinStateSpace();
	~ProteinStateSpace();
};

#endif // STATESPACE_H

