

#include "Tree.h"
#include "MappingParser.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

class BranchStat {

	SimpleStateSpace* sstatespace;
	int samplesize;
	double meansub;
	double branchlength;


public:

	BranchStat(SimpleStateSpace* instatespace) {
		sstatespace = instatespace;
		samplesize = 0;
		meansub = 0;
		branchlength = 0;
	};


	void addnewstat(string input) {
		samplesize++;
		SimpleMappingParser* m = new SimpleMappingParser(input, sstatespace);

		meansub += m->GetNbSub();
		branchlength += m->GetTotalLength();

		delete m;
	}

	void Normalize() {
		meansub /= samplesize;
		branchlength /= samplesize;
	}



	double GetMeanSub() {
		return meansub;
	}

	double GetBranchLength(){
		return branchlength;
	}





};


class ReadMapTreeAA : public Tree {
	SimpleStateSpace* sstatespace;
	BranchStat** branchstats;

	double** profiles;


public:

	ReadMapTreeAA(string filename) : Tree(filename) {
		sstatespace = new ProteinStateSpace();

		SetIndices();

		branchstats = new BranchStat*[GetNbranch()];
		for (int j = 1; j < GetNbranch(); j++) {
			branchstats[j] = new BranchStat(sstatespace);
		}

		profiles = new double*[GetNbranch()*2+1];
		for (int i = 1; i < GetNbranch()*2+1; i++) {
			profiles[i] = new double[sstatespace->GetNstate()];
			for (int j = 1; j < sstatespace->GetNstate(); j++) {
				profiles[i][j] = 0;
			}
		}

	};

	~ReadMapTreeAA() {
		for (int j = 1; j < GetNbranch(); j++) {
			delete branchstats[j];
		}
		delete[] branchstats;
	};

	void JsonToStream(ostream& os) {
		os << "{\n";
		JsonToStream(GetRoot(), "", os);
		os << "}\n";
	}

	string Rereadtree(Link * from, string input);


	void RecursiveNormalise(Link* from) {
		for (Link* link = from->Next(); link != from; link = link->Next()) {
				RecursiveNormalise(link->Out());
		}
		if (!from->isRoot()) {
				GetBranchStat(from)->Normalize();

		}

	}

	double GetMeanSub(){
		return RecursiveAddValue(GetRoot(), &BranchStat::GetMeanSub);
	}


private:



	BranchStat* GetBranchStat(Link* from) {
		if (from->isRoot()) {
				cerr << "The root should not have been call";
		}
		return branchstats[from->GetBranch()->GetIndex()];
	}


	double RecursiveAddValue(Link* from, double(BranchStat::*ptr)());


	void JsonToStream(Link* from, string tab, ostream& os);

};
