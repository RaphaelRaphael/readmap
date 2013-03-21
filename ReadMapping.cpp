
#include "Tree.h"
#include "CodonStateSpace.h"
#include "MappingParser.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;




class BranchStat {

	CodonStateSpace* cstatespace;
	int samplesize;
	double meansub;
	double meansyn;
	double meannonsyn;
	double sub1;
	double sub2;
	double sub3;
	double branchlength;


public:

	BranchStat(CodonStateSpace* instatespace) {
		cstatespace = instatespace;
		samplesize = 0;
		meansub = 0;
		meansyn = 0;
		meannonsyn = 0;
		sub1 = 0;
		sub2 = 0;
		sub3 = 0;
		branchlength = 0;
	};


	void addnewstat(string input) {
		samplesize++;
		CodonMappingParser* m = new CodonMappingParser(input, cstatespace);

		meannonsyn += m->GetNbNonSynSub();
		meansub += m->GetNbSub();
		meansyn += m->GetNbSynSub();
		sub1 += m->GetNbSingleSub();
		sub2 += m->GetNbDoubleSub();
		sub3 += m->GetNbTrippleSub();
		branchlength += m->GetTotalLength();

		delete m;
	}

	void Normalize() {
		meansub /= samplesize;
		meannonsyn /= samplesize;
		meansyn /= samplesize;
		sub1 /= samplesize;
		sub2 /= samplesize;
		sub3 /= samplesize;
		branchlength /= samplesize;
	}



	double GetMeanNonSyn() {
		return meannonsyn;
	}

	double GetMeanSyn() {
		return meansyn;
	}

	double GetMeanSub() {
		return meansub;
	}

	double GetMeanSingleSub() {
		return sub1;
	}

	double GetMeanDoubleSub() {
		return sub2;
	}

	double GetMeanTrippleSub() {
		return sub3;
	}

	double GetBranchLength(){
		return branchlength;
	}





};

class ReadMapTree : public Tree {
	CodonStateSpace* statespace;
	BranchStat** branchstats;


public:

	ReadMapTree(string filename) : Tree(filename) {
		statespace = new CodonStateSpace(Universal);

		SetIndices();

		branchstats = new BranchStat*[GetNbranch()];

		for (int j = 1; j < GetNbranch(); j++) {
				branchstats[j] = new BranchStat(statespace);

		}
	};

	~ReadMapTree() {
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

	string Rereadtree(Link * from, string input) {

		// Position of the last character of the mapping
		int end = 0;

		if (from->isLeaf()) {
			int last = 0;
			while ((input[end] != ',') and (input[end] != ')')) {
				if (input[end] == '_') {last = end; }
				end++;
			}

			// to remove the name of the leaf
			input = input.substr(last, -1);
			end -= last;
		} else {
			// to remove an open bracket
			input = input.substr(1, -1);

			//recursion
			for (Link* link = from->Next(); link != from; link = link->Next()) {
				input = Rereadtree(link->Out(), input);
			}

			// to find the last character of the mapping
			while ((input[end] != ',') and (input[end] != ')') and (input[end] != ';')) {
					end++;
			}
		}

		// send the maping 
		if (!from->isRoot()) {
			GetBranchStat(from)->addnewstat(input.substr(1, end - 1));
		}

		// to remove the maping, and a closing bracket or a coma
		input = input.substr(end + 1, -1);

		return input;
	}



	void RecursiveNormalise(Link* from) {
		for (Link* link = from->Next(); link != from; link = link->Next()) {
				RecursiveNormalise(link->Out());
		}
		if (!from->isRoot()) {
				GetBranchStat(from)->Normalize();

		}

	}

	double GetMeanNonSyn(){
		return RecursiveAddValue(GetRoot(), &BranchStat::GetMeanNonSyn);
	}

	double GetMeanSyn(){
		return RecursiveAddValue(GetRoot(), &BranchStat::GetMeanSyn);
	}

	double GetMeanSub(){
		return RecursiveAddValue(GetRoot(), &BranchStat::GetMeanSub);
	}

	double GetMeanSingleSub(){
		return RecursiveAddValue(GetRoot(), &BranchStat::GetMeanSingleSub);
	}

	double GetMeanDoubleSub(){
		return RecursiveAddValue(GetRoot(), &BranchStat::GetMeanDoubleSub);
	}

	double GetMeanTrippleSub(){
		return RecursiveAddValue(GetRoot(), &BranchStat::GetMeanTrippleSub);
	}


private:



	BranchStat* GetBranchStat(Link* from) {
		if (from->isRoot()) {
				cerr << "The root should not have been call";
		}
		return branchstats[from->GetBranch()->GetIndex()];
	}


	double RecursiveAddValue(Link* from, double(BranchStat::*ptr)()  ){
		double value = 0;
		for (Link* link = from->Next(); link != from; link = link->Next()) {
				value += RecursiveAddValue(link->Out(), ptr);
		}
		if (!from->isRoot()) {
				value += (GetBranchStat(from)->*ptr)();
		}
		return value;
	}


	void JsonToStream(Link* from, string tab, ostream& os) {
		tab += ' ';
		tab += ' ';
		if (!from->isRoot()) {
				os << tab << "\"length\": \"" << GetBranchStat(from)->GetBranchLength() << "\",\n";
				os << tab << "\"Sub\": \"" << GetBranchStat(from)->GetMeanSub() << "\",\n";
				os << tab << "\"SingleSub\": \"" << GetBranchStat(from)->GetMeanSingleSub() << "\",\n";
				os << tab << "\"DoubleSub\": \"" << GetBranchStat(from)->GetMeanDoubleSub() << "\",\n";
				os << tab << "\"TrippleSub\": \"" << GetBranchStat(from)->GetMeanTrippleSub() << "\",\n";
				os << tab << "\"dN\": \"" << GetBranchStat(from)->GetMeanNonSyn() << "\",\n";
				os << tab << "\"dS\": \"" << GetBranchStat(from)->GetMeanSyn() << "\",\n";
				os << tab << "\"dNdS\": \"" << GetBranchStat(from)->GetMeanNonSyn() / GetBranchStat(from)->GetMeanSyn() << "\",\n";
				os << tab << "\"mapping\": \"" << GetBranchName(from).c_str() << "\",\n";
				os << tab << "\"name\": \"" << GetBranchName(from).substr(0, 3);
				if (from->isLeaf()) {
					os << '_' << from->GetNode()->GetName();
				}
				os << '\"';
		} else {
				os << tab << "\"length\": \"0\",\n";
				os << tab << "\"nbSub\": \"na\",\n";
				os << tab << "\"name\": \"" << from->GetNode()->GetName() << '\"';
		}

		if (from->isLeaf()) {
				os << '\n';
		} else {
				os << ",\n";
				os << tab << "\"children\":\n";
				os << tab << '[';
				for (const Link* link = from->Next(); link != from; link = link->Next()) {
					os << "{\n";
					JsonToStream(link->Out(), tab, os);
					if (link->Next() != from) {
						os << tab << "},\n" << tab;
					} else {
						os << tab << "}]\n";
					}
				}
		}

	}

};

int main(int argc, char* argv[]) {

	string mapname = argv[1];

	ReadMapTree* tree = new ReadMapTree(mapname);

	ifstream is(mapname.c_str());
	if (!is) {
		cout << "cannot find file : " << mapname << '\n';
	}

	string line;
	while (!is.eof()) {
		is >> line;
		tree->Rereadtree(tree->GetRoot(), line);
		is >> line;
		is >> line;

	}

	tree->RecursiveNormalise(tree->GetRoot());

	ofstream os((mapname+".json").c_str());
	tree->JsonToStream(os);
	os.close();



	double dN = tree->GetMeanNonSyn();
	double dS = tree->GetMeanSyn();
	double omega = -1;
	if (dS != 0) {
		omega = dN / dS;
	}

	cout << "sub:\t" << tree->GetMeanSub() << '\n';
	cout << "dS:\t" << dS << '\n';
	cout << "dN:\t" << dN << '\n';
	cout << "omega:\t" << omega << '\n';
	cout << "Single:\t" << tree->GetMeanSingleSub() << '\n';
	cout << "Double:\t" << tree->GetMeanDoubleSub() << '\n';
	cout << "Trippl:\t" << tree->GetMeanTrippleSub() << '\n';




}
