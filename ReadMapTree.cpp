

#include "ReadMapTree.h"

string ReadMapTreeAA::Rereadtree(Link * from, string input) {

	// Position of the last character of the mappingM
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

void ReadMapTreeAA::JsonToStream(Link* from, string tab, ostream& os) {
	tab += ' ';
	tab += ' ';
	if (!from->isRoot()) {
			os << tab << "\"length\": \"" << GetBranchStat(from)->GetBranchLength() << "\",\n";
			os << tab << "\"Sub\": \"" << GetBranchStat(from)->GetMeanSub() << "\",\n";
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




double ReadMapTreeAA::RecursiveAddValue(Link* from, double(BranchStat::*ptr)()  ){
	double value = 0;
	for (Link* link = from->Next(); link != from; link = link->Next()) {
			value += RecursiveAddValue(link->Out(), ptr);
	}
	if (!from->isRoot()) {
			value += (GetBranchStat(from)->*ptr)();
	}
	return value;
}
