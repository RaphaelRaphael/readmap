
#include "ReadMapTree.h"






int main(int argc, char* argv[]) {

	string mapname = argv[1];

	ReadMapTreeAA* tree = new ReadMapTreeAA(mapname);

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



	cout << "sub:\t" << tree->GetMeanSub() << '\n';




}
