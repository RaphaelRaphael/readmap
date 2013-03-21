        
#include <iostream>
#include<sstream>
#include <cstdlib>
using namespace std;

#include "CodonStateSpace.h"
// #include "Exception.h"

CodonStateSpace::CodonStateSpace(GeneticCodeType type)	{
	
	nucstatespace = new DNAStateSpace;
	protstatespace = new ProteinStateSpace;

	code = type;
	if (code == Universal)	{
		Nstate = Ncodon - UniNStopCodons;

		CodonCodeWithStops = new int[Ncodon]; 	// stops included
		CodonCode = new int[Nstate]; 		// stops excluded

		CodonPos = new int*[Npos];
		for (int pos=0; pos<Npos; pos++)	{
			CodonPos[pos] = new int[Nstate]; // stops excluded
		}

		int k = 0;
		for (int i=0; i<Ncodon; i++)	{
			CodonCodeWithStops[i] = UniCodonCode[i];
			if (CodonCodeWithStops[i] != -1)	{
				CodonCode[k] = CodonCodeWithStops[i];
				for (int pos=0; pos<Npos; pos++)	{
					CodonPos[pos][k] = codonpos[pos][i];
				}
				k++;
			}
		}
	}
	else if (code == MtInv)	{
		Nstate = Ncodon - MtInvNStopCodons;

		CodonCodeWithStops = new int[Ncodon]; 	// stops included
		CodonCode = new int[Nstate]; 		// stops excluded

		CodonPos = new int*[Npos];
		for (int pos=0; pos<Npos; pos++)	{
			CodonPos[pos] = new int[Nstate]; // stops excluded
		}

		int k = 0;
		for (int i=0; i<Ncodon; i++)	{
			CodonCodeWithStops[i] = MtInvCodonCode[i];
			if (CodonCodeWithStops[i] != -1)	{
				CodonCode[k] = CodonCodeWithStops[i];
				for (int pos=0; pos<Npos; pos++)	{
					CodonPos[pos][k] = codonpos[pos][i];
				}
				k++;
			}
		}
	}
	else if (code == MtMam)	{
		Nstate = Ncodon - MtMamNStopCodons;

		CodonCodeWithStops = new int[Ncodon]; 	// stops included
		CodonCode = new int[Nstate]; 		// stops excluded

		CodonPos = new int*[Npos];
		for (int pos=0; pos<Npos; pos++)	{
			CodonPos[pos] = new int[Nstate]; // stops excluded
		}

		int k = 0;
		for (int i=0; i<Ncodon; i++)	{
			CodonCodeWithStops[i] = MtMamCodonCode[i];
			if (CodonCodeWithStops[i] != -1)	{
				CodonCode[k] = CodonCodeWithStops[i];
				for (int pos=0; pos<Npos; pos++)	{
					CodonPos[pos][k] = codonpos[pos][i];
				}
				k++;
			}
		}
	}
	/*
	else if (code == MtProt)	{

		Nstate = Ncodon - MtProtNStopCodons;

		CodonCodeWithStops = new int[Ncodon]; 	// stops included
		CodonCode = new int[Nstate]; 		// stops excluded

		CodonPos = new int*[Npos];
		for (int pos=0; pos<Npos; pos++)	{
			CodonPos[pos] = new int[Nstate]; // stops excluded
		}

		int k = 0;
		for (int i=0; i<Ncodon; i++)	{
			CodonCodeWithStops[i] = MtProtCodonCode[i];
			if (CodonCodeWithStops[i] != -1)	{
				CodonCode[k] = CodonCodeWithStops[i];
				for (int pos=0; pos<Npos; pos++)	{
					CodonPos[pos][k] = codonpos[pos][i];
				}
				k++;
			}
		}

	}
	else if (code == MtEch)	{

		Nstate = Ncodon - MtEchNStopCodons;

		CodonCodeWithStops = new int[Ncodon]; 	// stops included
		CodonCode = new int[Nstate]; 		// stops excluded

		CodonPos = new int*[Npos];
		for (int pos=0; pos<Npos; pos++)	{
			CodonPos[pos] = new int[Nstate]; // stops excluded
		}

		int k = 0;
		for (int i=0; i<Ncodon; i++)	{
			CodonCodeWithStops[i] = MtEchCodonCode[i];
			if (CodonCodeWithStops[i] != -1)	{
				CodonCode[k] = CodonCodeWithStops[i];
				for (int pos=0; pos<Npos; pos++)	{
					CodonPos[pos][k] = codonpos[pos][i];
				}
				k++;
			}
		}

	}
	*/
	else 	{
		cerr << "genetic code not recognised\n";
		cerr << type << '\n';
		exit(1);
	}
	
}

CodonStateSpace::~CodonStateSpace()	{

	delete[] CodonCode;
	delete[] CodonCodeWithStops;
	for (int pos=0; pos<Npos; pos++)	{
		delete[] CodonPos[pos];
	}
	delete[] CodonPos;

	delete nucstatespace;
	delete protstatespace;
}

string CodonStateSpace::GetState(int codon)	{
	ostringstream s;
	if (codon == -1)	{
		s << "---";
	}
	else	{
		s << DNAletters[GetCodonPosition(0,codon)] << DNAletters[GetCodonPosition(1,codon)] << DNAletters[GetCodonPosition(2,codon)];
	}
	if (s.str().length() != 3)	{
		cerr << "error in translation\n";
		exit(1);
	}
	return s.str();	
}

int CodonStateSpace::GetState(string word)	{
	return GetCodonFromDNA(GetDNAStateSpace()->GetState(word.substr(0,1)),GetDNAStateSpace()->GetState(word.substr(1,1)),GetDNAStateSpace()->GetState(word.substr(2,1)));
}

int CodonStateSpace::GetCodonFromDNA(int pos1, int pos2, int pos3)	{
	if ((pos1 == unknown) || (pos2 == unknown) || (pos3 == unknown))	{
		return unknown;
	}
	int l = 0;
	while ((l<GetNstate()) && ((pos1 != GetCodonPosition(0,l)) || (pos2 != GetCodonPosition(1,l)) || (pos3 != GetCodonPosition(2,l))))	{
		l++;
	}
	if (l == GetNstate())	{
		cerr << "error in CodonStateSpace::GetCodonFromDNA : out of bound : " << pos1 << '\t' << pos2 << '\t' << pos3 << '\n';
		if (code == Universal)	{
			cerr << "universal\n";
		}
		else if (code == MtMam)	{
			cerr << "mt mam\n";
		}
		else if (code == MtInv)	{
			cerr << "mt inv\n";
		}
		return -1;
		// throw;
		// throw new Exception;
	}
	return l;
	
}

int CodonStateSpace::GetDifferingPosition(int i, int j)	{

	// identical
	if ((GetCodonPosition(0,i) == GetCodonPosition(0,j)) && (GetCodonPosition(1,i) == GetCodonPosition(1,j)) && (GetCodonPosition(2,i) == GetCodonPosition(2,j)))	{
		return -1;
	}
	if (GetCodonPosition(0,i) != GetCodonPosition(0,j))	{
		if ((GetCodonPosition(1,i) == GetCodonPosition(1,j)) && (GetCodonPosition(2,i) == GetCodonPosition(2,j)))	{
			return 0;
		}
		else	{
			return 3;
		}
	}
	if (GetCodonPosition(1,i) != GetCodonPosition(1,j))	{
		if ((GetCodonPosition(0,i) == GetCodonPosition(0,j)) && (GetCodonPosition(2,i) == GetCodonPosition(2,j)))	{
			return 1;
		}
		else	{
			return 3;
		}
	}
	if (GetCodonPosition(2,i) != GetCodonPosition(2,j))	{
		if ((GetCodonPosition(1,i) == GetCodonPosition(1,j)) && (GetCodonPosition(0,i) == GetCodonPosition(0,j)))	{
			return 2;
		}
		else	{
			return 3;
		}
	}
	return 3;
}

/*
int CodonStateSpace::GetStateWithStops(string word)	{
	return GetCodonFromDNAWithStops(GetDNAStateSpace()->GetState(word.substr(0,1)),GetDNAStateSpace()->GetState(word.substr(1,1)),GetDNAStateSpace()->GetState(word.substr(2,1)));
}

int CodonStateSpace::GetCodonFromDNAWithStops(int pos1, int pos2, int pos3)	{
	if ((pos1 == unknown) || (pos2 == unknown) || (pos3 == unknown))	{
		return unknown;
	}
	int l = 0;
	while ((l<Ncodon) && ((pos1 != codonpos[0][l]) || (pos2 != codonpos[1][l]) || (pos3 != codonpos[2][l])))	{
		l++;
	}
	if (l == Ncodon)	{
		cerr << "error in CodonStateSpace::GetCodonFromDNAWithStops : out of bound : " << pos1 << '\t' << pos2 << '\t' << pos3 << '\n';
		exit(1);
	}
	return l;
}
string GetStateWithStops(int codon)	{

	if (codon == unknown)	{
		return "?";
	}
	int aa = TranslationWithStops(codon);
	if (aa == -1)	{
		return "O";
	}
	return GetProteinStateSpace()->GetState(aa);

}

string CodonStateSpace:TranslateDNASequenceWithStops(string dnaseq)	{

	if (dnaseq.length() % 3)	{
		cerr << "error in CodonStateSpace::Translate: dna sequence not multiple of three\n";
		cerr << "length is : " << dnaseq.length() << '\n';
		exit(1);
	}
	int N = dnaseq.length() / 3;

	ostringstream s;
	for (int i=0; i<N; i++)	{
		s << GetStateWithStops(TranslationWithStops(GetStateWithStops(dnaseq.substr(3*i,3))));
	}
	return s.str();
}
*/
		
		

