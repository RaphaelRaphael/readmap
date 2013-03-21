

#include "CodonStateSpace.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;


class CodonMappingParser {

	string map;
	int nbstates;
	int* states;
	double* times;
	CodonStateSpace* cstatespace;

	public:

	CodonMappingParser(string input, CodonStateSpace* incstatespace) {

		cstatespace =incstatespace;
		map = input;

		// On calcule le nombre d etats dans le mapping
		// Egal au nombre de ':' sur deux
		nbstates = 0;
		for(int j= 0; j < input.length(); j++) {
			if (input[j] == ':') {
				nbstates++;
			}
		}
		nbstates /= 2;

		// On cré les deux tables qui contiendront les états et les longueurs de branche
		states = new int[nbstates];
		times = new double[nbstates];

		//On remplie les deux tables
		input=input.substr(3, -1);
		for(int i = nbstates-1; i>-1; i--){
			input=input.substr(1, -1);
			int j;
			for(j=0;input[j] != ':';j++){}
			times[i] = atof(input.substr(0, j).c_str());
			states[i] = cstatespace->GetState(input.substr(j + 1, 3));
			input = input.substr(j+4, -1);
		}
	}

	~CodonMappingParser() {
		delete[] states;
		delete[] times;
	}


	// Donne le nombre de substitutions du mapping
	int GetNbSub(){
		return nbstates-1;
	}

	// Donne le nombre de substitutions synonymes du mapping
	int GetNbSynSub(){
		int res = 0;
		for (int i = 0; i < GetNbSub(); i++){
			if(cstatespace->Synonymous(states[i], states[i+1])){
				res++;
			}
		}
		return res;
	}

	// Donne le nombre de substitutions non-synonymes du mapping
	int GetNbNonSynSub(){
		return GetNbSub() - GetNbSynSub();
	}


	// Donne le nombre de substitutions simples du mapping
	int GetNbSingleSub(){
		int res = 0;
		for (int i = 0; i < GetNbSub(); i++){
			if(cstatespace->GetNbDiff(states[i], states[i+1]) ==1){
				res++;
			}
		}
		return res;
	}

	// Donne le nombre de substitutions doubles du mapping
	int GetNbDoubleSub(){
		int res = 0;
		for (int i = 0; i < GetNbSub(); i++){
			if(cstatespace->GetNbDiff(states[i], states[i+1]) == 2){
				res++;
			}
		}
		return res;
	}


	// Donne le nombre de substitutions tripples du mapping
	int GetNbTrippleSub(){
		int res = 0;
		for (int i = 0; i < GetNbSub(); i++){
			if(cstatespace->GetNbDiff(states[i], states[i+1]) ==3){
				res++;
			}
		}
		return res;
	}

	// Donne la taille de la branche
	double GetTotalLength(){
		double res = times[0];
		for (int i = 0; i < GetNbSub(); i++){
			res += times[i+1];
		}
		return res;
	}

	// Donne la taille d'un intervalle de la branche
	double GetLength(int i){
		return times[i];
	}

};



class SimpleMappingParser {

	string map;
	int nbstates;
	int* states;
	double* times;
	SimpleStateSpace* sstatespace;

	public:

	SimpleMappingParser(string input, SimpleStateSpace* insstatespace) {

		sstatespace = insstatespace;
		map = input;

		// On calcule le nombre d etats dans le mapping
		// Egal au nombre de ':' sur deux
		nbstates = 0;
		for(int j= 0; j < input.length(); j++) {
			if (input[j] == ':') {
				nbstates++;
			}
		}
		nbstates /= 2;

		// On cré les deux tables qui contiendront les états et les longueurs de branche
		states = new int[nbstates];
		times = new double[nbstates];

		//On remplie les deux tables
		input=input.substr(1, -1);
		for(int i = nbstates-1; i>-1; i--){
			input=input.substr(1, -1);
			int j;
			for(j=0;input[j] != ':';j++){}
			times[i] = atof(input.substr(0, j).c_str());
			states[i] = sstatespace->GetState(input.substr(j + 1, 1));
			input = input.substr(j+4, -1);
		}
	}

	~SimpleMappingParser() {
		delete[] states;
		delete[] times;
	}


	// Donne le nombre de substitutions du mapping
	int GetNbSub(){
		return nbstates-1;
	}


	// Donne la taille de la branche
	double GetTotalLength(){
		double res = times[0];
		for (int i = 0; i < GetNbSub(); i++){
			res += times[i+1];
		}
		return res;
	}



	// Donne la taille d'un intervalle de la branche
	double GetLength(int i){
		return times[i];
	}


};
