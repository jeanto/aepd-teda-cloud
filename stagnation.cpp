
#include "msp.hpp"

/*
	Check stagnation
*/
void has_stagnant(stag &stagi){
	for (int k = 0; k < DIM; k++){
		if (stagi.lambda[k] >= UN){
			stagi.gamma[k] = 1;
		}
		else{
			stagi.gamma[k] = 0;
		}
	}
}

void stagnation(stag &stagi, conv &convi, int ger){

	if (ger == 1){
		std::vector<int> gamajg(DIM, 0);
		std::vector<int> lambdajg(DIM, 0);
		stagi.gamma  = gamajg;
		stagi.lambda = lambdajg;
	}
	else{
		for (int k = 0; k < DIM; k++){
			if((convi.m[k] == stagi.m_minus1[k]) && (convi.std[k] == stagi.std_minus1[k])){
				stagi.lambda[k] = stagi.lambda[k] + 1;
			}
			else{
				stagi.lambda[k] = 0;
			}
		}
	}
}

void print_stag(stag &stagi){
	std::cout << "----- START LAMBDA STAG -----" << std::endl;
	for (int k = 0; k < DIM; k++){
		std::cout << stagi.lambda[k] << " ";
	}
	std::cout << "\n----- END LAMBDA STAG -----\n" << std::endl;	

	std::cout << "----- START GAMMA STAG -----" << std::endl;
	for (int k = 0; k < DIM; k++){
		std::cout << stagi.gamma[k] << " ";
	}
	std::cout << "\n----- END GAMMA STAG -----\n" << std::endl;
}