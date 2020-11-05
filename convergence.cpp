
#include "msp.hpp"

/*
	Check convergation
*/
void has_converged(conv &convi){
	for (int k = 0; k < DIM; k++){
		// Theta				
		if (convi.std[k] <= T)
			convi.theta[k] = (std::fabs(convi.m[k] - (convi.MR[k])) * T);
		else
			convi.theta[k] = T;

		// Omega
		convi.omega[k] = std::min(T, convi.theta[k]);

		// Gamma
		if (convi.std[k] <= convi.omega[k])
			convi.gamma[k] = 1;
		else
			convi.gamma[k] = 0;
	}
}

void convergence(std::vector<node> &pop, conv &convi, stag &stagi, int &ger){

	std::vector<double> mjg(DIM, 0); 	// media das variaveis dos individuos na jth dimensao na gth geracao 
	std::vector<double> stdjg(DIM, 0);   // desvio padrao das variaveis dos individuos na jth dimensao na gth geracao
	std::vector<int> gamajg(DIM, 0);	// gamma is set to 1 to indicate that the population has converged in the jth dimention
	std::vector<double> thetajg(DIM, 0); // |m - MR|*T if std <= T; T otherwise
	std::vector<double> omegajg(DIM, 1); // if std is not greater than omega, gamma is set to 1. min(T, theta)

	// soma
	for (int i = 0; i < SUBPOP_SIZE; i++){
		for (int j = 0; j < DIM; j++){
			mjg[j] = mjg[j] + pop[i].x[j];
		}
	}   

	// media
	for (int k = 0; k < DIM; k++){
		mjg[k] = mjg[k] / SUBPOP_SIZE;
	}

    // desvio padrao
	for (int i = 0; i < SUBPOP_SIZE; i++){
		for (int j = 0; j < DIM; j++){
			stdjg[j] = stdjg[j] + std::pow(pop[i].x[j] - mjg[j],2);
		}
	}

	for (int k = 0; k < DIM; k++){
		stdjg[k] = std::sqrt(stdjg[k] / SUBPOP_SIZE);
	}

	// get previous mean and std
	if (ger > 1){
		stagi.m_minus1   = convi.m;
		stagi.std_minus1 = convi.std;		
	}

	convi.m 	= mjg;
	convi.std 	= stdjg;
	convi.gamma = gamajg;
	convi.theta = thetajg;
	convi.omega = omegajg;

	// the initial value of MRj is set to the mean value mjg of the initialized population
	if (ger == 1){
		convi.MR = convi.m;
	}
}

void print_conv(int ger, conv convi, stag stagi, std::vector<node> pop){
	std::cout << "----- START POP -----" << std::endl;
	for (int i = 0; i < SUBPOP_SIZE; i++){
		for (int j = 0; j < DIM; j++){
			std::cout << pop[i].x[j] << " ";
		}
		std::cout << std::endl;
	}  
	std::cout << "----- END POP -----\n" << std::endl;

	std::cout << "----- START MEDIA -----" << std::endl;
	std::cout << "m : ";
	for (int k = 0; k < DIM; k++){
		std::cout << convi.m[k] << " ";
	}
	std::cout << std::endl;
	std::cout << "MR: ";
	for (int k = 0; k < DIM; k++){
		std::cout << convi.MR[k] << " ";
	}
	std::cout << std::endl;
	std::cout << "m-: ";
	for (int k = 0; k < DIM; k++){
		std::cout << stagi.m_minus1[k] << " ";
	}
	std::cout << "\n----- END MEDIA -----\n" << std::endl;

	std::cout << "----- START STD -----" << std::endl;
	std::cout << "std : ";
	for (int k = 0; k < DIM; k++){
		std::cout << convi.std[k] << " ";
	}
	std::cout << std::endl;
	std::cout << "std-: ";
	for (int k = 0; k < DIM; k++){
		std::cout << stagi.std_minus1[k] << " ";
	}
	std::cout << "\n----- END STD -----\n" << std::endl;
	
	std::cout << "----- START THETA -----" << std::endl;
	for (int k = 0; k < DIM; k++){
		printf("%.15f ", convi.theta[k]);
	}
	std::cout << "\n----- END THETA -----\n" << std::endl;
	
	std::cout << "----- START OMEGA -----" << std::endl;
	for (int k = 0; k < DIM; k++){
		printf("%.15f ", convi.omega[k]);
	}
	std::cout << "\n----- END OMEGA -----\n" << std::endl;

	std::cout << "----- START GAMMA CONV -----" << std::endl;
	for (int k = 0; k < DIM; k++){
		std::cout << convi.gamma[k] << " "; 
	}
	std::cout << "\n----- END GAMMA CONV -----\n" << std::endl;
}