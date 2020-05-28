
#include "msp.hpp"

/*
	Check if it is necessary to do rediversification
*/

bool compare_ls(node a, node b);

void is_rediversify(conv &convi, stag stagi, enhan &enhani, int nfe){

	std::vector<int> gamma(DIM, 0);
	enhani.gamma = gamma;

	// rediversify with a small probability
	enhani.zg = 0;
	if (rand_0_1() < c){
		enhani.zg = 1;
	}

	int sum_rjg = 0;
	//int sum_mig = 0;
	for (int k = 0; k < DIM; k++){
		if (convi.gamma[k] == 1 || stagi.gamma[k] == 1){
			sum_rjg++;
			enhani.gamma[k] = 1;
		}
	}

	// rediversify if all dimensions have converged or stagnation
	//if (sum_rjg == DIM){
	//if (sum_rjg >= DIM*0.5){	
	double conv_stag = sum_rjg/(double)DIM;
	double prob_mig  = (1.0-(nfe/(double)NUM_NFE));
	//std::cout << conv_stag << " " << prob_mig << " " << NUM_NFE << std::endl;
	if (conv_stag >= prob_mig){
		enhani.zg = 1;
	}

	// make migration if sum_rjg = 0
	enhani.mig = 0;
	if (sum_rjg == 0){
	 	enhani.mig = 1;
	}
}

// generate new j according AEPD (generate new individuo based on (18))
double newj(double mjg, int nfe){

	double lowjg = std::min(mjg, XLOW);
	double upjg  = std::max(mjg, XUP);

	double meanjg = ((mjg - (lowjg)) / (upjg - (lowjg)));

	double std_hat = std::max(meanjg, 1 - meanjg);

	double a = 5e-4;
	double stdjg = exp((-a*nfe)/DIM) * std_hat;

	if (stdjg < 1e-3){
		stdjg = 1e-3;
	}

	std::default_random_engine generator;
	std::normal_distribution<double> randn(meanjg, pow(stdjg,2));

	double xij = lowjg + (upjg - lowjg) * randn(generator);

	return xij;
}

// Rediversification strategies
//		Do rediversification based on (18) as to mjg - lowjg turn into (best - mjg) - lowjg
int do_rediversification(std::vector<node> &pop, std::vector<double> xr, enhan enhani, conv convi, node best, int nfe){
	int flagM = 0;
	for (int j = 0; j < DIM; j++){
		if (enhani.gamma[j] == 1){
			flagM = 1;
			for (int i = 0; i < SUBPOP_SIZE; i++){
				// is pop[i] is not the best individuo
				if (pop[i].node_id != best.node_id) {
					pop[i].x[j] = newj((xr[j] - convi.m[j]), nfe);
				}
			}
			// it changes MR when it is necessary to diversify
			convi.MR[j] = convi.m[j]; 
		}
	}

	return flagM;
}

// Do rediversification using (18) as to mjg - lowjg
int inland_rediversification(std::vector<node> &pop, node best, conv &convi, enhan enhani, int nfe){
	// Get rediversification
	int flagM = 0;
	for (int j = 0; j < DIM; j++){
		if (enhani.gamma[j] == 1){
			flagM = 1;
			for (int i = 0; i < SUBPOP_SIZE; i++){
				// is pop[i] is not the best individuo
				if (pop[i].node_id != best.node_id) {
					pop[i].x[j] = newj(convi.m[j], nfe);
				} 
			}
			// it changes MR when it is necessary to diversify
			convi.MR[j] = convi.m[j];
		}
	}

	return flagM;
}

// Do rediversification using migration where one random individuo is replaced by best individuo
void do_migration(std::vector<node> &pop, std::vector<double> xr, node best){
	int id = rand() % pop.size();
	while (best.node_id != pop[id].node_id) {
		id = rand() % pop.size();
	}

	// replace by the best individuo sent by master 
	pop[id].x = xr;
}

/*
	Generates productive diversity
		- aepd_s -> mean of source population
		- aepd_i -> mean of current population
		- zp	 -> vector that shows which j's was diversified
					if zp[j] = 0 -> it was diversified, else pass it away
*/
int do_productive_diversity(std::vector<double> &aepd_s, conv aepd_i, std::vector<double> &zp){

	float ecc, ecc_norm, thr_out;
	int k 	= 2;
	int m 	= 1; 		// number of standard deviations
	int sum = 0; 		// sum of how much j's was diversified
	int pass_away = 0;
	for (int j = 0; j < DIM; j++){
		// ecentricidade
		ecc = (1 / k) + ((pow((aepd_i.m[j] - aepd_s[j]), 2.0)) / (k * aepd_i.std[j]));
		// ecentricidade normalizada
		ecc_norm 	= ecc / 2.0;
		thr_out 	= ((pow(m, 2.0) + 1.0) / 4.0);

		// identifica outliers atraves da desigualdade de chebyshev
		// se ecc_norm > thr_out eh outlier -> tem diversidade
		if (ecc_norm > thr_out) {
			zp[j] 		= 0;

			// what about use DE mutation
			aepd_s[j] 	= aepd_i.m[j];
			sum 		= sum + 1;     
		}
	}

	//if (sum < DIM){
	if (sum > 0){
		pass_away = 1;
	}

	// if (sum == 0){
		// Usar alguma estratégia mais agressiva para gerar diversidade (juntar as duas ilhas
		// e reiniciar alguns individuos!)
	//}
	return pass_away;
}

// Function for comparing two nodes
bool compare_ls(node a, node b) {   
    // get best nodes
    return a.fitness < b.fitness;
} 

tedacloud local_search_teda(tedacloud teda, node ind, int igen, 
						int optimum, int &nfe, int maxevals, 
						std::string method = "mts"){

	int popsize_ls_teda = popsize_teda; 
	evaluatePopulation(teda.xk, optimum);
	teda.best_teda = get_best_node(teda.xk, false);

	ind.fitness = opt_func_cec2014(ind.x);
	if (ind.fitness < teda.best_teda.fitness)
		teda.best_teda = ind;

	// population of local search
	std::vector<node> popr_ls_teda;
	for (int size_LS = 0; size_LS < teda.xk.size(); size_LS++){
		//if (size_LS < popsize_ls_teda){
			popr_ls_teda.push_back(teda.xk[size_LS]);
		//}
		//else{
		//	break;
		//}
	}
	popsize_ls_teda = popr_ls_teda.size();
	//node best_ls = get_best_node(popr_ls_teda, false);

	// creating new point
	if (method.compare("gaussian") == 0){
		for (int gen_LS = 0; gen_LS < GenMaxSelectedTeda; gen_LS++){
			for (int i = 0; i < popsize_ls_teda; i++){
				popr_ls_teda[i] = ls_process(popr_ls_teda[i], igen, teda.best_teda);
			}
			evaluatePopulation(popr_ls_teda, optimum);
		}
	}
	else if (method.compare("mts") == 0){
		for (int i = 0; i < popsize_ls_teda; i++){
			popr_ls_teda[i] = mts_ls1(popr_ls_teda[i], maxevals);

			if (popr_ls_teda[i].fitness < teda.best_teda.fitness){
				teda.best_teda = popr_ls_teda[i];
			}
		}
	}

	for (int i = 0; i < popsize_ls_teda; i++){
		//int r_index = rand() % teda.xk.size();

		// Update those 10 random individuals from pop L-SHADE
		//if (new_point[i].fitness < teda.xk[r_index].fitness){
		//if (popr_ls_teda[i].fitness < teda.xk[i].fitness){
			//teda.xk[i] = popr_ls_teda[i];
			//teda.xk.push_back(new_point[i]);
			//teda.xk[r_index] = new_point[i];
		//} 
		// Update best individual L-SHADE
		if (popr_ls_teda[i].fitness < teda.best_teda.fitness){
			teda.best_teda = popr_ls_teda[i];
			//teda.xk[r_index] = new_point[i];
		}
		nfe = nfe + 1;
	}

	teda.xk = popr_ls_teda;

	//sort(new_point.begin(), new_point.end(), compare_ls);
	// first point is the best
	//best_ls 	 = new_point[0];
	//popr_ls_teda = new_point;

	return teda;
}


LSresult mts_ls1_improve_dim(std::vector<double> sol, 
							double best_fitness, unsigned i, 
							std::vector<double> SR){
  	std::vector<double> newsol = sol;
  	newsol[i] -= SR[i];

  	double initial_fit = best_fitness;

  	newsol[i] = std::min(newsol[i], (double)XUP);
 	newsol[i] = std::max(newsol[i], (double)XLOW);

  	double fitness_newsol = opt_func_cec2014(newsol);

  	unsigned evals = 1;

  	if(fitness_newsol < best_fitness){
    	best_fitness = fitness_newsol;
    	sol = newsol;
  	} else if( fitness_newsol > best_fitness ){
		newsol[i] = sol[i];
		newsol[i] += 0.5 * SR[i];

		newsol[i] = std::min(newsol[i], (double)XUP);
		newsol[i] = std::max(newsol[i], (double)XLOW);

		fitness_newsol = opt_func_cec2014(newsol);
		evals++;

		if( fitness_newsol < best_fitness ){
			best_fitness = fitness_newsol;
			sol = newsol;
		}
	}
	//printf("%1.20E -> %1.20E\n", initial_fit, best_fitness);
	return LSresult{sol, best_fitness, evals};
}

node mts_ls1(node ind, unsigned maxevals){
	//unsigned int n_dim = ndim;
	std::vector<double> sol = ind.x;
	unsigned totalevals = 0;
	double best_fitness = ind.fitness;

	std::vector<double> SR(ind.x.size());
	for(unsigned i = 0; i < SR.size(); i++){
		//SR[i] = (XUP - XLOW) * 0.01;
		SR[i] = (XUP - XLOW) * 0.2;
	}

	LSresult result;
	LSresult current_best = {sol, best_fitness, 0};

  	std::vector<double> improvement(DIM, 0.0);

  	std::vector<double> dim_sorted(DIM);
  	std::iota(dim_sorted.begin(), dim_sorted.end(), 0);

  	double improve;
  	//warm-up
  	if(totalevals < maxevals){
		std::next_permutation(dim_sorted.begin(), dim_sorted.end());
		for(auto it = dim_sorted.begin(); it != dim_sorted.end(); it++){
			result 				= mts_ls1_improve_dim(sol, best_fitness, *it, SR);
			totalevals 			+= result.evals;
			improve 			= std::max(current_best.fitness - result.fitness, 0.0);
			improvement[*it] 	= improve;

			if(improve > 0.0){
				std::cout 	<< "{1}" << current_best.fitness << " > " 
							<< result.fitness << ", improve ~ " 
							<< improve << std::endl;
				current_best = result;
			} else {
				SR[*it] /= 2.0;
				//std::cout 	<< "{1} Sem melhoria!" << std::endl;
			}
		}
  	}

  	std::iota(dim_sorted.begin(), dim_sorted.end(), 0);
  	std::sort(dim_sorted.begin(), dim_sorted.end(), [&](unsigned i1, unsigned i2) { return improvement[i1] > improvement[i2]; });

  	int i, d = 0, next_d, next_i;
  	while(totalevals < maxevals){
    	i 			= dim_sorted[d];
    	result 		= mts_ls1_improve_dim(current_best.solution, current_best.fitness, i, SR);
    	totalevals 	+= result.evals;
    	improve 	= std::max(current_best.fitness - result.fitness, 0.0);
    	improvement[i] 	= improve;
    	next_d 			= (d+1) % DIM;
    	next_i 			= dim_sorted[next_d];

    	if( improve > 0.0 ){
      		std::cout 	<< "{2}" << current_best.fitness << " > " 
			  			<< result.fitness << ", improve ~ " 
						<< improve << std::endl;
      		current_best = result;

      		if( improvement[i] < improvement[next_i] ){
        		std::iota(dim_sorted.begin(), dim_sorted.end(), 0);
        		std::sort(dim_sorted.begin(), dim_sorted.end(), [&](unsigned i1, unsigned i2) { return improvement[i1] > improvement[i2]; });
      		}
    	} else {
      		SR[i] /= 2.0;
      		d = next_d;
 
     		if(SR[i] < 1e-16){
        		SR[i] = (XUP - XLOW) * 0.4;
				//SR[i] = (XUP - XLOW) * 1e6;
      		}
			//std::cout 	<< "{2} Sem melhoria!" << std::endl;
    	}
  	}
	node new_ind_mts;
	new_ind_mts.x 		= current_best.solution;
	new_ind_mts.fitness = current_best.fitness;
  	return new_ind_mts;
}
