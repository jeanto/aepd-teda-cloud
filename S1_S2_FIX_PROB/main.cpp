
/*
  AEPD-BEST
	Island model. Send the best individual
  
  	Implemented by C++ on Real-Parameter Single 
	  Objective Optimization at CEC-2014.

  Version: 1.0   Date: 05/Jun/2020
  Written by Jean Nunes (jean.to[at]gmail.com)
*/

#include "msp.hpp"
#include <mpi.h>

using namespace std;

// DE standard parameters
double F  = 0.516; // differential factor for weighting the 
//						difference vector: mutation rate
double CR = 0.9;   // crossover rate

// Distributed algorithm parameters
int MASTER           	= 0;    		// master node

// General variables
double END_THRESHOLD 	= 1e-8;			// 1e-8; minimum error to exit the program
int GEN_THRESHOLD 		= 2;			// minimum threshold of generation to exit the program

int DIM;					 // problem dimension
int	SUBPOP_SIZE;			  // l-shade inits with 18*DIM
int	NUM_NFE;				  // number of function evaluations
int FUN;					  // function number
int MET;					  // migration method
int FIXED_TIME_MIG = 100;	  // as Apolloni (2008, 2014)
double MIGRATION_PROB = 0.05; // probabilist migration

// AEPD: Auto-Enhanced Population Diversity (AEPD)
double T	= 1e-3; 		// min(T,thetajg), where T = 10^-3
double c	= 1e-3; 		// the population will be also rediversified with a small 
							// 		probability in any dimension j where rjg = 1

// Change UN when SUBPOP_SIZE changes
int UN; 	// UN = SUBPOP_SIZE, the larger population size, the more 
			//		generations it will take for the population to 
			// 		enter a stable stagnation rate.

double XLOW = -100.0;	// upper bound
double XUP  = 100.0;	// lower bound

/* CEC 2014 */
double *OShift, *M, *y, *z, *x_bound;
int ini_flag = 0, n_flag, func_flag, *SS;

double INF = 1.0e99;
double EPS = 1.0e-14;
double E   = 2.7182818284590452353602874713526625;
double PI  = 3.1415926535897932384626433832795029;

// SHADE
int memory_size;
double arc_rate;

// ENS L-SHADE 
double freq_inti 	= 0.5;
double min_pop_size = 4;
double max_pop_size;

double optimum;

// For local search
double GenMaxSelected 	  	= 250; 
double GenMaxSelectedTeda 	= 5;
int popsize_teda 			= 5;
double G_Max 				= 0;
int popsize_LS 				= 10;
int maxevals;

// For teda
const double r0 		= 0.001;

// Log string
string log_name;

void shade(int rank, vector<node> &popr, vector<node> &popr_ls, 
			vector<node> &children, node_arc &archive, node &bestr, node &bestls, 
			int &nfe, int p_num, vector<double> &memory_cr, vector<double> &memory_sf,
			vector<double> &memory_freq, int &memory_pos, int &gg, int igen, 
			int G_Max, int &counter);

void de(vector<node> &popr, conv &convi, stag &stagi, enhan &enhani, 
			node &bestr, int GEN, int &nfe, int rank, double &effective_move, 
			int run, vector<int> &is_improved, vector<double> &best_fit);
			
bool compare(node a, node b);
bool compare_by_dist(node a, node b);

void aepd_teda(vector<node> &popr, conv &convi, 
				stag &stagi, enhan &enhani, int GEN, int nfe);

int create_new_individuo(int rank, int rank_source, 
						vector<double> ind, vector<double> pd);

void check_finish(int &GEN, int &exit_signal, int nfe, node bestr);
void save_log(int fun, int rank, double best, double error, double time, 
				int num_mig, int num_mig_new_inds, double eff_moves);
void save_1run(int fun, int rank, vector<int> gen_is_migrated, 
				vector<int> gen_is_improved, vector<int> improve_after_mig, 
				vector<double> best_fit);
void save_std(int fun, int rank, int GEN, vector<double> std);
void save_entropy(int fun, int rank, int GEN, double ent);
void save_nfes(int fun, int rank, int GEN, int nfe, double error);
double entropy_calc(vector<node> popr);

int main(int argc, char** argv) {

	FUN  	 	= atoi(argv[1]); // number of bench function
	int run  	= atoi(argv[2]); // 0 -> save generations that migrations is made
	DIM 	 	= atoi(argv[3]); // dimension
	MET 		= atoi(argv[4]); // method: 0 - aepd_best; 1 - aepd_teda_best; 2 - fix; 3 - prob
	NUM_NFE  	= DIM * 1e+4;
	SUBPOP_SIZE = 200; //18 * DIM;
	UN 			= SUBPOP_SIZE;		

	if (MET == 0) log_name = ("aepd_best");
	else if (MET == 1) log_name = ("aepd_teda_best");
	else if (MET == 2) log_name = ("fix");
	else if (MET == 3) log_name = ("prob");

	/******************* 
		SHADE variables
	********************/
	if (DIM == 10)
		G_Max = 2163;
	else if (DIM == 30)
		G_Max = 2745;
	else if (DIM == 50)
		G_Max = 3022;
	else if (DIM == 100)
		G_Max = 3401;

	cout << scientific << setprecision(8);	
	
	// Init SHADE parameters
	memory_size 		= 5;
	arc_rate 			= 1.4;
    max_pop_size 		= SUBPOP_SIZE;
    min_pop_size 		= 4.0;	
	double p_best_rate 	= 0.11;
	int arc_size 		= (int)round(SUBPOP_SIZE * arc_rate);	
	int gg 				= 0; 	// generation counter used For Sin
	int igen 			= 1; 	// generation counter used For LS
	int counter 		= 0;

	// optimum fitness
	optimum = FUN * 100.0; 

	// the contents of M_f and M_cr are all initialiezed 0.5
	vector<double> memory_sf(memory_size, 0.5);
	vector<double> memory_cr(memory_size, 0.5);
	vector<double> memory_freq(memory_size, freq_inti*1);

	// memory index counter
	int memory_pos = 0;

	vector<double> pop_sf;
	vector<double> pop_cr;
	vector<double> pop_freq;

	//for current-to-pbest/1
	int p_num;

	/************************ 
		END SHADE variables
	*************************/

	int rank, size; 		// process id and number of process

	// if first execution save which generations have improvements and migrations
	vector<int> is_migrated;
	vector<int> is_improved;
	vector<int> improve_after_mig; // accumullate when fitness improves after migration

	vector<double> best_fit;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int GEN = 1;					// current generation
	vector<node> popr;				// population of each process
	vector<node> popr_ls;			// population of local search
	vector<node> children;			// children of each process
	node_arc archive;				// archive struct
	conv conv_stats;				// convergence stats
	stag stag_stats;				// stagnation stats
	enhan enhan_stats;				// enhancement stats

	int nfe  				= 0;	// number of evaluation functions
	int exit_signal 		= 0;	// exit signal
	int num_mig  			= 0;	// number of migrations
	int num_mig_new_inds	= 0;	// number of migrations with new individuos
									//	(individuo that is different of the sent individuo)
	double effective_move 	= 0.0;	// sum the number of effective moves

	const clock_t begin_time = clock();
	srand(time(NULL) + rank);

	// compute the next node
	int next 		= (rank + 1) % size;
	int prev 		= (size + rank - 1) % size;
	
	node bestr;		// best from current run
	node bestp;		// best from previous run
	node bestls;	// best local search

	// for asyncronous communication
	MPI_Request myRequestSend[3];
	MPI_Request myRequestRecv[3];
	MPI_Request myRequest;

	int found_item;
	MPI_Request found_request = MPI_REQUEST_NULL;;

	MPI_Barrier(MPI_COMM_WORLD);

	// island keeps a broadcasting openned to receive signal to exit
	if (rank != MASTER){
		//MPI_Ibcast(&found_item, 1, MPI_INT, MASTER, MPI_COMM_WORLD, &found_request);
		MPI_Irecv(&found_item, 1, MPI_INT, MASTER, 3, MPI_COMM_WORLD, &found_request);
	}

	// save number of island that exits
	int num_island_exit = 0;

	// TEDA Cloud params
	vector<tedacloud> teda;
    int k_teda = 0;
	int found  = 0;
	int migs   = 0; // if already migrated?

	while(true){
		// end lock file config
		int flag_exit   = 0;	// flag to exit from evolution loop
		int flag_send   = 0;	// flag to check send of request
		int flag_recv   = 0;	// flat to check receive of request
		int flag_broad	= 0;	// flag to check if broadcast gets everyone
		int rank_source;		// rank of the source island
		vector<double> ind; 	// individuo from island
		vector<double> pd; 		// pd is the metric to generate
								//	productive diversity in migration

		ind.resize(DIM);
		pd.resize(DIM);
		pop_sf.resize(SUBPOP_SIZE);
		pop_cr.resize(SUBPOP_SIZE);
		pop_freq.resize(SUBPOP_SIZE);
				
		//cout << "[" << rank << "] nova rodada com nfe = " << nfe << endl;

		if (GEN == 1){
			// population with lenght of number of individuos
			popr.resize(SUBPOP_SIZE); 
			popr_ls.resize(popsize_LS);
			children.resize(SUBPOP_SIZE);

			// initialize the main population
			int node_counter = 0;
			for (int i = 0; i < SUBPOP_SIZE; i++) {
				popr[i] = rand_node(node_counter++);
			}

			/* LSHADE - INIT */
			// evaluate the initial population's fitness values
			evaluatePopulation(popr, optimum);

			// best of first run
			bestr = get_best_node(popr, false);

			// Initialize LS population to re-start them
			node_counter = 0;
			for (int i = 0; i < popsize_LS; i++) {
				popr_ls[i] = rand_node(node_counter++);
			}
			evaluatePopulation(popr_ls, optimum);
			// best of first run
			bestls = get_best_node(popr_ls, false);

			// nfe is equal to popsize in first run
			nfe = nfe + popsize_LS;// SUBPOP_SIZE;

			archive.NP 				= arc_size;
			archive.arc_ind_count 	= 0;
			archive.pop.resize(arc_size);
		}

		// fitness from previous generatio			
		double best_costp 	= bestr.fitness;
		p_num 				= round(SUBPOP_SIZE * p_best_rate);

		shade(rank ,popr, popr_ls, children, archive, 
			bestr, bestls, nfe, p_num, memory_cr, memory_sf, 
			memory_freq, memory_pos, gg, igen, G_Max, counter);

		//de(popr, conv_stats, stag_stats, enhan_stats, bestr, GEN, 
		//	nfe, rank, effective_move, run, is_improved, best_fit);

		// best fitness until there
		double best_costr = bestr.fitness;
		best_fit.push_back(best_costr);

		if (best_costr < best_costp){
			effective_move++;
			// save generation which improvement has succeeded
			//if (run == 0) is_improved.push_back(GEN);
			if (run == 0) is_improved.push_back(nfe);
		}

		//cout << "[" << rank << "] AEPD-TEDA..." << endl;
		enhan_stats.zg = 0;
		aepd_teda(popr, conv_stats, stag_stats, enhan_stats, GEN, nfe);
		if (MET == 0 || MET == 1){
			aepd_teda(popr, conv_stats, stag_stats, enhan_stats, GEN, nfe);
		}
		else if (MET == 2){
			if (GEN % FIXED_TIME_MIG == 0){
				enhan_stats.zg = 1;
			} else{
				enhan_stats.zg = 0;
			}
		}
		else if (MET == 3){
			if (rand_0_1() < MIGRATION_PROB){
				enhan_stats.zg = 1;
			} else{
				enhan_stats.zg = 0;
			}
		}

		if (run == 0){
			save_std(FUN, rank, GEN, conv_stats.std);
		}

		MPI_Send(&enhan_stats.zg, 1, MPI_INT, next, 0, MPI_COMM_WORLD);
		MPI_Send(conv_stats.m.data(), conv_stats.m.size(), MPI_DOUBLE, next, 0, MPI_COMM_WORLD);

		int mig = 0;
		MPI_Recv(&mig, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(ind.data(), ind.size(), MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// method: 0 - aepd_best; 1 - aepd_teda_best; 2 - fixo; 3 - prob
		if (mig == 1){
			node ind_to_send;

			// Generate new individuo with productive diversity using teda
			if (MET == 1){
				int pass_away = do_productive_diversity(ind, conv_stats, pd);	
				ind_to_send.node_id = 1;	
				ind_to_send.x 		= ind;
			}
			// send the best
			else{
				ind_to_send = bestr;
			}
			
			MPI_Send(&ind_to_send.node_id, 1, MPI_INT, prev, 0, MPI_COMM_WORLD);
			MPI_Send(ind_to_send.x.data(), bestr.x.size(), MPI_DOUBLE, prev, 0, MPI_COMM_WORLD); 
		}

//		cout 	<< "[" << rank << "][" << nfe  << "] [" << bestr.fitness << "] [" 
//				<< bestr.fitness - optimum << "] ["  << enhan_stats.zg << "]" << endl;
		if (enhan_stats.zg == 1) {
			vector<double> xr;
			xr.resize(DIM);
			int node_id;

			MPI_Recv(&node_id, 1, MPI_INT, next, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(xr.data(), xr.size(), MPI_DOUBLE, next, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			//cout << "\nIsland " << rank << " gets individual from " << next << " island." << endl;
			num_mig++;

			// if first fun, save the which generations the migrations happened 
			if (run == 0) is_migrated.push_back(GEN);
			//if (run == 0) is_migrated.push_back(nfe);

			// Avoid the best individual is replaced
			int id = rand() % popr.size();

			while (bestr.node_id == popr[id].node_id) {
				id = rand() % popr.size();
			}

			// replace with best individual of the island in front 
			popr[id].x = xr;

			node xr_fit;
			xr_fit.x 		= xr;
			xr_fit.fitness 	= opt_func_cec2014(xr);

			if (xr_fit.fitness < bestr.fitness) {
				bestr = xr_fit;
				if (run == 0){
					improve_after_mig.push_back(1);
				}
			}
			else{
				if (run == 0){
					improve_after_mig.push_back(0);
				}
			}
		}

		// check if stop criterion has been reached
		check_finish(GEN, flag_exit, nfe, bestr);

		// exit signal to exit of the loop
		MPI_Bcast(&flag_exit, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if (flag_exit == 1){
			// a stop criterian has been reached, the process has terminated
			cout << "[" << rank << "] saindo!" << endl;
			break;
		}

		if (run == 0){
			// Calculate entropy
			//double ent = entropy_calc(popr);
			//save_entropy(FUN, rank, GEN, ent);
			//save_nfes(FUN, rank, GEN, nfe, bestr.fitness);
		}
	}

	cout << "[" << rank << "] saindo..." << endl;
	MPI_Barrier(MPI_COMM_WORLD);
			
	const clock_t end_time = clock();
	//if (rank != MASTER){
	// save infos to first run
      	//if (run == 0) save_1run(FUN, rank, is_migrated, is_improved, improve_after_mig, best_fit);
	double timei 		= double(end_time - begin_time) / CLOCKS_PER_SEC;
	double eff_moves 	= effective_move / GEN;
	double errori 		= bestr.fitness - optimum; 	// error
	//save_log(FUN, rank, bestr.fitness, errori, timei, num_mig, num_mig_new_inds, eff_moves);
	//}

	MPI_Finalize();
	return 0;
}

// ~1 the largest uncertainty (high diversity index)
// ~0 no uncertainty
double entropy_calc(vector<node> popr){
	vector<double> fitness;		// save frequency of fitness
	for (int en = 0; en < popr.size(); en++){
		int e = popr[en].fitness;
		fitness.push_back(e);
	}

	vector<double> repf;	
	
	// fitness
	vector<int> fit_freq;
	for (int en = 0; en < fitness.size(); en++){
		if (!count(repf.begin(), repf.end(), fitness[en]))
			fit_freq.push_back(count(fitness.begin(), fitness.end(), fitness[en]));

		repf.push_back(fitness[en]);
	}

	double entropy_fitness = 0.0;
	for (int en = 0; en < fit_freq.size(); en++){
		double freq 	= (double)fit_freq[en]/fitness.size();
		entropy_fitness = entropy_fitness + ((freq * log(freq))/log(fitness.size()));
	}

	entropy_fitness = -entropy_fitness;

	return entropy_fitness;
}

// Function for comparing two nodes
bool compare(node a, node b) {   
    // get best nodes
    return a.fitness < b.fitness;
} 

// Function for comparing two nodes
bool compare_by_dist(node a, node b) {   
    // get best nodes
    return a.distance < b.distance;
} 

void shade(int rank, vector<node> &popr, vector<node> &popr_ls, 
			vector<node> &children, node_arc &archive, node &bestr, node &bestls, 
			int &nfe, int p_num, vector<double> &memory_cr, vector<double> &memory_sf, 
			vector<double> &memory_freq, int &memory_pos, int &gg, int igen, 
			int G_max, int &counter){

	vector<double> pop_sf;
	vector<double> pop_cr;
	vector<double> pop_freq;
	pop_sf.resize(SUBPOP_SIZE);
	pop_cr.resize(SUBPOP_SIZE);
	pop_freq.resize(SUBPOP_SIZE);	

	gg = gg + 1;
	sort(popr.begin(), popr.end(), compare);
	//cout << popr.size() << " - " << SUBPOP_SIZE << endl;
	//cout << " {1} " << popr.size() << " - " << SUBPOP_SIZE;
    for (int target = 0; target < popr.size(); target++) {

		// In each generation, CR_i and F_i used by 
		//	each individual x_i are generated by  
		//	first selecting an index r_i randomly from [1, H] 
		int random_selected_period = rand() % memory_size;
		double mu_sf 	= memory_sf[random_selected_period];
		double mu_cr 	= memory_cr[random_selected_period];
		double mu_freq 	= memory_freq[random_selected_period];

      	//generate CR_i and repair its value
      	if (mu_cr == -1) {
			pop_cr[target] = 0;
      	}
      	else {
			pop_cr[target] = gauss(mu_cr, 0.1);
			if (pop_cr[target] > 1) pop_cr[target] = 1;
			else if (pop_cr[target] < 0) pop_cr[target] = 0;	
      	}

      	// generate F_i and repair its value
		int cont1 = 0;
      	do {
			cont1++;
			pop_sf[target] = cauchy_g(mu_sf, 0.1);
			if (cont1 > 1000) cout << "{1} " << cont1;
      	} while (pop_sf[target] <= 0);

      	// generate F_i and repair its value
		int cont2 = 0;
      	do {	
			cont2++;
			pop_freq[target] = cauchy_g(mu_freq, 0.1);
			if (cont2 > 1000) cout << "{2} " << cont2;
      	} while (pop_freq[target] <= 0);		

      	if (pop_sf[target] > 1) pop_sf[target] = 1;
		if (pop_freq[target] > 1) pop_freq[target] = 1;

		vector<double> sf_freqi(popr.size(), 0);
      	if(nfe <= NUM_NFE/2)
			pop_sf[target] = sf_corr(freq_inti, pop_freq[target], G_Max, gg);

      	// p-best individual is randomly selected from 
		//  	the top pop_size *  p_i members
		int p_best_ind = rand() % p_num;

      	children[target] = operateCurrentToPBest1BinWithArchive(popr, target, 
		  	p_best_ind, pop_sf[target], pop_cr[target], archive);
    }
	//cout << " {2} passou 1.";
	
    // evaluate the children's fitness values
    evaluatePopulation(children, optimum);

	// nfe is equal to popsize in first run
	nfe = nfe + SUBPOP_SIZE;

	// update the best
	node bestg = get_best_node(children, false);
	if (bestg.fitness < bestr.fitness){
		bestr = bestg;
	}

  	// update the bsf-solution and check the current number 
	//  	of fitness evaluations if the current number of 
	// 		fitness evaluations over the max number of fitness 
	// 		evaluations, the search is terminated. So, this 
	// 		program is unconcerned about L-SHADE algorithm directly
	vector<double> dif, goodCR, goodSF, goodFreq, dif_val;
	for (int i = 0; i < popr.size(); i++) {
		double difi = fabs(popr[i].fitness - children[i].fitness);
		dif.push_back(difi);
		if (popr[i].fitness > children[i].fitness) {
			goodCR.push_back(pop_cr[i]);
			goodSF.push_back(pop_sf[i]);
			goodFreq.push_back(pop_freq[i]);
			dif_val.push_back(difi);

			int flag_dup = 0;
			// Method 2: Dont includes duplicate elements
			for (int j = 0; j < archive.pop.size(); j++){
				if (archive.pop[j].fitness == popr[i].fitness){
					flag_dup = 1;
					break;
				}
			}	
			if (flag_dup == 0){
				if (archive.arc_ind_count < archive.NP) {
					archive.pop[archive.arc_ind_count] = popr[i];
					archive.arc_ind_count++;
				}
				//Whenever the size of the archive exceeds, 
				//  randomly selected elements are deleted to
				// 	make space for the newly inserted elements
				else {
					int random_selected_arc_ind = rand() % archive.NP;
					archive.pop[random_selected_arc_ind] = popr[i];	    
				}				
			}		
		}
	}
	pop_cr.clear();
	pop_sf.clear();
	pop_freq.clear();

    // generation alternation
    for (int i = 0; i < popr.size(); i++) {
		if (children[i].fitness < popr[i].fitness) {
			popr[i] = children[i];
		}
	}

	//popr_old = popr;
	
    int num_success_params = goodCR.size();

    // if numeber of successful parameters > 0, 
	//  	historical memories are updated 
    if (num_success_params > 0) {      
      	memory_sf[memory_pos] 	= 0;
      	memory_cr[memory_pos] 	= 0;
      	double temp_sum_sf 		= 0;
      	double temp_sum_cr 		= 0;
		double temp_sum_freq 	= 0;
      	double sum 				= 0;
		double weight;
      
		for (auto& n : dif_val)
			sum += n;

    	// weighted lehmer mean
      	for (int i = 0; i < num_success_params; i++) {
			weight = dif_val[i] / sum;

			memory_sf[memory_pos] += weight * goodSF[i] * goodSF[i];
			temp_sum_sf += weight * goodSF[i];

			memory_cr[memory_pos] += weight * goodCR[i] * goodCR[i];
			temp_sum_cr += weight * goodCR[i];

			memory_freq[memory_pos] += weight * goodFreq[i] * goodFreq[i];
			temp_sum_freq += weight * goodFreq[i];
      	}

      	memory_sf[memory_pos] /= temp_sum_sf;

      	if (temp_sum_cr == 0 || memory_cr[memory_pos] == -1) 
		  	memory_cr[memory_pos] = -1;
      	else memory_cr[memory_pos] /= temp_sum_cr;

      	if (temp_sum_freq == 0 || memory_freq[memory_pos] == -1) 
		  	memory_freq[memory_pos] = -1;
      	else memory_freq[memory_pos] /= temp_sum_freq;

      	memory_pos++;
      	if (memory_pos >= memory_size) memory_pos = 0;

      	//clear out the S_F, S_CR and delta fitness
      	goodCR.clear();
      	goodSF.clear();
		goodFreq.clear();
      	dif_val.clear();
    }	

	// for resizing the population size
	if (nfe >= (NUM_NFE/2)){
		counter++;
		int plan_pop_size = round(((((min_pop_size + 1) - max_pop_size) / NUM_NFE) * nfe) + max_pop_size);
		
		int reduction_ind_num;
		if (SUBPOP_SIZE > plan_pop_size) {
			reduction_ind_num = SUBPOP_SIZE - plan_pop_size;
            if ((SUBPOP_SIZE - reduction_ind_num) < min_pop_size) 
				reduction_ind_num = SUBPOP_SIZE - min_pop_size;
			
			if (counter == 1){
				int count = 0;
				bool stop = false;
				vector<node> niche;
			
            	// Change here Niche-based reduction
            	//		Firstly exclude the best niche (Half of the individuals)
            	//		Step 1: sort according fitness to pick up the best individual			
			
            	// 		Step 2: find E-distance between best_mem and all others
            	// 			To Choose neighbourhood region to the best individual

				// 	sort(popr.begin(), popr.end(), compare);
				
				// calculate euclidean distance
				for (int i = 0; i < popr.size(); i++){
					popr[i].distance = distance(popr[i].x, bestr.x);
				}
				// Sort and chooose smallest distance to have higher diversity
				sort(popr.begin(), popr.end(), compare_by_dist);

				// Select the members of the best niche
				int best_niche_size = round(SUBPOP_SIZE/2);
				vector<node> best_niche;
				best_niche.resize(best_niche_size);

				for (int i = 0; i < best_niche_size; i++){
					best_niche[i] = popr[i];
					popr.erase(popr.begin()+i);
				}
				//popr = popr_old;

				int niche_size = 20;
				node best_mem;
				for (int i = 0; i < reduction_ind_num; i++){
					sort(popr.begin(), popr.end(), compare);
					best_mem = popr[0];

					// calculate euclidean distance
					for (int j = 0; j < popr.size(); j++){
						popr[j].distance = distance(popr[j].x, best_mem.x);
					}

					// Sort and chooose smallest distance to have higher diversity
					sort(popr.begin(), popr.end(), compare_by_dist);

					if (popr.size() < niche_size)
						niche_size = popr.size();

					niche.resize(niche_size);		

					for (int j = 0; j < niche_size; j++)
						niche[j] = popr[j];	

            		// Now remove half of them excluding the best 
					int del_vec = 0;
					for (int t = 0; t < (niche_size/2); t++){
						count++;
						del_vec++;
						if (count == reduction_ind_num) {
							stop = true;
							break;
						}
					}

					if (niche.size() > 0)
						niche.erase(niche.begin() + 1, niche.begin() + (del_vec+1)); 
					if (popr.size() > 0){
						popr.erase(popr.begin() + 1, popr.begin() + (del_vec+1));
					}
					
					if (stop == true)
						break;					
				}

				popr.insert(popr.end(), best_niche.begin(), best_niche.end());

				SUBPOP_SIZE = SUBPOP_SIZE - reduction_ind_num;
				UN 			= SUBPOP_SIZE;

			}
			// if counter is not 1
			else{
				SUBPOP_SIZE = SUBPOP_SIZE - reduction_ind_num;
				UN 			= SUBPOP_SIZE;

				for (int i = 0; i < reduction_ind_num; i++){
					sort(popr.begin(), popr.end(), compare);
					popr.erase(popr.end()-1);
				}
			}

			// update archive size
			archive.NP = round(arc_rate * SUBPOP_SIZE);
			
			if (archive.pop.size() > archive.NP){
				int tam_arc = archive.pop.size();
				//cout << " {4} " << archive.pop.size() << " " << archive.NP;
				int cont5 = 0;
				do {
					cont5++;
					if (cont5 > 1000) cout << "{5} " << cont5;
					int del = rand() % archive.pop.size();
					archive.pop.erase(archive.pop.begin()+del);	
					tam_arc--;
				} while (tam_arc > archive.NP);
			}
			if (archive.arc_ind_count > archive.NP) 
				archive.arc_ind_count = archive.NP;
		}	
	}

    // Call LS based on Gaussian works when NP 
	//  	is less than 20 for the first time
	if (SUBPOP_SIZE <= 20)
		counter = counter + 1;

	// Local Search
	bool flag_LS = false;
	if (counter == 1)
		flag_LS = true;
	else
		flag_LS = false;

	flag_LS = false;
	if (flag_LS == true){
		// Pick 10 random individuals from L-SHADE pop
		for (int gen_LS = 0; gen_LS < GenMaxSelected; gen_LS++){
			// creating new point
			vector<node> new_point;
			new_point.resize(popsize_LS);

			for (int i = 0; i < popsize_LS; i++){
				new_point[i] = ls_process(popr_ls[i], igen, bestls);
			}
			evaluatePopulation(new_point, optimum);

			for (int i = 0; i < popsize_LS; i++){
				int r_index = rand() % popr.size();
				// Update those 10 random individuals from pop L-SHADE
				if (new_point[i].fitness < popr[r_index].fitness){
					popr[r_index] = new_point[i];
				} 
				// Update best individual L-SHADE
				if (new_point[i].fitness < bestr.fitness){
					bestr = new_point[i];
				}
				nfe = nfe + 1;
			}

			sort(new_point.begin(), new_point.end(), compare);
			// first point is the best
			bestls = new_point[0];
			popr_ls = new_point;
		}
	}
	//cout << "[" << rank << "] {6} passou." << endl; 
}

void aepd_teda(vector<node> &popr, conv &convi, stag &stagi, 
				enhan &enhani, int GEN, int nfe){
	// check convergence
	convergence(popr, convi, stagi, GEN);
	has_converged(convi);
	
	// check stagnation
	stagnation(stagi, convi, GEN);
	has_stagnant(stagi);

	// rediversify with migration
	is_rediversify(convi, stagi, enhani, nfe);
}

void check_finish(int &GEN, int &exit_signal, int nfe, node bestr){
	GEN++; // increment generation

	// stop criterion has been reached -> set exit signal
	double error = bestr.fitness - optimum;		// compute error
	if (error <= END_THRESHOLD || nfe >= NUM_NFE) {
		exit_signal = 1; // set exit signal setado
	}
}

int create_new_individuo(int rank, int rank_source, vector<double> ind, vector<double> pd){

	//cout << "[" << rank << "] recebeu uma solicitacao de [" << rank_source << "]..." << endl;

	MPI_Recv(ind.data(), DIM, MPI_DOUBLE, rank_source, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(pd.data(), DIM, MPI_DOUBLE, rank_source, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	int pass_away = 1;
	int tag_send = 2;
	//MPI_Send(&rank, 1, MPI_INT, rank_source, tag_send, MPI_COMM_WORLD);
	MPI_Send(&pass_away, 1, MPI_INT, rank_source, tag_send, MPI_COMM_WORLD);
	MPI_Send(ind.data(), DIM, MPI_DOUBLE, rank_source, tag_send, MPI_COMM_WORLD);
	MPI_Send(pd.data(), DIM, MPI_DOUBLE, rank_source, tag_send, MPI_COMM_WORLD);
	
	return 0;
}

/*
	--- LOG FILE FORMAT ----

	dimension  function  rank  best_fitness  error  time  
	num_mig		num_mig_new_inds	eff_moves  
*/
void save_log(int fun, int rank, double best, double error, double time, 
				int num_mig, int num_mig_new_inds, double eff_moves){

	string filePath = "../../logs/" + log_name + "/log_" + to_string(DIM) + "_" + to_string(fun) + "_" + to_string(rank) + ".csv";

    ofstream ofs(filePath.c_str(), ios_base::out | ios_base::app);
    ofs << DIM << ';' << fun << ';' << rank << ';' 
					<< best << ';' << error << ';'
					<< time << ';' << num_mig << ';'
					<< num_mig_new_inds << ';' 
					<< eff_moves << '\n';
    ofs.close();
}

/*
	--- LOG FILE FORMAT TO FIRST RUN ----

	dimension  function  process  	\n 
	generation_with_migration		\n
	generation_with_improvement

*/
void save_1run(int fun, int rank, vector<int> gen_is_migrated, 
				vector<int> gen_is_improved, vector<int> improve_after_mig,
				vector<double> best_fit){
	string filePath = "../../logs/" + log_name + "/one_log_" + to_string(DIM) + "_" + to_string(fun) + "_" + to_string(rank) + ".csv";

	ofstream ofs(filePath.c_str(), ios_base::out | ios_base::ate);
    ofs << DIM << ';' << fun << ';' << rank << '\n' << '\n';
	for (const auto &e : gen_is_migrated) ofs << e << " ";
	ofs << '\n' << '\n';
	//for (const auto &e : gen_is_improved) ofs << e << " ";
	//ofs << '\n' << '\n';
	for (const auto &e : improve_after_mig) ofs << e << " ";
	// ofs << '\n' << '\n';
	// for (const auto &e : best_fit) ofs << e << " ";
    ofs.close();
}

void save_std(int fun, int rank, int GEN, vector<double> std){
	string filePath = "../../logs/" + log_name + "/std_log_" + to_string(DIM) + "_" + to_string(fun) + "_" + to_string(rank) + ".csv";
	ofstream ofs(filePath.c_str(), ios_base::out | ios_base::app);
    ofs << DIM << ';' << fun << ';' << rank << '\n' << '\n';
	ofs << GEN << " ";
	for (const auto &e : std) ofs << e << " ";
	ofs << '\n';
	ofs.close();
}

void save_entropy(int fun, int rank, int GEN, double ent){
	string filePath = "../../logs/" + log_name + "/entropy_log_" + to_string(DIM) + "_" + to_string(fun) + "_" + to_string(rank) + ".csv";

    ofstream ofs(filePath.c_str(), ios_base::out | ios_base::app);
    ofs << DIM << ';' << fun << ';' << ';' << rank << ';' << GEN << ';' 
		<< ent << '\n';
    ofs.close();	
}

void save_nfes(int fun, int rank, int GEN, int nfe, double error){
	string filePath = "../../logs/" + log_name + "/nfes_log_" + to_string(DIM) + "_" + to_string(fun) + "_" + to_string(rank) + ".csv";

    ofstream ofs(filePath.c_str(), ios_base::out | ios_base::app);
    ofs << DIM << ';' << fun << ';' << rank << ';' << GEN << ';' << nfe << ';' << error << '\n';
    ofs.close();	
}

void de(vector<node> &popr, conv &convi, stag &stagi, enhan &enhani, 
			node &bestr, int GEN, int &nfe, int rank, double &effective_move, 
			int run, vector<int> &is_improved, vector<double> &best_fit){

	// fitness from previous generation
	double best_costp = 0;
	if (GEN > 1) best_costp = bestr.fitness; //opt_func_cec2014(diff_vector(bestr.x, ZEROES));

	// every single process makes DE steps 
	vector<node> v = mutation(popr);			// mutation
	vector<node> u = crossover(v, popr);		// crossover
	popr  		   = selection(u, popr);		// selection
	
	// get the best individuo of each subpopulation
	bestr = get_best_node(popr, false);

	nfe = nfe + SUBPOP_SIZE;
}
