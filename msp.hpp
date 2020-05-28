#ifndef MSP_HPP
#define MSP_HPP

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <time.h>
#include <math.h>
#include <fstream>
#include <iomanip> 
#include <random>
#include <numeric>
#include <cstdlib>
#include "fun.hpp"
//#include "lshade.hpp"

// DE parameters
extern double F; 	// differential factor for weighting the difference vector: mutation rate
extern double CR;  // crossover rate

// Distributed algorithm parameters
extern int MASTER;    			// master node
extern double MIGRATION_PROB; 	// migration probability

// General variables
extern double END_THRESHOLD;
extern int DIM;
extern int SUBPOP_SIZE;
extern int NUM_NFE;   		// number of function evaluations
extern int FUN;           	// 1: SPHERE, 2: ELLIPTIC, 3: RASTRIGIN.

// AEPD: Auto-Enhanced Population Diversity (AEPD)
extern double T; 		// 
extern double c; 		// the population will be also rediversified with a small probability in any dimension j where rjg = 1
extern int UN; 

extern int GEN;
extern double XLOW;
extern double XUP;

extern double GenMaxSelectedTeda;
extern int popsize_LS;
extern int popsize_teda;

/* CEC 2014 */
//extern double *OShift, *M, *y, *z, *x_bound;
//extern int ini_flag, n_flag, func_flag, *SS;
//extern double INF, EPS, E, PI;


// struct convergence 
struct conv {
	std::vector<double> m;	 		// mean vector
	std::vector<double> std;	 	// std vector
	std::vector<double> MR;	 		// mean vector, it is the value of m just before the last diversity
										// enhancement operation
	std::vector<int> gamma;	 		// gamma is set to 1 to indicate that the population has converged 
										// in the jth dimention
	std::vector<double> theta; 		// |m - MR|*T if std <= T; T otherwise
	std::vector<double> omega; 		// if std is not greater than omega, gamma is set to 1. min(T, theta)
};

// struct stagnation
struct stag {
	std::vector<double> m_minus1;    	// mean vector (G - 1)
	std::vector<double> std_minus1;  	// std vector (G - 1)
	std::vector<int> lambda;        	// number of successive generations where the values of mjg 
											// and stdjg remain unchanged
	std::vector<int> gamma;   	    	// denote whether the population has stagnated in the jth 
											// dimension at the gth generation
};

// struct enhancement
struct enhan {
	std::vector<int> gamma;	// flag to denote whether the population needs to be rediversified in the jth dimension at the gth generation
	int zg;				    // flag to denote whether the population is rediversified at the gth generation
	int mig;				// flat ot denote whether it is necessary to make migration
};

struct node {
	int node_id;
	std::vector<double> x;
	double fitness;
	double distance;
};

struct node_arc {
	int NP;					// number of individuals
	int arc_ind_count;		// number of filled positions
	std::vector<node> pop;	// archived individuals
};

struct tedacloud {
    int id;
    double sk;          		// quantidade de amostras;
    double vk;          		// variância;
    std::vector<double> uk;  	// centro ou média;
    std::vector<node> xk;  		// amostras que pertencem ao microgrupo
	node best_teda;
};

struct LSresult{
	std::vector<double> solution;
  	double fitness;
  	unsigned evals;
};


// DE and MSP functions
node rand_node(int node_id);
std::vector<node> mutation(std::vector<node> input);
std::vector<node> crossover(std::vector<node> v, std::vector<node> x);
std::vector<node> selection(std::vector<node> u, std::vector<node> x);
std::vector<double> sum_vector(std::vector<double> a, std::vector<double> b);
std::vector<double> diff_vector(std::vector<double> a, std::vector<double> b);
std::vector<double> mult_vector(double s, std::vector<double> v);
node get_best_node(std::vector<node>, bool maximum);
double rand_0_1();
double rand_interval();

// Auto-enhanced Population Diversity
void convergence(std::vector<node> &pop, conv &convi, stag &stagi, int &ger);
void has_converged(conv &convi);
void print_conv(int ger, conv convi, stag stagi, std::vector<node> pop);

// Stagnation
void stagnation(stag &stagi, conv &convi, int ger);
void has_stagnant(stag &stagi);
void print_stag(stag &stagi);

/* Rediversification */
// Check if it is necessary to rediversify
void is_rediversify(conv &convi, stag stagi, enhan &enhani, int nfe);

// Do rediversification based on (18) as to mjg - lowjg turn into (best - mjg) - lowjg
int do_rediversification(std::vector<node> &pop, std::vector<double> xr, enhan enhani, conv convi, node best, int nfe);

// Do rediversification using migration where one random individuo is replaced by best individuo
void do_migration(std::vector<node> &pop, std::vector<double> xr, node best);

// Do rediversification using (18) as to mjg - lowjg
int inland_rediversification(std::vector<node> &pop, node best, conv &convi, enhan enhani, int nfe);

// Productive diversification migration
int do_productive_diversity(std::vector<double> &aepd_s, conv aepd_i, std::vector<double> &zp);

// Generate new individuo based on (18)
double newj(double mjg, int nfe);


// LocalSearch
tedacloud local_search_teda(tedacloud teda, node ind, int igen, 
						int optimum, int &nfe, int maxevals, 
						std::string method);

LSresult mts_ls1_improve_dim(std::vector<double> sol, 
							double best_fitness, unsigned i, 
							std::vector<double> SR);

node mts_ls1(node ind, unsigned maxevals);


/* Rediversification */

// Optimization functions
double opt_func_cec2014(std::vector<double> input);

/* LSHADE functions */
void evaluatePopulation(std::vector<node> &pop, double optimum);
void sortIndexWithQuickSort(double array[], int first, int last, int index[]);
double gauss(double mu, double sigma);
double cauchy_g(double mu, double gamma);
node operateCurrentToPBest1BinWithArchive(std::vector<node> pop, 
        int target, int p_best_individual, double scaling_factor, 
        double cross_rate, node_arc archive);
double sf_corr(double freq_const, double freq, double G_Max, double gg);
//double distance(node popr, node bestr);
double distance(std::vector<double> x, std::vector<double> uk);
node ls_process(node pop_ls, int igen, node best_point);
double bound_checking(double x);

void modifySolutionWithParentMedium(node &child, node parent);

#endif