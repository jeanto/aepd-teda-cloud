
#include "msp.hpp"

using namespace std;

// following the rules of CEC 2014 real parameter competition, 
//      if the gap between the error values of the best solution 
//      found and the optimal solution was 10^{−8} or smaller,
//      the error was treated as 0
void evaluatePopulation(vector<node> &pop, double optimum) {
	for (int i = 0; i < pop.size(); i++) {
		pop[i].fitness = opt_func_cec2014(pop[i].x);
        if ((pop[i].fitness - optimum) < END_THRESHOLD) 
            pop[i].fitness = optimum; 
    }
}

void sortIndexWithQuickSort(double array[], int first, int last, int index[]) {
    double x = array[(first + last) / 2];
    int i = first;
    int j = last;
    double temp_var = 0;
    int temp_num = 0;

    while (true) {
      while (array[i] < x) i++;    
      while (x < array[j]) j--;      
      if (i >= j) break;

      temp_var = array[i];
      array[i] = array[j];
      array[j] = temp_var;

      temp_num = index[i];
      index[i] = index[j];
      index[j] = temp_num;

      i++;
      j--;
    }

    if (first < (i -1)) sortIndexWithQuickSort(array, first, i - 1, index);  
    if ((j + 1) < last) sortIndexWithQuickSort(array, j + 1, last, index);    
}

node operateCurrentToPBest1BinWithArchive(vector<node> pop,
        int target, int p_best_individual, 
        double scaling_factor, double cross_rate, 
        node_arc archive) {
  
    node child = pop[target];

    int r1, r2;
    
    // is there risk of infinite loop?
    //cout << " {3} " << pop.size();
    do {
        r1 = rand() % pop.size();
    } while (r1 == target);
    do {
        r2 = rand() % (pop.size() + archive.arc_ind_count);
    } while ((r2 == target) || (r2 == r1));

    int random_variable = rand() % DIM;

    // cout << "-----aki-----" << r1 << "---" << r2 << "---" << scaling_factor << endl;
  
    if (r2 >= pop.size()) {
        r2 -= pop.size();
        for (int i = 0; i < DIM; i++) {
            if ((rand_0_1() < cross_rate) || (i == random_variable)) {
	            child.x[i] = pop[target].x[i] + (scaling_factor * 
                            (pop[p_best_individual].x[i] - pop[target].x[i])) + 
                            (scaling_factor * (pop[r1].x[i] - archive.pop[r2].x[i]));
            }
            else {
	            child.x[i] = pop[target].x[i];
            }
        }
    }
    else {        
        for (int i = 0; i < DIM; i++) {
            if ((rand_0_1() < cross_rate) || (i == random_variable)) {
                child.x[i] = pop[target].x[i] + (scaling_factor * 
                            (pop[p_best_individual].x[i] - pop[target].x[i])) + 
                            (scaling_factor * (pop[r1].x[i] - pop[r2].x[i]));
            }
            else {
	            child.x[i] = pop[target].x[i];
            }
        }
    }
    
    for (int j = 0; j < DIM; j++) {
        if (child.x[j] < XLOW) {
            child.x[j]= (XLOW + pop[target].x[j]) / 2.0;
        }
        else if (child.x[j] > XUP) {
            child.x[j]= (XUP + pop[target].x[j]) / 2.0;
        }
    }

    // If the mutant vector violates bounds, the bound handling method is applied
    //modifySolutionWithParentMedium(child,  pop[target]);

    return child;
}

void modifySolutionWithParentMedium(node &child, node parent) {
    int l_problem_size = DIM;
    double l_min_region = XLOW;
    double l_max_region = XUP;

    for (int j = 0; j < l_problem_size; j++) {
        if (child.x[j] < l_min_region) {
            child.x[j]= (l_min_region + parent.x[j]) / 2.0;
        }
        else if (child.x[j] > l_max_region) {
            child.x[j]= (l_max_region + parent.x[j]) / 2.0;
        }
    }
}

double gauss(double mu, double sigma){
    return mu + sigma * std::sqrt(-2.0 * std::log(rand_0_1())) * std::sin(2.0 * PI * rand_0_1());
}

double cauchy_g(double mu, double gamma) {
    return mu + gamma * std::tan(PI*(rand_0_1() - 0.5));
}

// Voltar aqui
double sf_corr(double freq_const, double freq, double G_Max, double gg){
    double c = rand_0_1();
    double sf;
    if (c < 0.5)
        sf = 0.5 * (std::sin(2*PI*freq_const*PI) * ((G_Max-gg)/G_Max) + 1);
    else
        sf = 0.5 * (std::sin(2*PI*freq*gg) * (gg/G_Max) + 1);

    return sf;
}

double distance(vector<double> x, vector<double> uk){

    double sum = 0;
    double distance = 0;

    for (int j = 0; j < DIM; j++){
        sum = sum + std::pow((x[j] - uk[j]),2.0);
    }
    distance = std::sqrt(sum);  
    return distance;  
}

// double distance(node popr, node bestr){

//     double sum      = 0;
//     double distance = 0;

//     for (int j = 0; j < DIM; j++){
//         sum = sum + std::pow((popr.x[j]-bestr.x[j]),2.0);
//     }
//     distance = sqrt(sum);  
//     return distance;  
// }

node ls_process(node pop_ls, int igen, node best_point){

    node createPoint;
    for (int i = 0; i < DIM; i++){
        std::default_random_engine generator1;
        std::normal_distribution<double> distribution(best_point.x[i],(std::log(igen)/igen)*(std::fabs(pop_ls.x[i]-best_point.x[i])));        
        double value = distribution(generator1);

        std::normal_distribution<double> randn1(0,1); 

        std::default_random_engine generator2;
        std::default_random_engine generator3;
        float add = randn1(generator2)*best_point.x[i] - randn1(generator3)*pop_ls.x[i];

        value = value + add;

        createPoint.x.push_back(bound_checking(value));
    }

    return createPoint;
}

double bound_checking(double x){
    if (x > XUP || x < XLOW){
        x = (XUP - XLOW) * rand_0_1() + XLOW; 
    }
    return x;
}