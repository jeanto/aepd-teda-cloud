
#include "msp.hpp"

using namespace std;

vector<node> mutation(vector<node> input) {
	// escolhe tres indices randomicos
	vector<node> result = input;
	for (int i = 0; i < input.size(); i++) {
		int r1, r2, r3;
		r1 = r2 = r3 = i;
		while (r1 == i)	r1 = rand() % input.size();
		while (r2 == i)	r2 = rand() % input.size();
		while (r3 == i)	r3 = rand() % input.size();

		result[i].x = (sum_vector(input[r1].x, mult_vector(F, diff_vector(input[r2].x, input[r3].x))));
	}
	return result;
}

vector<node> crossover(vector<node> v, vector<node> x) {
	vector<node> result = x;

	for (int i = 0; i < x.size(); i++) {
		int j_rand = rand() % DIM;
		node n;
		n.node_id = x[i].node_id;
		for (int j = 0; j < DIM; j++) {
			if (rand_0_1() <= CR || j_rand == j) {
				n.x.push_back(v[i].x[j]);
			}
			else {
				n.x.push_back(x[i].x[j]);
			}
		}
		result[i] = n;
	}
	return result;
}

vector<node> selection(vector<node> u, vector<node> x) {
	vector<node> result = x;

	for (int i = 0; i < x.size(); i++) {
		u[i].fitness = opt_func_cec2014(u[i].x);
		x[i].fitness = opt_func_cec2014(x[i].x);
		if (u[i].fitness <= x[i].fitness) {
			result[i] 			= u[i];
		}
		else {
			result[i] 			= x[i];
		}
		
	}
	return result;
}

vector<double> sum_vector(vector<double> a, vector<double> b) {
	vector<double> result;
	if (a.size() != b.size()) {
		cerr << "SUM OF VECTOR IS INVALID!" << endl;
		exit(1);
	}

	for (int i = 0; i < a.size(); i++) {
		result.push_back(a[i] + b[i]);
	}
	return result;
}

vector<double> mult_vector(double s, vector<double> v) {
	vector<double> result;

	for (int i = 0; i < v.size(); i++) {
		result.push_back(s * v[i]);
	}
	return result;
}

vector<double> diff_vector(vector<double> a, vector<double> b) {
	vector<double> result;
	if (a.size() != b.size()) {
		cerr << "DIFFERENCE OF VECTOR IS INVALID! " << endl;
		exit(1);
	}

	for (int i = 0; i < a.size(); i++) {
		result.push_back(a[i] - b[i]);
	}
	return result;
}

node rand_node(int node_id) {
	node new_node;
	new_node.node_id = node_id;
    
	for (int i = 0; i < DIM; i++) {
		// generate without concerning with bounds
		//new_node.x.push_back(rand()); 

		// generate considering the bounds
		new_node.x.push_back(rand_interval());
	}

	return new_node;
}

node get_best_node(vector<node> input, bool maximum = true) {
	node best_max = input[0];
	node best_min = input[0];

	double best_value_max = input[0].fitness;
	double best_value_min = input[0].fitness;

	for (int i = 1; i < input.size(); i++) {
		double value = input[i].fitness; 
		if (value > best_value_max) {
			best_value_max 	= value;
			best_max 		= input[i];
		}

		if (value < best_value_min) {
			best_value_min 	= value;
			best_min 		= input[i];
		}
	}
	if (maximum) {
		return best_max;
	}
	else {
		return best_min;
	}
}

double opt_func_cec2014(std::vector<double> input) {

    double f;
	cec14(input.data(), &f, (int) input.size(), FUN);  

	return f;
}

double rand_interval(){
	return (XLOW + 1) + (((double) rand()) / (double) RAND_MAX) * (XUP - (XLOW + 1));
}

double rand_0_1() {
	double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	return r;
}