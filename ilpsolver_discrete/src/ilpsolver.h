#ifndef __ILP_H__
#define __ILP_H__

#include "gurobi_c++.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

typedef map<string, double> MSD;
typedef pair<string, double> PSD;

class ilpsolver
{
public:
	ilpsolver(const string &alphabet_file, const string &spectrum_file);
	~ilpsolver();

public:
	// instance
	int slots;
	double ubound;
	double max_weight;
	vector<double> spectrum;
	vector<string> aa_list;
	vector<double> aa_mass;

	// members for ILP
	GRBModel * model;
	GRBEnv * env;

	vector< vector<GRBVar> > xvars;		// amino acid variables
	vector< vector<GRBVar> > lvars;		// lower endpoints variables
	vector< vector<GRBVar> > uvars;		// upper endpoints variables
	vector< vector< vector<GRBVar> > > mvars;		// region variables (M replaces U and L)
	vector< vector<GRBVar> > rvars;		// range variables
	vector< vector<GRBVar> > ovars;		// map location  variables
	vector< vector< vector<GRBVar> > > svars;		// map location acit variables
	vector<GRBVar> evars;				// error variables

	// results
	vector<int> xassign;
	vector<int> lassign;
	vector<int> uassign;
	vector<double> wassign;
	vector<double> eassign;
	vector< vector<int> > sjassign;
	vector< vector<int> > siassign;
	vector< vector<int> > massign;
	
public:
	int solve();
	int print();
	int write(const string &file);
    int reset();
    int greedy_warm_start();
    int graph_greedy_warm_start();
    
private:
	// read input files
	int read_alphabet(const string &file);
	int read_spectrum(const string &file);
	int compute_upper_bound();

	// add variables
	int add_amino_acid_variables();
	int add_lower_endpoints_variables();
	int add_upper_endpoints_variables();
	int add_range_variables();
	int add_error_variables();
	int add_set_map_variables();
        int add_set_acid_map_variables();
        int add_range_map_variables();
        


	// add constraints
	int add_amino_acid_constraints();
	int add_lower_endpoints_constraints();
	int add_upper_endpoints_constraints();
	int add_range_constraints();
	int add_error_constraints();
        int add_set_map_lbound_constraints();
        int add_set_map_ubound_constraints();
        int add_set_map_constraints();
        int add_set_acid_map_constraints();
        int add_set_location_map_ubound_constraints();
        int add_set_location_map_lbound_constraints();
        int add_range_map_constraints();
        int add_error_constraints_mvars();
        int add_order_constraints();
        int add_anchor();
        
    //warm start
    int set_mvars();
    int set_xvars();
    double graph_greedy_warm_start_helper(int edges, bool complete, int sp, vector< vector<double> > spectrum_sort, vector<int> assigned, vector<int> connectable, vector< vector <double> > connections, vector< vector<int> > direct_edges, vector<int> assignment, double total_error, vector<int> my_xassign, vector<int> my_lassign, vector<int> my_uassign);
    double best_total_error;
    

	// cutting planes (extra constraints)
	int add_ordering_cutting_planes();

	// set objective
	int set_objective();

	// collect results
	int collect_results();
};

#endif
