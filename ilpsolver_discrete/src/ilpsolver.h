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
	vector< vector<GRBVar> > rvars;		// range variables
	vector<GRBVar> evars;				// error variables

	// results
	vector<int> xassign;
	vector<int> lassign;
	vector<int> uassign;
	vector<double> wassign;
	vector<double> eassign;

public:
	int solve();
	int print();
	int write(const string &file);

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

	// add constraints
	int add_amino_acid_constraints();
	int add_lower_endpoints_constraints();
	int add_upper_endpoints_constraints();
	int add_range_constraints();
	int add_error_constraints();

	// cutting planes (extra constraints)
	int add_ordering_cutting_planes();

	// set objective
	int set_objective();

	// collect results
	int collect_results();
};

#endif
