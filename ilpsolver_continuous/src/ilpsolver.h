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
	ilpsolver(const string &spectrum_file);
	~ilpsolver();

public:
	// instance
	int slots;
	double ubound;
	vector<double> spectrum;

	// members for ILP
	GRBModel * model;
	GRBEnv * env;

	vector<GRBVar> yvars;		// amino acid variables
	vector< vector<GRBVar> > lvars;		// lower endpoints variables
	vector< vector<GRBVar> > uvars;		// upper endpoints variables
	vector< vector<GRBVar> > rvars;		// range variables
	vector<GRBVar> evars;				// error variables

	// results
	vector<int> yassign;
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
	int read_spectrum(const string &file);
	int compute_upper_bound();

	// add variables
	int add_distance_variables();
	int add_lower_endpoints_variables();
	int add_upper_endpoints_variables();
	int add_range_variables();
	int add_error_variables();

	// add constraints
	int add_lower_endpoints_constraints();
	int add_upper_endpoints_constraints();
	int add_range_constraints();
	int add_error_constraints();

	// set objective
	int set_objective();

	// collect results
	int collect_results();
};

#endif
