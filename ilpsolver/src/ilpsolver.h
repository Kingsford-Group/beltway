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

public:
	int solve();
	int print();

private:
	// read input files
	int read_alphabet(const string &file);
	int read_spectrum(const string &file);

	// add variables
	int add_amino_acid_variables();
	int add_lower_endpoints_variables();
	int add_upper_endpoints_variables();
	int add_range_variables();
	int add_error_variables();

};

#endif
