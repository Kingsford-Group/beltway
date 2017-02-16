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


public:
	MSD aa2m;
	int num_aa;
	vector<double> spectrum;

	// members for ILP
	GRBModel * model;
	GRBEnv * env;

	vector<GRBVar> vars;					// gene variables

public:
	int solve();
	int read_alphabet(const string &file);
	int read_spectrum(const string &file);
	int print();
};

#endif
