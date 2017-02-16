#ifndef __ILP_H__
#define __ILP_H__

#include "gurobi_c++.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

class ilp
{
public:
	// members for ILP
	GRBModel * model;
	GRBEnv * env;

	vector<GRBVar> vars;					// gene variables

public:
	int solve();
};

#endif
