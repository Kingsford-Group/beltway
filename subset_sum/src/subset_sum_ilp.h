#ifndef __SUB_SUM_ILP_H__
#define __SUB_SUM_ILP_H__

#include "gurobi_c++.h"
#include <vector>
#include "subset_sum.h"

using namespace std;

class SubsetSumILP{
    public:
        SubsetSumILP(SubsetSum* s);
        ~SubsetSumILP();
        
        void solve();
        void print();
        void write(const string &file);
        
    private:
        int slots;
        SubsetSum* subset_sum;
        
        // members for ILP
	    GRBModel * model;
	    GRBEnv * env;

	    vector<GRBVar> tvars;		// spectrum variables
	    vector<GRBLinExpr> texpr;		// spectrum variables
	    vector<GRBVar> svars;		// edge variables
	    vector<GRBVar> cvars;		// location variables
	    
	    vector<bool> tassign;
	    vector<bool> sassign;
	    vector<bool> cassign;
	    
	    void add_spectrum_variables();
	    void add_edge_variables();
	    void add_location_variables();
	    
	    void add_spectrum_edge_constraints();
	    void add_edge_location_constraints();
	    void add_location_constraints();
	    
        // set objective
        void set_objective();

        // collect results
        void collect_results();
	    
};

#endif