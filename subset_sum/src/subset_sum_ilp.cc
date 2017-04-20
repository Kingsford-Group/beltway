#include "subset_sum_ilp.h"

SubsetSumILP::SubsetSumILP(SubsetSum* s){
    slots = s->k;
    subset_sum = s;
    
    env = new GRBEnv();
	model = new GRBModel(*env);
	
}

void SubsetSumILP::solve(){
    try{
        add_spectrum_variables();
        add_location_variables();

        add_edge_location_constraints();
        add_spectrum_edge_constraints();
        add_location_constraints();

        set_objective();

        //model->getEnv().set(GRB_DoubleParam_TimeLimit, ilp_time_limit);

        model->update();
        
        model->write("temp.lp");

		model->optimize();
    }
	catch(GRBException &e)
	{
		printf("Error: %s\n", e.getMessage().c_str());
		exit(-1);
	}
		
}

SubsetSumILP::~SubsetSumILP(){
    svars.clear();
    tvars.clear();
    texpr.clear();
    cvars.clear();
    
    tassign.clear();
    sassign.clear();
    cassign.clear();
}

void SubsetSumILP::add_spectrum_variables(){
	tvars.clear();
	texpr.clear();
	for(int i = 0; i < subset_sum->spectrum.size(); i++){
        stringstream s;
        s << "T_" << i;
        GRBVar var = model->addVar(0, 1, 0, GRB_BINARY, s.str());
        tvars.push_back(var);
        GRBLinExpr expr;
        texpr.push_back(expr);
	}
	model->update();
}

void SubsetSumILP::add_edge_variables(){

}

void SubsetSumILP::add_location_variables(){
	cvars.clear();
	for(int i = 0; i < subset_sum->sum_array.size(); i++){
        stringstream s;
        s << "C_" << i;
        GRBVar var = model->addVar(0, 1, 0, GRB_BINARY, s.str());
        cvars.push_back(var);
	}
	
	model->addConstr(cvars[0], GRB_EQUAL, 1);
	model->addConstr(cvars[cvars.size()-1], GRB_EQUAL, 1);
	model->update();
}

void SubsetSumILP::add_edge_location_constraints(){
    svars.clear();
    for(int i = 0; i < subset_sum->sum_array.size(); i++){
        for(int j=0; j<subset_sum->sum_array[i].size(); j++){
            stringstream s_var;
            s_var << "S_" << i << "_" << subset_sum->sum_array[i][j].next_index << "_" << subset_sum->sum_array[i][j].spectrum_index;
            GRBVar var = model->addVar(0, 1, 0, GRB_BINARY, s_var.str());
            svars.push_back(var);
            
            texpr[subset_sum->sum_array[i][j].spectrum_index] += var;
            
            stringstream s_spec_edge;
            s_spec_edge << "SE_" << i << "_" << subset_sum->sum_array[i][j].next_index;
            GRBLinExpr expr = cvars[i] + cvars[subset_sum->sum_array[i][j].next_index];
            model->addConstr(2 * var, GRB_LESS_EQUAL, expr, s_spec_edge.str());
            s_spec_edge << "x";
            model->addConstr(var, GRB_LESS_EQUAL, cvars[i], s_spec_edge.str());
            s_spec_edge << "x";
            model->addConstr(var, GRB_LESS_EQUAL, cvars[subset_sum->sum_array[i][j].next_index], s_spec_edge.str());
            s_spec_edge << "x";
            model->addConstr(var, GRB_GREATER_EQUAL, expr - 1, s_spec_edge.str());
        }
    }
    model->update();
}

void SubsetSumILP::add_location_constraints(){
    GRBLinExpr expr;
    for(int i = 0; i < subset_sum->sum_array.size(); i++){
        expr += cvars[i];
    }
    model->addConstr(expr, GRB_LESS_EQUAL, (slots+1), "K_CONST");
    
}

void SubsetSumILP::add_spectrum_edge_constraints(){
    for(int i = 0; i < subset_sum->spectrum.size(); i++){
        stringstream s;
        s << "TE_" << i;
        model->addConstr(tvars[i], GRB_LESS_EQUAL, texpr[i], s.str());
    }
    model->update();
}

void SubsetSumILP::set_objective(){
    GRBLinExpr expr;
	for(int i = 0; i < tvars.size(); i++) expr += tvars.at(i);
	
	model->setObjective(expr, GRB_MAXIMIZE);
}

void SubsetSumILP::collect_results(){
    tassign.clear();
    for(int i = 0; i < tvars.size(); i++){
        tassign.push_back(tvars[i].get(GRB_DoubleAttr_X)>0.5);
    }

    sassign.clear();
    for(int i = 0; i < svars.size(); i++){
        sassign.push_back(svars[i].get(GRB_DoubleAttr_X)>0.5);
    }

    cassign.clear();
    for(int i = 0; i < cvars.size(); i++){
        cassign.push_back(cvars[i].get(GRB_DoubleAttr_X)>0.5);
    }

}

void SubsetSumILP::print(){
    try{
        collect_results();
        
        int used_count = 0;
        for(int i = 0; i < tvars.size(); i++) if(tassign[i]) used_count++;
        
        printf("Used %d of %d spectrum elements (%.2lf%%)\n",used_count,(int)tvars.size(),(100*used_count/(double)tvars.size()));
        
        int spectrum_id = 1;
        int last_found = 0;
        for(int i=1; i<cvars.size(); i++){
            if(cassign[i]){
                printf("Peptide Location %d: %d\n",spectrum_id++,(i - last_found));
                last_found = i;
            }
        }
    }
	catch(GRBException &e)
	{
		printf("Error: %s\n", e.getMessage().c_str());
		exit(-1);
	}
}