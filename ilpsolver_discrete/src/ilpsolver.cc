#include "ilpsolver.h"
#include "config.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <cmath> 

int less_than(vector<double> i, vector<double> j){ return (i[0]<j[0]);} 

ilpsolver::ilpsolver(const string &alphabet_file, const string &spectrum_file)
{
	read_alphabet(alphabet_file);
	read_spectrum(spectrum_file);
	compute_upper_bound();
	
	//spectrum = slice(spectrum,0,spectrum.size()-slots);

	env = new GRBEnv();
	model = new GRBModel(*env);
}

ilpsolver::~ilpsolver()
{
	if(model != NULL) delete model;
	if(env != NULL) delete env;
}

int ilpsolver::solve()
{
	try
	{
		add_amino_acid_variables();
		if(use_mvars) add_range_map_variables();
		if(!use_mvars) add_lower_endpoints_variables();
		if(!use_mvars) add_upper_endpoints_variables();
		if(no_infinity_contraints){
			add_set_map_variables();
			add_set_acid_map_variables();		
		}else{
			add_range_variables();
		}
		add_error_variables();
		
		add_amino_acid_constraints();
		if(use_mvars) add_range_map_constraints();
		if(!use_mvars) add_lower_endpoints_constraints();
		if(!use_mvars) add_upper_endpoints_constraints();
		if(no_infinity_contraints){
			add_set_map_lbound_constraints();
			add_set_map_ubound_constraints();
			add_set_map_constraints();
			add_set_acid_map_constraints();
			add_set_location_map_ubound_constraints();
			add_set_location_map_lbound_constraints();
		}else{
			add_range_constraints();
		}
		add_order_constraints();
		if(!use_mvars) add_error_constraints();
		if(use_mvars) add_error_constraints_mvars();
		//add_unique_map_constraints();
		add_anchor();
		//add_ordering_cutting_planes();
		
		//if(lp_relax) 
		add_m_linear_constraints();
		//add_s_cutting_planes();
		
		set_objective();

		model->getEnv().set(GRB_DoubleParam_TimeLimit, ilp_time_limit);

		model->update();
		
		if(lassign.size()>0 && uassign.size()>0){
		    set_mvars();
		}
        
		if(xassign.size()>0){
		    set_xvars();
		}
        
		model->write("temp.lp");

		model->optimize();

		//if(!lp_relax) 
		collect_results();
	}
	catch(GRBException &e)
	{
		printf("Error: %s\n", e.getMessage().c_str());
		exit(-1);
	}
	return 0;
}

int ilpsolver::read_alphabet(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) 
	{
		cout << "open " << file.c_str() << " error" << endl;
		return 0;
	}

	string line;
	char a[1024];
	double w;
	while(getline(fin, line))
	{
		stringstream sstr(line);
		sstr >> a >> w;
		aa_list.push_back(a);
		aa_mass.push_back(w);
		/*
		string s(a);
		if(aa2m.find(s) == aa2m.end()) aa2m.insert(PSD(s, w));
		else aa2m[s] = w;
		*/
	}

	fin.close();
	return 0;
}

int ilpsolver::read_spectrum(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) 
	{
		cout << "open " << file.c_str() << " error" << endl;
		return 0;
	}

	string line;
	char buf[10240];
	while(getline(fin, line))
	{
		if(line == "") continue;
		stringstream sstr(line);
		sstr >> buf >> buf >> slots;
		break;
	}

	spectrum.clear();
	double m;
	while(getline(fin, line))
	{
		if(line == "") continue;
		stringstream sstr(line);
		sstr >> m;
		spectrum.push_back(m);
	}

	fin.close();
	
	/*vector< vector<double> > spectrum_sort;
	spectrum_sort.clear();
	spectrum_sort.resize(spectrum.size());
	for(int p=0;p<spectrum.size();p++){
	    spectrum_sort[p].push_back(spectrum[p]);
	    spectrum_sort[p].push_back((double)p);
	}
	
	std::sort (spectrum_sort.begin(), spectrum_sort.end(), less_than);
	
	for(int sp=0;sp<spectrum.size();sp++){
	    int p = (int) spectrum_sort[sp][0];
	    printf("sorted spectrum (%d) = (%lf,%lf)\n",p,spectrum_sort[sp][0],spectrum_sort[sp][1]);
	    spectrum[sp] =p;
	}
	*/
	return 0;
}

int ilpsolver::compute_upper_bound()
{
	max_weight = 0;
	for(int i = 0; i < spectrum.size(); i++)
	{
		if(spectrum[i] > max_weight) max_weight = spectrum[i];
	}
	ubound = max_weight * slots;
	return 0;
}

int ilpsolver::add_amino_acid_variables()
{
	xvars.clear();
	xvars.resize(slots);
	for(int i = 0; i < slots; i++)
	{
		for(int k = 0; k < aa_list.size(); k++)
		{
		    stringstream s;
            s << "X_" << i << "_" << k;
			GRBVar var = model->addVar(0, 1, 0, ((lp_relax)?GRB_CONTINUOUS:GRB_BINARY), s.str());
			xvars[i].push_back(var);
		}
	}
	model->update();
	return 0;
}

int ilpsolver::add_lower_endpoints_variables()
{
	lvars.clear();
	lvars.resize(spectrum.size());
	for(int i = 0; i < spectrum.size(); i++)
	{
		for(int k = 0; k < slots; k++)
		{
		    stringstream s;
            s << "L_" << i << "_" << k;
			GRBVar var = model->addVar(0, 1, 0, ((lp_relax)?GRB_CONTINUOUS:GRB_BINARY), s.str());
			lvars[i].push_back(var);
		}
	}
	model->update();
	return 0;
}

int ilpsolver::add_upper_endpoints_variables()
{
	uvars.clear();
	uvars.resize(spectrum.size());
	for(int i = 0; i < spectrum.size(); i++)
	{
		for(int k = 0; k < slots; k++)
		{
		    stringstream s;
            s << "U_" << i << "_" << k;
			GRBVar var = model->addVar(0, 1, 0, ((lp_relax)?GRB_CONTINUOUS:GRB_BINARY), s.str());
			uvars[i].push_back(var);
		}
	}
	model->update();
	return 0;
}

int ilpsolver::add_set_map_variables(){
	ovars.clear();
	ovars.resize(slots);
	for(int i=0; i<slots; i++){
		for(int p=0; p<spectrum.size();p++){
		    stringstream s;
            s << "O_" << i << "_" << p;
			GRBVar var = model->addVar(0,1,0,((lp_relax)?GRB_CONTINUOUS:GRB_BINARY), s.str());
			ovars[i].push_back(var);
		}
	}
	model->update();
	return 0;
}

int ilpsolver::add_range_map_variables(){
    mvars.clear();
    mvars.resize(slots);
    for(int i=0;i<slots;i++){
		mvars[i].clear();
		mvars[i].resize(slots);
		for(int j=0;j<slots;j++){
			for(int p=0;p<spectrum.size();p++){
			    stringstream s;
                s << "M_" << i << "_" << j << "_" << p;
				GRBVar var = model->addVar(0,1,0,((lp_relax)?GRB_CONTINUOUS:GRB_BINARY), s.str());
				mvars[i][j].push_back(var);
			}
		}
	}
	return 0;
}

int ilpsolver::add_set_acid_map_variables(){
	svars.clear();
	svars.resize(slots);
	for(int i=0;i<slots;i++){
		svars[i].clear();
		svars[i].resize(aa_list.size());
		for(int j=0;j<aa_list.size();j++){
			for(int p=0;p<spectrum.size();p++){
			    stringstream s;
                s << "S_" << i << "_" << j << "_" << p;
				GRBVar var = model->addVar(0,1,0,((lp_relax)?GRB_CONTINUOUS:GRB_BINARY), s.str());
				svars[i][j].push_back(var);
			}
		}
	}
	return 0;
}

int ilpsolver::add_range_variables()
{
	rvars.clear();
	rvars.resize(slots);
	for(int i = 0; i < slots; i++)
	{
		for(int k = 0; k < slots; k++)
		{
		    stringstream s;
            s << "R_" << i << "_" << k;
			GRBVar var = model->addVar(0, ubound, 0, GRB_CONTINUOUS,s.str());
			rvars[i].push_back(var);
		}
	}
	model->update();
	return 0;
}

int ilpsolver::add_error_variables()
{
	evars.clear();
	for(int i = 0; i < spectrum.size(); i++)
	{
	    stringstream s;
            s << "Z_" << i;
		GRBVar var = model->addVar(0, ubound, 1, GRB_CONTINUOUS,s.str());
		evars.push_back(var);
	}
	model->update();
	return 0;
}

int ilpsolver::add_amino_acid_constraints()
{
	for(int i = 0; i < slots; i++)
	{
		GRBLinExpr expr;
		for(int k = 0; k < aa_list.size(); k++)
		{
			expr += xvars[i][k];
		}
		model->addConstr(expr, GRB_EQUAL, 1);
	}
	return 0;
}

int ilpsolver::add_lower_endpoints_constraints()
{
	for(int i = 0; i < spectrum.size(); i++)
	{
		GRBLinExpr expr;
		for(int k = 0; k < slots; k++)
		{
			expr += lvars[i][k];
		}
		model->addConstr(expr, GRB_EQUAL, 1);
	}
	return 0;
}

int ilpsolver::add_upper_endpoints_constraints()
{
	for(int i = 0; i < spectrum.size(); i++)
	{
		GRBLinExpr expr;
		for(int k = 0; k < slots; k++)
		{
			expr += uvars[i][k];
		}
		model->addConstr(expr, GRB_EQUAL, 1);
	}
	return 0;
}

int ilpsolver::add_set_map_lbound_constraints(){
	for(int k=0;k<slots;k++){
		for(int l=0;l<slots;l++){
			for(int p=0;p<spectrum.size();p++){
			    GRBLinExpr expr;
				if(!use_mvars) expr = lvars[p][k] + uvars[p][l] - 1;
				else expr = mvars[k][l][p];
				if(k<=l){
					for(int i=k;i<=l;i++){
						model->addConstr(ovars[i][p], GRB_GREATER_EQUAL, expr);
					}
				}else{
					for(int i=k;i<slots;i++){
						model->addConstr(ovars[i][p], GRB_GREATER_EQUAL, expr);
					}
					for(int i=0;i<=l;i++){
						model->addConstr(ovars[i][p], GRB_GREATER_EQUAL, expr);
					}
				}
			}
		}
	}
	return 0;
}

int ilpsolver::add_set_map_ubound_constraints(){
        for(int k=0;k<slots;k++){
                for(int l=0;l<slots;l++){
                        for(int p=0;p<spectrum.size();p++){
                                GRBLinExpr expr;
                                if(!use_mvars)  expr = 2 - lvars[p][k] - uvars[p][l];
                                else expr = 1 - mvars[k][l][p];
                                if((l+1)%slots!=k){
					if(l<k){
                                	        for(int i=l+1;i<k;i++){
                        	                        model->addConstr(ovars[i][p], GRB_LESS_EQUAL, expr);
                	                        }
        	                        }else{
	                                        for(int i=l+1;i<slots;i++){
                                        	        model->addConstr(ovars[i][p], GRB_LESS_EQUAL, expr);
                                	        }
                        	                for(int i=0;i<k;i++){
                	                                model->addConstr(ovars[i][p], GRB_LESS_EQUAL, expr);
        	                                }
	                                }
				}
                        }
                }
        }
        return 0;
}

int ilpsolver::add_set_acid_map_constraints(){
	for(int i=0;i<slots;i++){
		for(int j=0;j<aa_list.size();j++){
			for(int p=0;p<spectrum.size();p++){
				model->addConstr(svars[i][j][p], GRB_LESS_EQUAL, xvars[i][j]);
			}
		}
	}
	return 0;
}

int ilpsolver::add_set_location_map_ubound_constraints(){
        for(int i=0;i<slots;i++){
                for(int j=0;j<aa_list.size();j++){
                        for(int p=0;p<spectrum.size();p++){
                                model->addConstr(svars[i][j][p], GRB_LESS_EQUAL, ovars[i][p]);
                        }
                }
        }
        for(int i=0;i<slots;i++){
            for(int p=0;p<spectrum.size();p++){
                GRBLinExpr expr;
                for(int j=0;j<aa_list.size();j++){
                    expr += svars[i][j][p];
                }
                model->addConstr(expr, GRB_EQUAL, ovars[i][p]);
            }
        }
        
        
        return 0;
}

int ilpsolver::add_set_location_map_lbound_constraints(){
        for(int i=0;i<slots;i++){
        	for(int p=0;p<spectrum.size();p++){
			GRBLinExpr expr;
			for(int j=0;j<aa_list.size();j++){
				expr += svars[i][j][p];
                        }
			model->addConstr(expr, GRB_GREATER_EQUAL, ovars[i][p]);
                }
        }
        return 0;
}

int ilpsolver::add_set_map_constraints(){
	for(int p=0;p<spectrum.size();p++){
		GRBLinExpr expr;
		for(int i=0;i<slots;i++){
			for(int j=0;j<aa_list.size();j++){
				expr += svars[i][j][p];
			}
		}
		model->addConstr(expr, GRB_GREATER_EQUAL, 1);
	}
	return 0;
}

int ilpsolver::add_s_cutting_planes(){
	for(int p=0;p<spectrum.size();p++){
		for(int i=0;i<slots;i++){
			for(int j=0;j<aa_list.size();j++){
			    GRBLinExpr expr = xvars[i][j] + ovars[i][p] - 1;
				model->addConstr(svars[i][j][p], GRB_GREATER_EQUAL, expr);
			}
		}
		
	}
	return 0;
}

int ilpsolver::add_range_map_constraints(){
    for(int p=0;p<spectrum.size();p++){
        GRBLinExpr expr;
        for(int k = 0; k < slots; k++){
            for(int l = 0; l < slots; l++){
                expr += mvars[k][l][p];
            }
        }
        model->addConstr(expr, GRB_EQUAL, 1);
    }
}


int ilpsolver::add_m_linear_constraints(){

    for(int p=0;p<spectrum.size();p++){
        for(int i=0;i<slots;i++){
            GRBLinExpr expr;
            for(int kp=0;kp<slots;kp++){
                for(int lp=0;lp<slots;lp++){
                    int k = (i - kp + slots) % slots;
                    int l = (i + lp) % slots;
                    expr += mvars[k][l][p];
                }
            }
            model->addConstr(ovars[i][p], GRB_LESS_EQUAL, expr);
            
            GRBLinExpr expr2;
            for(int kp=1;kp<slots-1;kp++){
                int k = (i + kp) % slots;
                for(int lp=0;lp<(slots-kp-1);lp++){
                    int l = (k + lp) % slots;
                    expr += mvars[k][l][p];
                }
            }
            model->addConstr(ovars[i][p], GRB_GREATER_EQUAL, expr2);
        }
    }

    /*for(int k=0;k<slots;k++){
        for(int l=0;l<slots;l++){
            for(int p=0;p<spectrum.size();p++){
                GRBLinExpr expr = 1 - mvars[k][l][p];
                if((l+1)%slots!=k){
                    if(l<k){
                        for(int i=l+1;i<k;i++){
                            model->addConstr(ovars[i][p], GRB_LESS_EQUAL, expr);
                        }
                    }else{
                        for(int i=l+1;i<slots;i++){
                            model->addConstr(ovars[i][p], GRB_LESS_EQUAL, expr);
                        }
                        for(int i=0;i<k;i++){
                            model->addConstr(ovars[i][p], GRB_LESS_EQUAL, expr);
                        }
                    }
                }
            }
        }
    }*/
        
    return 0;
}

int ilpsolver::add_range_constraints()
{
	for(int k = 0; k < slots; k++)
	{
		for(int l = 0; l < slots; l++)
		{
			GRBLinExpr expr;
			if(k <= l)
			{
				for(int i = k; i <= l; i++)
				{
					for(int j = 0; j < aa_list.size(); j++)
					{
						expr += xvars[i][j] * aa_mass[j];
					}
				}
			}
			else
			{
				for(int i = k; i < slots; i++)
				{
					for(int j = 0; j < aa_list.size(); j++)
					{
						expr += xvars[i][j] * aa_mass[j];
					}
				}
				for(int i = 0; i <= l; i++)
				{
					for(int j = 0; j < aa_list.size(); j++)
					{
						expr += xvars[i][j] * aa_mass[j];
					}
				}
			}
			model->addConstr(rvars[k][l], GRB_EQUAL, expr);
		}
	}
	return 0;
}

int ilpsolver::add_error_constraints()
{
	if(no_infinity_contraints){
		for(int p=0;p<spectrum.size();p++){
			GRBLinExpr expr1 = spectrum[p];
                        GRBLinExpr expr2;
			for(int i=0;i<slots;i++){
				for(int j=0;j<aa_list.size();j++){
					expr1 -= svars[i][j][p] * aa_mass[j];
					expr2 += svars[i][j][p] * aa_mass[j];
				}
			}
                	expr2 -= spectrum[p];
			model->addConstr(evars[p], GRB_GREATER_EQUAL, expr1);
                        model->addConstr(evars[p], GRB_GREATER_EQUAL, expr2);
		}
	}else{
		for(int p = 0; p < spectrum.size(); p++)
		{
			for(int k = 0; k < slots; k++)
			{
				for(int l = 0; l < slots; l++)
                        	{
					double bound = 0;
					if(k <= l) bound = (l - k + 1) * max_weight;
					else bound = (l + 1 + slots - k) * max_weight;
					GRBLinExpr expr1 = rvars[k][l] - spectrum[p] + bound * (lvars[p][k] + uvars[p][l] - 2);
					GRBLinExpr expr2 = spectrum[p] - rvars[k][l] + bound * (lvars[p][k] + uvars[p][l] - 2);
					model->addConstr(evars[p], GRB_GREATER_EQUAL, expr1);
					model->addConstr(evars[p], GRB_GREATER_EQUAL, expr2);
				}
			}
		}
	}
	return 0;
}


int ilpsolver::add_error_constraints_mvars()
{
	if(no_infinity_contraints){
		for(int p=0;p<spectrum.size();p++){
			GRBLinExpr expr1 = spectrum[p];
                        GRBLinExpr expr2;
			for(int i=0;i<slots;i++){
				for(int j=0;j<aa_list.size();j++){
					expr1 -= svars[i][j][p] * aa_mass[j];
					expr2 += svars[i][j][p] * aa_mass[j];
				}
			}
                	expr2 -= spectrum[p];
			model->addConstr(evars[p], GRB_GREATER_EQUAL, expr1);
                        model->addConstr(evars[p], GRB_GREATER_EQUAL, expr2);
		}
	}else{
		for(int p = 0; p < spectrum.size(); p++)
		{
			for(int k = 0; k < slots; k++)
			{
				for(int l = 0; l < slots; l++)
                        	{
					double bound = 0;
					if(k <= l) bound = (l - k + 1) * max_weight;
					else bound = (l + 1 + slots - k) * max_weight;
					GRBLinExpr expr1 = rvars[k][l] - spectrum[p] + bound * mvars[k][l][p];
					GRBLinExpr expr2 = spectrum[p] - rvars[k][l] + bound * mvars[k][l][p];
					model->addConstr(evars[p], GRB_GREATER_EQUAL, expr1);
					model->addConstr(evars[p], GRB_GREATER_EQUAL, expr2);
				}
			}
		}
	}
	return 0;
}

int ilpsolver::add_order_constraints(){
    
    for(int j=0; j<aa_list.size(); j++){
        for(int i=0;i<slots;i++){
            for(int k=0; k<j;k++){
                stringstream s;
                s << "OC_" << j << "_" << i << "_" << k;
                model->addConstr(1 - xvars[0][j], GRB_GREATER_EQUAL, xvars[i][k],s.str());
            }
        }
    }
}

int ilpsolver::add_ordering_cutting_planes()
{
    if(no_infinity_contraints){
		for(int pi1=0;pi1<spectrum.size();pi1++){
		    for(int pi2=pi1+1;pi2<spectrum.size();pi2++){
                if(spectrum[pi1] != spectrum[pi2]){
                    int p1 = pi1;
                    int p2 = pi2;
                    if(spectrum[p1] < spectrum[p2]){
                        p1 = pi2;
                        p2 = pi1;
                    }
                    GRBLinExpr expr;
                    for(int i=0;i<slots;i++){
                        for(int j=0;j<aa_list.size();j++){
                            expr += svars[i][j][p1] * aa_mass[j];
					        expr -= svars[i][j][p2] * aa_mass[j];
                        }
                    }
                    stringstream s;
                    s << "CPR_" << pi1 << "_" << pi2 ;
                    model->addConstr(expr, GRB_GREATER_EQUAL, 0, s.str());
                }
            }
        }
    }else{
        for(int p1 = 0; p1 < spectrum.size(); p1++)
        {
            for(int p2 = p1 + 1; p2 < spectrum.size(); p2++)
            {
                for(int k1 = 0; k1 < slots; k1++)
                {
                    for(int l1 = 0; l1 < slots; l1++)
                    {
                        for(int k2 = 0; k2 < slots; k2++)
                        {
                            for(int l2 = 0; l2 < slots; l2++)
                            {
                                if(spectrum[p1] > spectrum[p2])
                                {
                                    double bound = 0;
                                    if(k2 <= l2) bound = (l2 - k2 + 1) * max_weight;
                                    else bound = (l2 + 1 + slots - k2) * max_weight;
                                    GRBLinExpr expr1 = rvars[k1][l1] - rvars[k2][l2];   
                                    GRBLinExpr expr2;
                                    if(!use_mvars) expr2 = (lvars[p1][k1] + uvars[p1][l1] + lvars[p2][k2] + uvars[p2][l2] - 4) * bound;
                                    else expr2 = (mvars[k1][l1][p1] + mvars[k2][l2][p2] - 2) * bound;
                                    model->addConstr(expr1, GRB_GREATER_EQUAL, expr2);
                                }
                                else if(spectrum[p2] > spectrum[p1])
                                {
                                    double bound = 0;
                                    if(k1 <= l1) bound = (l1 - k1 + 1) * max_weight;
                                    else bound = (l1 + 1 + slots - k1) * max_weight;
                                    GRBLinExpr expr1 = rvars[k2][l2] - rvars[k1][l1];
                                    GRBLinExpr expr2;
                                    if(!use_mvars) expr2 = (lvars[p1][k1] + uvars[p1][l1] + lvars[p2][k2] + uvars[p2][l2] - 4) * bound;
                                    else expr2 = (mvars[k1][l1][p1] + mvars[k2][l2][p2] - 2) * bound;
                                    model->addConstr(expr1, GRB_GREATER_EQUAL, expr2);
                                }
                            }
                        }
                    }
                }
            }
		}
	}
	return 0;
}


int ilpsolver::add_anchor()
{
	GRBLinExpr expr;
	for(int l = 0; l < slots; l++){
	    expr += mvars[0][l][0];
	}
	model->addConstr(expr, GRB_EQUAL, 1);
	
	return 0;
}

int ilpsolver::add_unique_map_constraints(){
        for(int p1=0;p1<spectrum.size();p1++){
		    for(int p2=p1+1;p2<spectrum.size();p2++){
		        for(int k=0;k<slots;k++){
		            for(int l=0;l<slots;l++){
		                GRBLinExpr expr = mvars[k][l][p1] + mvars[k][l][p2];
		                model->addConstr(expr, GRB_LESS_EQUAL, 1);
		            }
		        }
		    }
		}
        return 0;
}

int ilpsolver::set_objective()
{
	GRBQuadExpr expr;
	for(int i = 0; i < evars.size(); i++) expr += evars.at(i);
	
	if(lp_relax){
	    for(int p = 0; p < spectrum.size(); p++){
			for(int j = 0; j < aa_list.size(); j++){
				for(int i = 0; i < slots; i++){
				    //expr -=  svars[i][j][p]*svars[i][j][p] - svars[i][j][p];
				}
			}
		}
	}
	
	model->setObjective(expr, GRB_MINIMIZE);
	return 0;
}

int ilpsolver::collect_results()
{
	xassign.clear();
	for(int i = 0; i < slots; i++)
	{
		int k = -1;
		double k_v = -1;
		for(int j = 0; j < aa_list.size(); j++)
		{
			if(lp_relax){
                if(xvars[i][j].get(GRB_DoubleAttr_X) > k_v || k==-1){
                    k = j;
                    k_v = xvars[i][j].get(GRB_DoubleAttr_X);
                }
			}else{
                if(xvars[i][j].get(GRB_DoubleAttr_X) <= 0.5) continue;
                assert(k == -1);
                k = j;
			}
		}
		assert(k >= 0);
		xassign.push_back(k);
	}

	if(no_infinity_contraints){
		sjassign.clear();
		sjassign.resize(spectrum.size());
		siassign.clear();
		siassign.resize(spectrum.size());
		for(int p = 0; p < spectrum.size(); p++){
			for(int j = 0; j < aa_list.size(); j++){
				for(int i = 0; i < slots; i++){
					if(svars[i][j][p].get(GRB_DoubleAttr_X) <= 0.5) continue;
					sjassign[p].push_back(j);
					siassign[p].push_back(i);
				}
			}
		}
	}

	lassign.clear();
	massign.clear();
    uassign.clear();
    
	if(!use_mvars){
        for(int p = 0; p < spectrum.size(); p++)
        {
            int k = -1;
            for(int j = 0; j < slots; j++)
            {
                if(lvars[p][j].get(GRB_DoubleAttr_X) <= 0.5) continue;
                assert(k == -1);
                k = j;
            }
            assert(k >= 0);
            lassign.push_back(k);
        }
        
        for(int p = 0; p < spectrum.size(); p++)
        {
            int k = -1;
            for(int j = 0; j < slots; j++)
            {
                if(uvars[p][j].get(GRB_DoubleAttr_X) <= 0.5) continue;
                assert(k == -1);
                k = j;
            }
            assert(k >= 0);
            uassign.push_back(k);
        }
    }else{
	    for(int p = 0; p < spectrum.size(); p++){
            int ks = -1;
            int ls = -1;
            double vs = -1;
            for(int k = 0; k < slots; k++){
                for(int l = 0; l < slots; l++){
                    if(lp_relax){
                        if(mvars[k][l][p].get(GRB_DoubleAttr_X) > vs){
                            vs = mvars[k][l][p].get(GRB_DoubleAttr_X);
                            ks = k;
                            ls = l;
                        }
                    }else{
                        if(mvars[k][l][p].get(GRB_DoubleAttr_X) <= 0.5) continue;
                        assert(ks == -1);
                        assert(ls == -1);
                        ks = k;
                        ls = l;
                    }
                }
            }
            lassign.push_back(ks);
            uassign.push_back(ls);
        }
    }
    
	if(!no_infinity_contraints){
		wassign.clear();
		for(int p = 0; p < spectrum.size(); p++)
		{
			int l = lassign[p];
			int u = uassign[p];
			assert(l >= 0 && l < slots);
			assert(u >= 0 && u < slots);
			double w = rvars[l][u].get(GRB_DoubleAttr_X);
			wassign.push_back(w);
		}
	}

	eassign.clear();
	for(int p = 0; p < spectrum.size(); p++)
	{
		double e = evars[p].get(GRB_DoubleAttr_X);
		eassign.push_back(e);
	}

	return 0;
}

int ilpsolver::print_relaxed(){

    vector<double> slot_assign;
    slot_assign.clear();
	for(int i = 0; i < slots; i++)
	{
	    double sum = 0.0;
		for(int j = 0; j < aa_list.size(); j++)
		{
			double k = xvars[i][j].get(GRB_DoubleAttr_X);
			if(k > 0){
			    sum += (k * aa_mass[j]);
			    //printf("X[%d][%d] = %.4lf\n",i,j,k);
			    printf("slot %d is assigned amino acid %d (%.2lf%%): %s -> %.3lf\n", i, j, k, aa_list[j].c_str(), aa_mass[j]);
			}
		}
		printf("slot %d assigned total weight of %.3lf\n", i, sum);
		slot_assign.push_back(sum);
	}

    for(int i = 0; i < spectrum.size(); i++){
        if(no_infinity_contraints){
            for(int k = 0; k < slots; k++){
                for(int l = 0; l < slots; l++){
                    //printf("(k,l,i) = (%d,%d,%d)\n",k,l,i);
                    double v = mvars[k][l][i].get(GRB_DoubleAttr_X);
                    double w = 0;
                    if(v>0){  
                        if(k<=l){
                            for(int m=k;m<=l;m++) w += slot_assign[m];
                        }else{
                            for(int m=k;m<slots;m++) w += slot_assign[m];
                            for(int m=0;m<=l;m++) w += slot_assign[m];
                        }
                        double e = evars[i].get(GRB_DoubleAttr_X);
                        printf("spectrum %d with mass %.3lf is assigned to interval [%d, %d], with mass %.3lf and weight %.3lf (actual %.3lf)\n", i, spectrum[i], k, l, w, v, (v*abs(w-spectrum[i])));
                    }
                }
            }
        }
        
    }

    return 0;
}

int ilpsolver::print()
{
	
	//if(lp_relax) return print_relaxed();
	/*
	for(MSD::iterator it = aa2m.begin(); it != aa2m.end(); it++)
	{
		string s = it->first;
		double m = it->second;
		printf("amino acid %s -> %.3lf\n", s.c_str(), m);
	}
	*/
	assert(aa_list.size() == aa_mass.size());
	for(int i = 0; i < aa_list.size(); i++)
	{
		printf("amino acid %d: %s -> %.3lf\n", i, aa_list[i].c_str(), aa_mass[i]);
	}
	printf("number of amino acid = %d\n", slots);
	for(int i = 0; i < spectrum.size(); i++)
	{
		printf("spectrum: %.3lf\n", spectrum[i]);
	}
	printf("max weight = %.3lf, upper bound = %.3lf\n", max_weight, ubound);

	for(int i = 0; i < xassign.size(); i++)
	{
		int k = xassign[i];
		printf("slot %d is assigned amino acid %d: %s -> %.3lf\n", i, k, aa_list[k].c_str(), aa_mass[k]);
	}
	for(int i = 0; i < spectrum.size(); i++)
	{
		int l = lassign[i];
		int u = uassign[i];
		double w = 0;
		if(no_infinity_contraints){
			if(l<=u){
				for(int i=l;i<=u;i++) w += aa_mass[xassign[i]];
			}else{
				for(int i=l;i<slots;i++) w += aa_mass[xassign[i]];
				for(int i=0;i<=u;i++) w += aa_mass[xassign[i]];
			}
		}else w = wassign[i];
		double e = -10;
        if(eassign.size()>i) e = eassign[i];
		printf("spectrum %d with mass %.3lf is assigned to interval [%d, %d], with mass %.3lf and error %.3lf (actual %.3lf)\n", i, spectrum[i], l, u, w, e, (w-spectrum[i]));
		
		if(false && no_infinity_contraints){
			printf("\t");
			double sum = 0;
			for(int c = 0; c<sjassign[i].size(); c++){
				int j = sjassign[i][c];
				if(c!=0) printf(",");
				printf("%d (%s)",j, aa_list[j].c_str());
				sum += aa_mass[j];
			}
			printf(" = %.3lf\n",sum);

			sum = 0;
			printf("\t");
                        for(int c = 0; c<siassign[i].size(); c++){
                                int j = siassign[i][c];
                                if(c!=0) printf(",");
                                printf("%d (%s)",j, aa_list[xassign[j]].c_str());
                                sum += aa_mass[xassign[j]];
                        }
                        printf(" = %.3lf\n",sum);
		}
	}
	return 0;
}

int ilpsolver::write(const string &file)
{
	ofstream fout(file.c_str());
	if(fout.fail()) return 0;

	for(int i = 0; i < xassign.size(); i++)
	{
		int k = xassign[i];
		fout << aa_list[k].c_str() << endl;
	}

	fout.close();
	return 0;
}

int ilpsolver::set_mvars(){
    for(int p = 0; p < uassign.size(); p++){
	    int ks = -1;
	    int ls = -1;
	    for(int k = 0; k < slots; k++){
			for(int l = 0; l < slots; l++){
			    //printf("Setting M_%d_%d_%d to 0\n",k,l,p);
                mvars[k][l][p].set(GRB_DoubleAttr_Start,0);
            }
        }
        
	    printf("Setting M_%d_%d_%d to 1\n",lassign[p],uassign[p],p);
        mvars[lassign[p]][uassign[p]][p].set(GRB_DoubleAttr_Start,1.0);
	}
	return 0;
}

int ilpsolver::set_xvars(){
	for(int i = 0; i < slots; i++){
		for(int j = 0; j < aa_list.size(); j++){
		    if(xassign[i] == j){
		        printf("Setting X_%d_%d to 1\n",i,j);
		        xvars[i][j].set(GRB_DoubleAttr_Start, 1.0);
		    }else{
		        //printf("Setting X_%d_%d to 0\n",i,j);
		        xvars[i][j].set(GRB_DoubleAttr_Start, 0);
		    }
		}
	}
	return 0;
}

int ilpsolver::reset(){
    env = new GRBEnv();
	model = new GRBModel(*env);
}



int ilpsolver::greedy_warm_start(){
    
    xassign.clear();
	for(int i = 0; i < slots; i++)
	{
        xassign.push_back(-1);
	}
	
	lassign.clear();
	uassign.clear();
	
	double last_assign = 0;
	
	vector< vector<double> > spectrum_sort;
	spectrum_sort.clear();
	spectrum_sort.resize(spectrum.size());
	for(int p=0;p<spectrum.size();p++){
	    spectrum_sort[p].push_back(spectrum[p]);
	    spectrum_sort[p].push_back((double)p);
	
	    lassign.push_back(-1);
	    uassign.push_back(-1);
	}
	printf("In Greedy\n");
	
	vector< vector<double> > ranges;
	std::sort (spectrum_sort.begin(), spectrum_sort.end(), less_than);
	int num_assigned = 0;
	ranges.clear();
	for(int sp=0;sp<spectrum.size();sp++){
	    int p = (int) spectrum_sort[sp][1];
	    printf("sorted spectrum (%d) = (%lf,%lf)\n",p,spectrum_sort[sp][0],spectrum_sort[sp][1]);
	    
	    int min_range_j = -1;
	    double min_range_dist = -1;
	    if(ranges.size()>0){
            min_range_dist = abs((double)(ranges[0][0] - spectrum_sort[sp][0]));
            min_range_j = 0;
            int min_j = 0;
            for(int j=0;j<ranges.size();j++){
                double dist = abs((ranges[j][0] - spectrum_sort[sp][0]));
                printf("j: %d\tDist: %.0lf\tMin Dist: %.0lf(%d)\n",j,dist,min_range_dist, min_range_j); 
                if(min_range_dist > dist){
                    min_range_j = j;
                    min_range_dist = dist;
                }
            }
	    }
	    
	    
        double min_slot_dist = abs((double)(aa_mass[0] - spectrum_sort[sp][0]));
        int min_slot_j = 0;
        
	    //printf("Num Assigned: %d\tSlots: %d\n",num_assigned,slots); 
	    if(num_assigned < slots){
	    
	        //printf("Num Assigned: %d\tSlots: %d\n",num_assigned,slots); 
	        for(int j=0;j<aa_mass.size();j++){
	            double dist = abs((aa_mass[j] - spectrum_sort[sp][0]));
	            
	            //printf("j: %d\tDist: %.0lf\tMin Dist: %.0lf(%d)\n",j,dist,min_slot_dist, min_slot_j); 
	            if(min_slot_dist > dist){
	                min_slot_j = j;
	                min_slot_dist = dist;
	            }
	        }
	      
	        
	    }else{
	        min_slot_dist = -1;
	        min_slot_j = -1;
	    }
	    
	    printf("ranges (%d,%.0lf)\tslots (%d,%.0lf)\n",min_range_j,min_range_dist,min_slot_j,min_slot_dist);
	    if((min_range_dist == -1 || min_slot_dist < min_range_dist) && min_slot_dist != -1){
	          for(int i = 0; i < slots && lassign[p] == -1; i++){
	            if(xassign[i] == min_slot_j){
	                lassign[p] = i;
	                uassign[p] = i;
	                if(last_assign <= aa_mass[xassign[i]]){
	                    last_assign = aa_mass[xassign[i]];
	                }
	                else{
	                    printf("error 1\n");
	                    exit(-1);
	                }
	            }
	        }
	        if(uassign[p] == -1){
	            xassign[num_assigned] = min_slot_j;
	            lassign[p] = num_assigned;
                uassign[p] = num_assigned;
                num_assigned++;
                
                if(last_assign <= aa_mass[min_slot_j]){
                    last_assign = aa_mass[min_slot_j];
                }
                else{
                    printf("error 2\n");
                    exit(-1);
                }
                
                double sum = 0;
                for(int k=num_assigned-1;k>=0;k--){
                    sum += aa_mass[xassign[k]];
                
                    vector<double> temp;
                    temp.push_back(sum);
                    temp.push_back(k);
                    temp.push_back((num_assigned-1));
                    ranges.push_back(temp);
                    printf("Adding range %.0lf (%d,%d)\n",sum,k,(num_assigned-1));
	            }
	        }
	        
	        
	        if(num_assigned == slots){
	            double total = 0;
	            for(int i = 0; i < slots; i++){
	                total += aa_mass[xassign[i]];
	            }
	            vector<double> temp;
	            temp.push_back(total);
	            temp.push_back(0);
	            temp.push_back(slots-1);
	            ranges.push_back(temp);
	            
	            for(int k=0;k<slots;k++){
	                for(int l=k;l<slots;l++){
	                    double sum = 0;
	                    for(int i=k;i<=l;i++){
	                        sum += aa_mass[xassign[i]];
	                    }
	                    
                        vector<double> temp;
                        temp.push_back(sum);
                        temp.push_back(k);
                        temp.push_back(l);
	                    ranges.push_back(temp);
	                    printf("Adding range %.0lf (%d,%d)\n",sum,k,l);
	                    vector<double> temp2;
                        temp2.push_back(total-sum+aa_mass[xassign[k]]+aa_mass[xassign[l]]);
                        temp2.push_back(l);
                        temp2.push_back(k);
	                    if(k!=l){ranges.push_back(temp2);}
	                }
	            }
	        }
	    }else{	        
            lassign[p] = ranges[min_range_j][1];
            uassign[p] = ranges[min_range_j][2];
            
            if(last_assign <= ranges[min_range_j][0]){
                last_assign = ranges[min_range_j][0];
            }
            else{
                printf("Assigned range (%d) %.0lf, last assigned %.0lf\n",min_range_j,ranges[min_range_j][0],last_assign);
                exit(-1);
            }
	                
	    }
	    
	}
	
	int min_x_index = aa_mass.size();
	int min_x_slot = -1;
	for(int i=0;i<xassign.size();i++){
	    if(xassign[i] < min_x_index){
	        min_x_index = xassign[i];
	        min_x_slot = i;
	    }
	}
	
	vector<int> xassign_save ;
	xassign_save = xassign;
	for(int i=0;i<xassign.size();i++){
	    printf("Shifting %d->%d (%d)\n",xassign[i],xassign_save[(i+min_x_slot)%slots],i);
	    xassign[i] = xassign_save[(i+min_x_slot)%slots];
	}
	for(int p=0;p<spectrum.size();p++){
	    lassign[p] = (lassign[p]+min_x_slot)%slots;
        uassign[p] = (uassign[p]+min_x_slot)%slots;
	}
	uassign.clear();
	lassign.clear();
	return 0;
}

double ilpsolver::graph_greedy_warm_start_helper(int edges, bool complete, int sp, vector< vector<double> > spectrum_sort, vector<int> assigned, vector<int> connectable, vector< vector <double> > connections, vector< vector<int> > direct_edges, vector<int> assignment, double total_error, vector<int> my_xassign, vector<int> my_lassign, vector<int> my_uassign){
    
    double equal_precision = 0.0001;
    
    if((best_total_error != -1 && best_total_error < (total_error - equal_precision)) || best_total_error == 0){
        return best_total_error;
    }

    if(sp == spectrum_sort.size()){
        printf("Total Error: %.0lf\tBest Error: %.0lf\n",total_error,best_total_error);
        if(total_error < (best_total_error - equal_precision) || best_total_error == -1){
            best_total_error = total_error;
            xassign = my_xassign;
            lassign = my_lassign;
            uassign = my_uassign;
        }
        return best_total_error;   
    }
    
    int p = (int) spectrum_sort[sp][1];
    
    
     if(!complete && edges==slots-1 && assigned.size() == slots){
        complete = true;
        
        //printf("Completed with assigned.size() = %d\tuassign.size() = %d\tlassign.size() = %d\n",assigned.size(),uassign.size(),lassign.size());
        int min_x_index = aa_mass.size();
        int min_x_slot = -1;
        for(int i=0;i<xassign.size();i++){
            //printf("Xassign[%d] = %d\n",i,my_xassign[i]);
            if(my_xassign[i] < min_x_index){
                min_x_index = my_xassign[i];
                min_x_slot = i;
            }
        }
        
        for(int i=0; i<assigned.size(); i++){
            for(int j=0; j<assigned.size(); j++){
                if(min_x_slot == i && min_x_slot == j) printf("*");
                printf("%.0lf (%d)\t",connections[i][j],direct_edges[i][j]);
            }
            printf("\n");
        }

        printf("Connectable: ");
        for(int i=0; i<connectable.size(); i++){
            printf("%d\t",connectable[i]);
        }
        printf("\n");
        
        vector<int> new_assignment;
        int i=min_x_slot;
        int i_prev = -1;
        new_assignment.push_back(i);
        while(new_assignment.size()<assignment.size()){
            bool found_new_j = false;
            for(int j=0;j<slots && !found_new_j;j++){
                //printf("i: %d\tj:%d\ti_prev: %d\tDE: %d\n",i,j,i_prev,direct_edges[i][j]);
                if(direct_edges[i][j] && j!=i_prev){
                    new_assignment.push_back(j);
                    i_prev = i;
                    i = j;
                    found_new_j = true;
                    break;
                }
                if(j!=i_prev && i!=j && connectable[i] == 1 && connectable[j]==1){
                    new_assignment.push_back(j);
                    i_prev = i;
                    i = j;
                    found_new_j = true;
                    break;
                }
            }
            assert(found_new_j);
            //printf("New Assignment %d -> %d\n",new_assignment[new_assignment.size()-1],new_assignment.size()-1);
        }
        
        for(int i=0;i<my_uassign.size();i++){
            if(my_uassign[i] != -1){
                int k = -1;
                for(int j=0;j<new_assignment.size();j++){
                    if(my_uassign[i] == new_assignment[j]) k = j;
                }
                printf("Shifting uassign[%d]: %d->%d\n",i,my_uassign[i],k);
                assert(k!=-1);
                my_uassign[i] = k;
            }
        }
        for(int i=0;i<my_lassign.size();i++){
        
            if(my_lassign[i] != -1){
                int k = -1;
                for(int j=0;j<new_assignment.size();j++){
                    if(my_lassign[i] == new_assignment[j]) k = j;
                }
                printf("Shifting lassign[%d]: %d->%d\n",i,my_uassign[i],k);
                assert(k!=-1);
                my_lassign[i] = k;
            }
        }
        
        string my_cycle = "";
        vector<int> save_xassign = my_xassign;
        for(int i=0;i<my_xassign.size();i++){
            //printf("New Xassign: %d\tOld: %d\t(%d->%d)\n",my_xassign[i],save_xassign[new_assignment[i]],i,new_assignment[i]);
            my_xassign[i] = save_xassign[new_assignment[i]];
            my_cycle += aa_list[my_xassign[i]];
            connections[i][i] = aa_mass[my_xassign[i]];
        }
        
        bool found_cycle = false;
        for(int i=0;i<greedy_search_cycles.size();i++){
            if(my_cycle == greedy_search_cycles[i]){
                found_cycle = true;
                if(greedy_search_cycles_error[i] > total_error){
                    greedy_search_cycles_error[i] = total_error;
                }else{
                    //return total_error;
                }
            }
        }
        if(!found_cycle){
            greedy_search_cycles.push_back(my_cycle);
            greedy_search_cycles_error.push_back(total_error);
        }
        
        printf("My Cycle: %s\tMy Error: %.0lf\n",my_cycle.c_str(),total_error);
        
        //do reverse range calculations
        double sum = 0;
        for(int i=0; i<assigned.size(); i++){
            sum += connections[i][i];
        }
        
        for(int i=0;i<slots;i++){
            double sum = 0;
            for(int k=0;k<slots;k++){
                int j = (i+k) % slots;
                sum += aa_mass[my_xassign[j]];
                connections[i][j] = sum;
                
            //for(int j=i;j<slots;j++){
                /*int in_sum = 0;
                for(int k=i;k<=j;k++) in_sum += aa_mass[my_xassign[k]];
                connections[i][j] = in_sum;
                if(i!=j) connections[j][i] = sum - in_sum + aa_mass[my_xassign[i]] + aa_mass[my_xassign[j]];*/
                
            }
        }
        printf("After Recalc:\n");
        for(int i=0; i<slots; i++){
            for(int j=0; j<slots; j++){
                printf("%.0lf (%d)\t",connections[i][j],direct_edges[i][j]);
            }
            printf("\n");
        }
        
        
        for(int i=0;i<my_lassign.size();i++){
            if(my_lassign[i] != -1 && my_uassign[i] != -1){
                if(abs(connections[my_lassign[i]][my_lassign[i]] - spectrum[i]) > abs(connections[my_uassign[i]][my_lassign[i]] - spectrum[i])){
                    int temp = my_lassign[i];
                    my_lassign[i] = my_uassign[i];
                    my_uassign[i] = temp;
                }
                printf("Shifted assign[%d]: %d,%d\n",i,my_lassign[i],my_uassign[i]);
            }else if(my_lassign[i] != -1 || my_uassign[i] != -1){
                printf("Assigned only 1!! assign[%d]: %d,%d\n",i,my_lassign[i],my_uassign[i]); 
                exit(20);
            }
        }
        
    }
    
    
    
    double min_slot_dist = abs((double)(aa_mass[0] - spectrum_sort[sp][0]));
    vector<int> min_slot_j;
    min_slot_j.clear();
    
    if(assigned.size() < slots){
    
        for(int j=0;j<aa_mass.size();j++){
            double dist = abs((aa_mass[j] - spectrum_sort[sp][0]));
            
            if(min_slot_dist >= (dist - equal_precision)){
                if(min_slot_dist > (dist + equal_precision)){
                    min_slot_j.clear();
                    min_slot_dist = dist;
                }
                min_slot_j.push_back(j);
            }
        }
    }else{
        min_slot_dist = -1;
        min_slot_j.clear();
    }
    
    vector<int> min_range_i;
    min_range_i.clear();
    vector<int> min_range_j;
    min_range_j.clear();
    double min_range_dist = -1;
    
    vector<int> min_new_edge_i;
    vector<int> min_new_edge_j;
    vector<int> min_new_edge_k;
    vector<int> min_new_edge_l;
    min_new_edge_i.clear();
    min_new_edge_j.clear();
    min_new_edge_k.clear();
    min_new_edge_l.clear();
    double min_new_edge_dist = -1;
    
    vector<int> new_aa_edge_a;
    vector<int> new_aa_edge_i;
    vector<int> new_aa_edge_j;
    new_aa_edge_a.clear();
    new_aa_edge_i.clear();
    new_aa_edge_j.clear();
    double new_aa_edge_dist = -1;
    
    if(assigned.size()>0){
        for(int i = 0; i < assigned.size(); i++){
            for(int j=0;j<assigned.size();j++){
                if(connections[i][j]!=-1){
                    double dist = abs(connections[i][j] - spectrum_sort[sp][0]);
                    if(min_range_dist >= (dist - equal_precision) || min_range_dist == -1){
                        if(min_range_dist > (dist + equal_precision) || min_range_dist == -1){
                            min_range_dist = dist;
                            min_range_i.clear();
                            min_range_j.clear();
                        }
                        min_range_i.push_back(i);
                        min_range_j.push_back(j);
                    }
                }
            }
        }
        
        if(!complete){
            for(int i=0; i<assigned.size(); i++){
                for(int j=0;j<assigned.size();j++){
                    for(int k=0;k<assigned.size();k++){
                        for(int l=0;l<assigned.size();l++){
                            //printf("i,k: %.0lf\tl,j: %.0lf\tk: %d\tl: %d\tk,l: %.0lf\ti,j: %.0lf\n",connections[i][k],connections[l][j],connectable[k],connectable[l],connections[k][l],connections[i][j]);
                            if(connections[i][k] != -1 && connections[l][j] != -1 && connectable[k]>0 && connectable[l]>0 && connections[k][l] == -1 && connections[i][j]==-1){
                                double dist = abs(connections[i][k] + connections[l][j] - spectrum_sort[sp][0]);
                                if(min_new_edge_dist >= (dist - equal_precision) || min_new_edge_dist == -1){
                                    if(min_new_edge_dist > (dist + equal_precision) || min_new_edge_dist == -1){
                                        min_new_edge_dist = dist;
                                            min_new_edge_i.clear();
                                            min_new_edge_j.clear();
                                            min_new_edge_k.clear();
                                            min_new_edge_l.clear();
                                    }
                                    min_new_edge_i.push_back(i);
                                    min_new_edge_j.push_back(j);
                                    min_new_edge_k.push_back(k);
                                    min_new_edge_l.push_back(l);
                                }
                            }
                        }
                    }
                }
            }
            if(assigned.size() < slots){
                for(int i=0; i<assigned.size(); i++){
                    for(int j=0;j<assigned.size();j++){
                        if(connectable[j]>0 && connections[i][j]!=-1){
                            for(int k=0;k<aa_mass.size();k++){
                                double dist = abs(aa_mass[k] + connections[i][j] - spectrum_sort[sp][0]);
                                //printf("new_aa_edge_dist: %.0lf\tk: %d\ti: %d\tj: %d\tdist: %.0lf\n",new_aa_edge_dist,k,i,j,dist);
                                if(new_aa_edge_dist >= (dist - equal_precision) || new_aa_edge_dist == -1){
                                    if(new_aa_edge_dist > (dist + equal_precision) || new_aa_edge_dist == -1){
                                        new_aa_edge_dist = dist;                                
                                        new_aa_edge_a.clear();
                                        new_aa_edge_i.clear();
                                        new_aa_edge_j.clear();
                                    }
                                    new_aa_edge_a.push_back(k);
                                    new_aa_edge_i.push_back(i);
                                    new_aa_edge_j.push_back(j); 
                                }
                            }
                        }
                    }
                }
            }
        }
    }//Assigned if
    
    //printf("SP: %d\tRange: (%0.lf,%d)\tSlot: (%.0lf,%d)\tEdge: (%.0lf,%d)\tAA Edge: (%0.lf,%d)\tComplete: %d\tEdges: %d\tSize: %.0lf\n",sp,min_range_dist,min_range_j.size(),min_slot_dist,min_slot_j.size(),min_new_edge_dist,min_new_edge_i.size(),new_aa_edge_dist,new_aa_edge_a.size(),complete,edges,spectrum_sort[sp][0]);
    for(int i=0; i<assigned.size(); i++){
        //printf("%d:(%d,%d)\t",i,connectable[i],assignment[i]);
    }
    /*for(int i=0; i<assigned.size(); i++){
        for(int j=0; j<assigned.size(); j++){
            printf("%.0lf (%d)\t",connections[i][j],direct_edges[i][j]);
        }
        printf("\n");
    }*/
    //printf("\n");
    if(min_range_dist != -1 && 
        (min_range_dist <= (min_slot_dist + equal_precision) || min_slot_dist == -1) &&
        (min_range_dist <= (min_new_edge_dist + equal_precision)  || min_new_edge_dist == -1) &&
        (min_range_dist <= (new_aa_edge_dist + equal_precision)  || new_aa_edge_dist == -1)){
        //for(int s=0;s<min_range_j.size();s++){
        for(int s=0;s<=0;s++){
            assert(min_range_i[s] != -1);
            
            vector<int> pass_lassign = my_lassign;
            vector<int> pass_uassign = my_uassign;
            pass_lassign[p] = min_range_i[s]; 
            pass_uassign[p] = min_range_j[s];
        
            graph_greedy_warm_start_helper(edges,complete,sp+1,spectrum_sort,assigned,connectable,connections,direct_edges,assignment,total_error+min_range_dist,my_xassign,pass_lassign,pass_uassign);
        }
    }
    if(min_new_edge_dist != -1 &&
        (min_new_edge_dist <= (min_range_dist + equal_precision)  || min_range_dist == -1) &&
        (min_new_edge_dist <= (min_slot_dist + equal_precision)  || min_slot_dist == -1) &&
        (min_new_edge_dist <= (new_aa_edge_dist + equal_precision)  || new_aa_edge_dist == -1)){
        
        for(int s=0;s<min_new_edge_i.size();s++){
            assert(min_new_edge_i[s] != -1);
            assert(min_new_edge_j[s] != -1);
            assert(min_new_edge_k[s] != -1);
            assert(min_new_edge_l[s] != -1);
            
            
            vector<int> pass_lassign = my_lassign;
            vector<int> pass_uassign = my_uassign;
            vector<int> pass_connectable = connectable;
            vector< vector <double> > pass_connections = connections;
            vector< vector<int> > pass_direct_edges = direct_edges;
        
            pass_uassign[p] = min_new_edge_i[s]; 
            pass_lassign[p] = min_new_edge_j[s];
        
            pass_connectable[min_new_edge_k[s]]--;
            pass_connectable[min_new_edge_l[s]]--;
            pass_direct_edges[min_new_edge_k[s]][min_new_edge_l[s]] = 1;
            pass_direct_edges[min_new_edge_l[s]][min_new_edge_k[s]] = 1;
            //printf("New Edge\ti:%d\tj:%d\tk:%d\tl:%d\n",min_new_edge_i,min_new_edge_j,min_new_edge_k,min_new_edge_l);
            for(int i=0; i<assigned.size(); i++){
                for(int j=0;j<assigned.size();j++){
                    if(i!=j && pass_connections[i][min_new_edge_k[s]] != -1 && pass_connections[min_new_edge_l[s]][j] != -1 && pass_connections[i][j]==-1){
                        pass_connections[i][j] = pass_connections[i][min_new_edge_k[s]] + pass_connections[min_new_edge_l[s]][j];
                        //printf("Assign (A): (%d,%d)[%.0lf] from (%d,%d)[%.0lf] + (%d,%d)[%.0lf]\n",i,j,pass_connections[i][j],i,min_new_edge_k,pass_connections[i][min_new_edge_k],min_new_edge_l,j,pass_connections[min_new_edge_l][j]);
                    }//else printf("Not Assign: (%d,%d)[%.0lf] from (%d,%d)[%.0lf] + (%d,%d)[%.0lf]\n",i,j,pass_connections[i][j],i,min_new_edge_k,pass_connections[i][min_new_edge_k],min_new_edge_l,j,pass_connections[min_new_edge_l][j]);
                    if(i!=j && pass_connections[j][min_new_edge_k[s]] != -1 && pass_connections[min_new_edge_l[s]][i] != -1 && pass_connections[j][i]==-1){
                        pass_connections[j][i] = pass_connections[j][min_new_edge_k[s]] + pass_connections[min_new_edge_l[s]][i];
                        //printf("Assign (B): (%d,%d)[%.0lf] from (%d,%d)[%.0lf] + (%d,%d)[%.0lf]\n",j,i,pass_connections[j][i],j,min_new_edge_k,pass_connections[j][min_new_edge_k],min_new_edge_l,i,pass_connections[min_new_edge_l][i]);
                    }//else printf("Not Assign: (%d,%d)[%.0lf] from (%d,%d)[%.0lf] + (%d,%d)[%.0lf]\n",j,i,pass_connections[j][i],j,min_new_edge_k,pass_connections[j][min_new_edge_k],min_new_edge_l,i,pass_connections[min_new_edge_l][i]);
                
                    if(i!=j && pass_connections[i][min_new_edge_l[s]] != -1 && pass_connections[min_new_edge_k[s]][j] != -1 && pass_connections[i][j]==-1){
                        pass_connections[i][j] = pass_connections[i][min_new_edge_l[s]] + pass_connections[min_new_edge_k[s]][j];
                        //printf("Assign (C): (%d,%d)[%.0lf] from (%d,%d)[%.0lf] + (%d,%d)[%.0lf]\n",i,j,pass_connections[i][j],i,min_new_edge_l,pass_connections[i][min_new_edge_l],min_new_edge_k,j,pass_connections[min_new_edge_k][j]);
                    
                    }//else printf("Not Assign: (%d,%d)[%.0lf] from (%d,%d)[%.0lf] + (%d,%d)[%.0lf]\n",i,j,pass_connections[i][j],i,min_new_edge_l,pass_connections[i][min_new_edge_l],min_new_edge_k,j,pass_connections[min_new_edge_k][j]);
                    if(i!=j && pass_connections[j][min_new_edge_l[s]] != -1 && pass_connections[min_new_edge_k[s]][i] != -1 && pass_connections[j][i]==-1){
                        pass_connections[j][i] = pass_connections[j][min_new_edge_l[s]] + pass_connections[min_new_edge_k[s]][i];
                        //printf("Assign (D): (%d,%d)[%.0lf] from (%d,%d)[%.0lf] + (%d,%d)[%.0lf]\n",j,i,pass_connections[j][i],j,min_new_edge_l,pass_connections[j][min_new_edge_l],min_new_edge_k,i,pass_connections[min_new_edge_k][i]);
                    }//else printf("Not Assign: (%d,%d)[%.0lf] from (%d,%d)[%.0lf] + (%d,%d)[%.0lf]\n",j,i,pass_connections[j][i],j,min_new_edge_l,pass_connections[j][min_new_edge_l],min_new_edge_k,i,pass_connections[min_new_edge_k][i]);
                
                }
            }
        
            graph_greedy_warm_start_helper(edges+1,complete,sp+1,spectrum_sort,assigned,pass_connectable,pass_connections,pass_direct_edges,assignment,total_error+min_new_edge_dist,my_xassign,pass_lassign,pass_uassign);
        }
    }
    if(min_slot_dist != -1 && 
        (min_slot_dist <= (min_range_dist + equal_precision)  || min_range_dist == -1) &&
        (min_slot_dist <= (min_new_edge_dist + equal_precision)  || min_new_edge_dist == -1) &&
        (min_slot_dist <= (new_aa_edge_dist + equal_precision)  || new_aa_edge_dist == -1)){
        for(int s=0;s<min_slot_j.size();s++){
            
            assert(min_slot_j[s] != -1);
            
            vector<int> pass_assigned = assigned;
            vector<int> pass_assignment = assignment;
            vector< vector <double> > pass_connections = connections;
            vector<int> pass_uassign = my_uassign;
            vector<int> pass_lassign = my_lassign;
            vector<int> pass_xassign = my_xassign;
        
            pass_assigned.push_back(1);
            pass_assignment.push_back(min_slot_j[s]);
            
            int k = pass_assigned.size()-1;
            pass_connections[k][k] = aa_mass[min_slot_j[s]];
            pass_uassign[p] = k; 
            pass_lassign[p] = k;
            pass_xassign[k] = min_slot_j[s];
            printf("pass_lassign: ");
            for(int i=0;i<pass_lassign.size();i++){
                printf("%d\t",pass_lassign[i]);
            }
            printf("\n");
    
            graph_greedy_warm_start_helper(edges,complete,sp+1,spectrum_sort,pass_assigned,connectable,pass_connections,direct_edges,pass_assignment,total_error+min_slot_dist,pass_xassign,pass_lassign,pass_uassign);
        }
    }
    if(new_aa_edge_dist != -1 && 
        (new_aa_edge_dist <= (min_slot_dist + equal_precision)  || min_slot_dist == -1) &&
        (new_aa_edge_dist <= (min_new_edge_dist + equal_precision)  || min_new_edge_dist == -1) &&
        (new_aa_edge_dist <= (min_range_dist + equal_precision)  || min_range_dist == -1)){
        
        for(int s=0;s<new_aa_edge_a.size();s++){
            assert(new_aa_edge_i[s] != -1);
            assert(new_aa_edge_j[s] != -1);
            assert(new_aa_edge_a[s] != -1);
        
            vector<int> pass_assigned = assigned;
            vector<int> pass_assignment = assignment;
            vector< vector <double> > pass_connections = connections;
            vector<int> pass_connectable = connectable;
            vector<int> pass_uassign = my_uassign;
            vector<int> pass_lassign = my_lassign;
            vector<int> pass_xassign = my_xassign;
            vector< vector<int> > pass_direct_edges = direct_edges;
        
            assert(new_aa_edge_dist != -1);
            pass_assigned.push_back(1);
            pass_assignment.push_back(new_aa_edge_a[s]);
            int k = pass_assigned.size()-1;
            pass_connections[k][k] = aa_mass[new_aa_edge_a[s]];
        
            pass_connectable[k]--;
            pass_connectable[new_aa_edge_j[s]]--;
            pass_direct_edges[k][new_aa_edge_j[s]] = 1;
            pass_direct_edges[new_aa_edge_j[s]][k] = 1;
        
            pass_lassign[p] = new_aa_edge_i[s]; 
            pass_uassign[p] = k;
        
            pass_xassign[k] = new_aa_edge_a[s];
        
            for(int i=0; i<pass_assigned.size(); i++){
                if(i!=k && pass_connections[i][new_aa_edge_j[s]]!=-1){
                    pass_connections[i][k] = pass_connections[i][new_aa_edge_j[s]] + pass_connections[k][k];
                    pass_connections[k][i] = pass_connections[i][new_aa_edge_j[s]] + pass_connections[k][k];
                }
            }
            graph_greedy_warm_start_helper(edges+1,complete,sp+1,spectrum_sort,pass_assigned,pass_connectable,pass_connections,pass_direct_edges,pass_assignment,total_error+new_aa_edge_dist,pass_xassign,pass_lassign,pass_uassign);
        }
    }
    
    
    
    return total_error;
}

int ilpsolver::graph_greedy_warm_start(){
    vector<int> assigned;
    vector<int> connectable;
    vector< vector <double> > connections;
    vector< vector<int> > direct_edges;
    vector<int> assignment;

    connections.clear();
    connections.resize(slots);
    connectable.resize(slots);
    direct_edges.clear();
    direct_edges.resize(slots);
    for(int i = 0; i < slots; i++){
        connectable[i] = 2;
        connections[i].resize(slots);
        direct_edges[i].resize(slots);
        for(int j=0;j<slots;j++){
            connections[i][j] = -1;
            direct_edges[i][j] = 0;
        }
    }

    xassign.clear();
	for(int i = 0; i < slots; i++)
	{
        xassign.push_back(-1);
	}
	
	
	vector< vector<double> > spectrum_sort;
	spectrum_sort.clear();
	spectrum_sort.resize(spectrum.size());
	lassign.clear();
	uassign.clear();
	for(int p=0;p<spectrum.size();p++){
	    spectrum_sort[p].push_back(spectrum[p]);
	    spectrum_sort[p].push_back((double)p);
	
	    lassign.push_back(-1);
	    uassign.push_back(-1);
	}
	
	std::sort (spectrum_sort.begin(), spectrum_sort.end(), less_than);
	int num_assigned = 0;
    best_total_error = (max_error_allowed>=0)?(max_error_allowed+1):max_error_allowed;
    greedy_search_cycles.clear();
    greedy_search_cycles_error.clear();
    graph_greedy_warm_start_helper(0, false, 0, spectrum_sort, assigned, connectable, connections, direct_edges, assignment, 0.0, xassign, lassign, uassign);


    
	int min_x_index = aa_mass.size();
	int min_x_slot = -1;
	for(int i=0;i<xassign.size();i++){
	    if(xassign[i] < min_x_index){
	        min_x_index = xassign[i];
	        min_x_slot = i;
	    }
	}
	assert(min_x_index!=-1);
	vector<int> xassign_save ;
	xassign_save = xassign;
	for(int i=0;i<xassign.size();i++){
	    printf("Shifting %d->%d (%d)\n",xassign[i],xassign_save[(i+min_x_slot)%slots],i);
	    xassign[i] = xassign_save[(i+min_x_slot)%slots];
	}
	for(int p=0;p<spectrum.size();p++){
	    lassign[p] = (lassign[p]+min_x_slot)%slots;
        uassign[p] = (uassign[p]+min_x_slot)%slots;
	}
	//uassign.clear();
	//lassign.clear();
	return best_total_error;
	
}
