#include "ilpsolver.h"
#include <iostream>
#include <fstream>
#include <sstream>

ilpsolver::ilpsolver(const string &alphabet_file, const string &spectrum_file)
{
	read_alphabet(alphabet_file);
	read_spectrum(spectrum_file);
	compute_upper_bound();

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
	add_amino_acid_variables();
	add_lower_endpoints_variables();
	add_upper_endpoints_variables();
	add_range_variables();
	add_error_variables();

	add_amino_acid_constraints();
	add_lower_endpoints_constraints();
	add_upper_endpoints_constraints();
	add_range_constraints();
	add_error_constraints();

	set_objective();

	model->update();
	model->optimize();

	collect_results();
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
	while(getline(fin, line))
	{
		if(line == "") continue;
		stringstream sstr(line);
		sstr >> slots;
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
	return 0;
}

int ilpsolver::compute_upper_bound()
{
	ubound = 0;
	for(int i = 0; i < spectrum.size(); i++)
	{
		if(spectrum[i] > ubound) ubound = spectrum[i];
	}
	ubound = ubound * slots;
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
			GRBVar var = model->addVar(0, 1, 0, GRB_BINARY);
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
			GRBVar var = model->addVar(0, 1, 0, GRB_BINARY);
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
			GRBVar var = model->addVar(0, 1, 0, GRB_BINARY);
			uvars[i].push_back(var);
		}
	}
	model->update();
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
			GRBVar var = model->addVar(0, ubound, 0, GRB_CONTINUOUS);
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
		GRBVar var = model->addVar(0, ubound, 1, GRB_CONTINUOUS);
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
	for(int p = 0; p < spectrum.size(); p++)
	{
		for(int k = 0; k < slots; k++)
		{
			for(int l = 0; l < slots; l++)
			{
				GRBLinExpr expr1 = rvars[k][l] - spectrum[p] + ubound * (lvars[p][k] + uvars[p][l] - 2);
				GRBLinExpr expr2 = spectrum[p] - rvars[k][l] + ubound * (lvars[p][k] + uvars[p][l] - 2);
				model->addConstr(evars[p], GRB_GREATER_EQUAL, expr1);
				model->addConstr(evars[p], GRB_GREATER_EQUAL, expr2);
			}
		}
	}
	return 0;
}

int ilpsolver::set_objective()
{
	GRBLinExpr expr;
	for(int i = 0; i < evars.size(); i++) expr += evars.at(i);
	model->setObjective(expr, GRB_MINIMIZE);
	return 0;
}

int ilpsolver::collect_results()
{
	xassign.clear();
	for(int i = 0; i < slots; i++)
	{
		int k = -1;
		for(int j = 0; j < aa_list.size(); j++)
		{
			if(xvars[i][j].get(GRB_DoubleAttr_X) <= 0.5) continue;
			assert(k == -1);
			k = j;
		}
		assert(k >= 0);
		xassign.push_back(k);
	}

	lassign.clear();
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

	uassign.clear();
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

	eassign.clear();
	for(int p = 0; p < spectrum.size(); p++)
	{
		double e = evars[p].get(GRB_DoubleAttr_X);
		eassign.push_back(e);
	}

	return 0;
}

int ilpsolver::print()
{
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
	printf("upper bound = %.3lf\n", ubound);

	for(int i = 0; i < xassign.size(); i++)
	{
		int k = xassign[i];
		printf("slot %d is assigned amino acid %d: %s -> %.3lf\n", i, k, aa_list[k].c_str(), aa_mass[k]);
	}
	for(int i = 0; i < spectrum.size(); i++)
	{
		int l = lassign[i];
		int u = uassign[i];
		double w = wassign[i];
		double e = eassign[i];
		printf("spectrum %d with mass %.3lf is assigned to interval [%d, %d], with mass %.3lf and error %.3lf\n", i, spectrum[i], l, u, w, e);
	}
	return 0;
}
