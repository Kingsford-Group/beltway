#include "ilpsolver.h"
#include "config.h"
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
	try
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

		model->getEnv().set(GRB_DoubleParam_TimeLimit, ilp_time_limit);

		model->update();
		//model->write("temp.lp");

		model->optimize();

		collect_results();
	}
	catch(GRBException &e)
	{
		printf("%s\n", e.getMessage().c_str());
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
	return 0;
}

int ilpsolver::add_ordering_cutting_planes()
{
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
								GRBLinExpr expr2 = (lvars[p1][k1] + uvars[p1][l1] + lvars[p2][k2] + uvars[p2][l2] - 4) * bound;
								model->addConstr(expr1, GRB_GREATER_EQUAL, expr2);
							}
							else if(spectrum[p2] > spectrum[p1])
							{
								double bound = 0;
								if(k1 <= l1) bound = (l1 - k1 + 1) * max_weight;
								else bound = (l1 + 1 + slots - k1) * max_weight;
								GRBLinExpr expr1 = rvars[k2][l2] - rvars[k1][l1];
								GRBLinExpr expr2 = (lvars[p1][k1] + uvars[p1][l1] + lvars[p2][k2] + uvars[p2][l2] - 4) * bound;
								model->addConstr(expr1, GRB_GREATER_EQUAL, expr2);
							}
						}
					}
				}
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
		double w = wassign[i];
		double e = eassign[i];
		printf("spectrum %d with mass %.3lf is assigned to interval [%d, %d], with mass %.3lf and error %.3lf (actual %.3lf)\n", 
				i, spectrum[i], l, u, w, e, (w-spectrum[i]));
		
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
